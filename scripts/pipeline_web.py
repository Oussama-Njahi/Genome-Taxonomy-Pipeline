#!/usr/bin/env python3
# pipeline_web.py - Parameterized version of pipeline.py for the web interface.
#
# Same 9 steps, same tools, same logic as pipeline.py. The only difference is
# that genome folder, protein folder and results folder are passed as
# command-line arguments instead of being hardcoded, so several independent
# jobs (one per web upload) can run without touching each other's files.
#
# Usage:
#   python pipeline_web.py --genomes-dir <dir> --output-dir <dir> [--job-id <id>] [--threads 4] [--steps qc,ani,aai]
#                          [--outgroup <genome>]
#
# Progress is reported on stdout as lines "PROGRESS <n>/9 <label>" so a caller
# (e.g. the FastAPI backend) can parse them to drive a progress bar.

import argparse
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
import uuid

STEPS = [
    "Quality control",
    "CheckM2 (completeness/contamination)",
    "Protein prediction (Prodigal)",
    "ANI (FastANI)",
    "dDDH (BLASTN, GBDP formula 2)",
    "AAI (EzAAI)",
    "POCP (DIAMOND)",
    "Phylogenomic tree (EasyCGTree)",
    "Figures (R)",
]

# Cles utilisees avec --steps, dans le meme ordre que STEPS ci-dessus
STEPS_KEYS = ["qc", "checkm2", "proteines", "ani", "dddh", "aai", "pocp", "arbre", "figures"]

# Genomes en entree : extensions acceptees (sans tenir compte de la casse), nombre minimum (l'arbre d'EasyCGTree
# en demande 5, les comparaisons par paires 2), caracteres permis dans un nom (un nom sert d'etiquette dans l'arbre
# Newick, ou parentheses, virgules, deux-points et espaces ont un sens), lettres d'ADN permises (code IUPAC).
EXTENSIONS = (".fas", ".fasta", ".fna")
MIN_GENOMES = 2
MIN_GENOMES_ARBRE = 5
NOM_PERMIS = re.compile(r"^[A-Za-z0-9._-]+$")
IUPAC = set("ACGTURYSWKMBDHVN-")
TAILLE_SUSPECTE = 500000        # en dessous, prevenir : un genome bacterien fait rarement moins de 0.5 Mb

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPTS)
EASYCG = BASE + "/EasyCGTree4/EasyCGTree.v4.2-Linux"


# Etapes choisies avec --steps, remplie par main() : la progression est
# numerotee parmi elles (un run "ani,dddh" affiche 1/2 puis 2/2, pas 4/9)
SELECTION = list(STEPS_KEYS)


def announce(step_index, extra=""):
    label = STEPS[step_index - 1]
    key = STEPS_KEYS[step_index - 1]
    position = SELECTION.index(key) + 1 if key in SELECTION else step_index
    print("PROGRESS %d/%d %s%s" % (position, len(SELECTION), label, (" - " + extra) if extra else ""))
    sys.stdout.flush()


def fail(message):
    print("STEP FAILED: %s" % message)
    sys.stdout.flush()
    sys.exit(1)


def run(cmd, **kwargs):
    """Run a subprocess and stop the pipeline with a clear message if it fails."""
    print("    $ " + " ".join(str(c) for c in cmd))
    sys.stdout.flush()
    result = subprocess.run(cmd, **kwargs)
    if result.returncode != 0:
        fail("command failed (exit code %d): %s" % (result.returncode, " ".join(str(c) for c in cmd)))


def etapes_demandees(steps_arg):
    """Traduit la valeur de --steps (texte) en liste de noms d'etapes valides."""
    if steps_arg == "all":
        return list(STEPS_KEYS)
    choisies = [s.strip() for s in steps_arg.split(",") if s.strip()]
    for s in choisies:
        if s not in STEPS_KEYS:
            fail("unknown step '%s' (valid steps: %s)" % (s, ", ".join(STEPS_KEYS)))
    return choisies


def lire_fasta(chemin):
    """Verifie qu'un fichier est bien un FASTA d'ADN et renvoie ses statistiques :
    (longueur totale, nombre de G+C, nombre de A/C/G/T, nombre de sequences, empreinte SHA-1)."""
    nom = os.path.basename(chemin)
    empreinte = hashlib.sha1()
    total = gc = acgt = sequences = 0
    vu_entete = False
    with open(chemin, "rb") as fh:
        for numero, brut in enumerate(fh, 1):
            empreinte.update(brut)
            ligne = brut.decode("ascii", "replace").strip()
            if not ligne:
                continue
            if ligne.startswith(">"):
                vu_entete = True
                sequences += 1
                continue
            if not vu_entete:
                fail("%s is not a FASTA file: its first line should start with '>'" % nom)
            ligne = ligne.upper()
            inconnues = set(ligne) - IUPAC
            if inconnues:
                fail("%s, line %d: characters that are not DNA (%s) - is it really a genome assembly?"
                     % (nom, numero, "".join(sorted(inconnues))[:10]))
            total += len(ligne)
            gc += ligne.count("G") + ligne.count("C")
            acgt += sum(ligne.count(b) for b in "ACGT")
    if total == 0:
        fail("%s contains no sequence (empty file or headers only)" % nom)
    return total, gc, acgt, sequences, empreinte.hexdigest()


def lister_genomes(genomes_dir, selection):
    """Toujours necessaire : trouve les genomes (extension sans tenir compte de la casse), verifie que chacun est un
    vrai FASTA d'ADN au nom utilisable, sans doublon, et qu'il y en a assez. Renvoie {nom: (chemin, stats)}."""
    fichiers = sorted(f for f in glob.glob(genomes_dir + "/*") if os.path.splitext(f)[1].lower() in EXTENSIONS)
    if not fichiers:
        fail("no .fas/.fasta/.fna genome files found in %s" % genomes_dir)
    minimum = MIN_GENOMES_ARBRE if "arbre" in selection else MIN_GENOMES
    if len(fichiers) < minimum:
        if minimum == MIN_GENOMES_ARBRE:
            fail("the phylogenomic tree needs at least %d genomes (EasyCGTree cannot build one with fewer); got %d. "
                 "Add genomes, or leave out the 'arbre' step." % (minimum, len(fichiers)))
        fail("at least %d genomes are needed to compare them; got %d" % (minimum, len(fichiers)))

    genomes = {}
    for chemin in fichiers:
        nom = os.path.splitext(os.path.basename(chemin))[0]
        if not NOM_PERMIS.match(nom):
            fail("'%s': genome names may only contain letters, digits, '.', '_' and '-' (the name is used as a "
                 "label in the tree file). Rename the file, e.g. '%s'"
                 % (os.path.basename(chemin), re.sub(r"[^A-Za-z0-9._-]+", "_", nom).strip("._") or "genome"))
        if nom.lower() in (n.lower() for n in genomes):
            fail("two files give the same genome name '%s' (%s and %s): rename one of them"
                 % (nom, os.path.basename(genomes[next(n for n in genomes if n.lower() == nom.lower())][0]),
                    os.path.basename(chemin)))
        genomes[nom] = (chemin, lire_fasta(chemin))
        if genomes[nom][1][0] < TAILLE_SUSPECTE:
            print("    !! WARNING: %s is only %d bp - is it a complete bacterial genome?"
                  % (os.path.basename(chemin), genomes[nom][1][0]))
    print("    %d genomes checked" % len(genomes))
    return genomes


def preparer_dossiers(genomes, output_dir, results_dir, proteins_dir):
    """Les resultats d'une etape restent d'un run a l'autre (on peut lancer 'qc' puis 'checkm2' plus tard), mais
    seulement pour LES MEMES genomes : sinon les figures melangeraient d'anciens resultats avec les nouveaux.
    results/genomes.tsv note les genomes analyses ; s'ils ont change, les anciens results/ et proteins/ sont
    deplaces dans previous_runs/<date>/. Les genomes sont ensuite exposes dans results/input_genomes/ sous des noms
    uniformes (<nom>.fna), seule liste lue par les etapes."""
    manifeste = results_dir + "/genomes.tsv"
    actuel = "".join("%s\t%d\t%s\n" % (nom, stats[0], stats[4]) for nom, (_, stats) in sorted(genomes.items()))
    ancien = open(manifeste).read() if os.path.exists(manifeste) else None
    reste = [d for d in (results_dir, proteins_dir) if os.path.isdir(d) and os.listdir(d)]
    if reste and ancien != actuel:
        archive = os.path.join(output_dir, "previous_runs", time.strftime("%Y-%m-%d_%H%M%S"))
        os.makedirs(archive, exist_ok=True)
        for d in reste:
            shutil.move(d, os.path.join(archive, os.path.basename(d)))
        print("    !! The genomes differ from the previous run: its results were moved to %s" % archive)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(proteins_dir, exist_ok=True)
    with open(manifeste, "w") as fh:
        fh.write(actuel)

    entree = results_dir + "/input_genomes"
    shutil.rmtree(entree, ignore_errors=True)
    os.makedirs(entree)
    fichiers = []
    for nom, (chemin, _) in sorted(genomes.items()):
        lien = os.path.join(entree, nom + ".fna")
        os.symlink(os.path.abspath(chemin), lien)
        fichiers.append(lien)
    return entree, fichiers


def lire_checkm2(results_dir, fichiers, threads):
    print("    (running CheckM2, this may take a few minutes...)")

    groupes = {}
    for chemin in fichiers:
        ext = os.path.splitext(chemin)[1].lstrip(".")
        groupes.setdefault(ext, []).append(chemin)

    qualite = {}
    for ext, genomes_du_groupe in groupes.items():
        sous_dossier = results_dir + "/checkm2_input_" + ext
        os.makedirs(sous_dossier, exist_ok=True)
        for chemin in genomes_du_groupe:
            shutil.copy(chemin, sous_dossier)

        sortie = results_dir + "/checkm2_out_" + ext
        run(["conda", "run", "-n", "checkm2env",
             "checkm2", "predict",
             "--input", sous_dossier,
             "--output-directory", sortie,
             "-x", ext, "--threads", str(threads), "--force"])

        rapport = sortie + "/quality_report.tsv"
        with open(rapport) as f:
            next(f)
            for ligne in f:
                cols = ligne.strip().split("\t")
                nom, completude, contamination = cols[0], cols[1], cols[2]
                qualite[nom] = (completude, contamination)

    return qualite


def etape_qc(genomes, results_dir):
    """Stats de base uniquement (taille, GC%, contigs) - rapide, pas de CheckM2. Les statistiques viennent de la
    verification des genomes (lire_fasta). Le %GC est calcule sur les seules bases A/C/G/T, pour que les N d'un
    assemblage ne le fassent pas baisser. Completeness/Contamination sont ecrits a "NA" ; lancer aussi l'etape
    'checkm2' pour les remplir (voir etape_checkm2 ci-dessous)."""
    announce(1)
    with open(results_dir + "/qc.txt", "w") as rapport:
        rapport.write("Genome\tTaille(pb)\tGC(%)\tContigs\tCompleteness(%)\tContamination(%)\n")
        for nom, (_, (total, gc, acgt, sequences, _)) in sorted(genomes.items()):
            taux = "%.1f" % (100.0 * gc / acgt) if acgt else "NA"
            rapport.write("%s\t%d\t%s\t%d\tNA\tNA\n" % (nom, total, taux, sequences))
    print("    -> qc.txt created (Completeness/Contamination = NA; also run 'checkm2' step to fill them in)")


def etape_checkm2(results_dir, fichiers, threads):
    """Lance CheckM2 et met a jour les colonnes Completeness/Contamination
    de qc.txt (qui doit deja exister, ecrit par l'etape 'qc' - verifie par
    l'appelant dans main(), comme pour la dependance proteines -> pocp)."""
    announce(2)
    qc_path = results_dir + "/qc.txt"
    qualite = lire_checkm2(results_dir, fichiers, threads)
    with open(qc_path) as f:
        lignes = f.readlines()
    nouvelles_lignes = [lignes[0]]
    for ligne in lignes[1:]:
        cols = ligne.rstrip("\n").split("\t")
        comp, cont = qualite.get(cols[0], ("NA", "NA"))
        cols[4], cols[5] = comp, cont
        nouvelles_lignes.append("\t".join(cols) + "\n")
    with open(qc_path, "w") as f:
        f.writelines(nouvelles_lignes)
    print("    -> qc.txt updated with Completeness and Contamination")


def etape_proteines(fichiers, proteins_dir):
    announce(3)
    os.makedirs(proteins_dir, exist_ok=True)
    for chemin in fichiers:
        nom = os.path.splitext(os.path.basename(chemin))[0]
        run(["prodigal", "-i", chemin, "-a", proteins_dir + "/" + nom + ".faa", "-o", os.devnull, "-q"])
    print("    -> proteins written to %s" % proteins_dir)


def etape_ani(fichiers, results_dir, threads):
    announce(4)
    liste = results_dir + "/liste.txt"
    with open(liste, "w") as f:
        for g in fichiers:
            f.write(g + "\n")
    run(["fastANI", "--ql", liste, "--rl", liste,
         "-o", results_dir + "/ani_resultats.txt", "-t", str(threads)])
    print("    -> ani_resultats.txt created")


def etape_dddh(genomes_dir, results_dir, threads):
    announce(5)
    run(["python", SCRIPTS + "/dddh.py", genomes_dir, results_dir + "/dddh_out", str(threads)])
    print("    -> dddh_out/dddh_distance.tsv and dddh_similarite.tsv created")


def etape_aai(fichiers, results_dir, threads):
    announce(6)
    aai_db = results_dir + "/aai_db"
    os.makedirs(aai_db, exist_ok=True)
    for chemin in fichiers:
        nom = os.path.splitext(os.path.basename(chemin))[0]
        run(["EzAAI", "extract", "-i", chemin, "-o", aai_db + "/" + nom + ".db", "-l", nom])
    run(["EzAAI", "calculate", "-i", aai_db, "-j", aai_db,
         "-o", results_dir + "/aai_resultats.tsv", "-t", str(threads)])
    print("    -> aai_resultats.tsv created")


def etape_pocp(proteins_dir, results_dir, threads):
    announce(7)
    run(["python", SCRIPTS + "/pocp.py", proteins_dir, results_dir + "/pocp_out", str(threads)])
    print("    -> pocp_out/pocp_matrice.tsv created")


def nom_easycgtree(nom):
    """EasyCGTree renomme ses fichiers d'entree avant de les lire (sub rename de EasyCGTree.pl) : espaces et '_'
    deviennent '-', et ce nom-la sert d'etiquette dans l'arbre. Les autres caracteres permis (NOM_PERMIS) restent."""
    return re.sub(r"[\s_]", "-", nom)


def etape_arbre(fichiers, results_dir, job_id, threads, outgroup=None):
    announce(8)
    noms = [os.path.splitext(os.path.basename(f))[0] for f in fichiers]
    vrai_nom = {}
    for nom in noms:
        autre = vrai_nom.setdefault(nom_easycgtree(nom), nom)
        if autre != nom:
            fail("EasyCGTree turns '_' into '-', so '%s' and '%s' would get the same name in the tree: "
                 "rename one of them" % (autre, nom))
    if not os.path.exists(EASYCG + "/EasyCGTree.pl"):
        fail("EasyCGTree4 not found at %s - it must be installed manually "
             "(not available via conda), see the README's Installation section" % EASYCG)
    input_name = "job_" + job_id
    entree = EASYCG + "/" + input_name
    os.makedirs(entree, exist_ok=True)
    for g in fichiers:
        shutil.copy(g, entree)
    run(["perl", "EasyCGTree.pl", "-input", input_name,
         "-hmm", "bac120", "-tree", "sm", "-tree_app", "fasttree",
         "-thread", str(threads)], cwd=EASYCG)
    arbre = EASYCG + "/" + input_name + ".bac120.120.supermatrix.fasttree.tree"
    if not os.path.exists(arbre):
        fail("tree file not found (%s) - EasyCGTree may have failed" % arbre)
    # les feuilles reprennent les noms des genomes (ceux des matrices), pas ceux qu'EasyCGTree leur a donnes
    with open(arbre) as fh:
        newick = fh.read()
    feuilles = re.findall(r"(?<=[(,])([^(),:;]+)(?=[:),;])", newick)
    inconnues = sorted(set(feuilles) - set(vrai_nom))
    if inconnues:
        fail("unexpected leaf names in EasyCGTree's tree: %s" % ", ".join(inconnues))
    newick = re.sub(r"(?<=[(,])([^(),:;]+)(?=[:),;])", lambda m: vrai_nom[m.group(1)], newick)
    # FastTree ecrit un arbre non enracine : on le garde tel quel, et arbre.tree
    # (lu par les figures et l'interface) est sa version enracinee, au point
    # median ou sur le groupe externe donne avec --outgroup.
    brut = results_dir + "/arbre_non_enracine.tree"
    with open(brut, "w") as fh:
        fh.write(newick)
    cleanup_easycgtree(input_name)
    cmd = ["python", SCRIPTS + "/root_tree.py", brut, results_dir + "/arbre.tree"]
    if outgroup:
        cmd += ["--outgroup", outgroup]
    run(cmd)
    print("    -> arbre.tree (rooted) and arbre_non_enracine.tree copied to results")


def cleanup_easycgtree(input_name):
    """Remove this job's working files from the shared EasyCGTree tool directory."""
    for path in glob.glob(EASYCG + "/" + input_name + "*"):
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                os.remove(path)
            except OSError:
                pass


def etape_figures(results_dir):
    announce(9)
    run(["Rscript", SCRIPTS + "/plot_all_web.R", results_dir])
    print("    -> PNG figures created in %s" % results_dir)


def main():
    parser = argparse.ArgumentParser(description="Run the taxonomy pipeline on an arbitrary genome folder.")
    parser.add_argument("--genomes-dir", required=True, help="Folder containing the genome files")
    parser.add_argument("--output-dir", required=True, help="Folder where proteins/ and results/ will be created")
    parser.add_argument("--job-id", default=None, help="Unique id for this run (default: random uuid)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                         help="Threads for the tools that can use several (default: all CPUs, %(default)s here)")
    parser.add_argument("--steps", type=str, default="all",
                         help="Comma-separated steps to run (" + ", ".join(STEPS_KEYS) + "). Default: all")
    parser.add_argument("--outgroup", type=str, default=None,
                         help="Genome name(s) to root the tree on, comma-separated, without the file extension "
                              "(default: midpoint rooting)")
    args = parser.parse_args()

    selection = etapes_demandees(args.steps)
    # dans l'ordre d'execution, pour que la progression suive les etapes lancees
    SELECTION[:] = [k for k in STEPS_KEYS if k in selection]

    job_id = args.job_id or uuid.uuid4().hex[:12]
    genomes_dir = os.path.abspath(args.genomes_dir)
    proteins_dir = os.path.abspath(args.output_dir) + "/proteins"
    results_dir = os.path.abspath(args.output_dir) + "/results"

    print("===== PIPELINE START (job %s) =====" % job_id)
    print("Steps selected: %s" % ", ".join(selection))

    genomes = lister_genomes(genomes_dir, selection)   # toujours necessaire, rapide
    genomes_dir, fichiers = preparer_dossiers(genomes, os.path.abspath(args.output_dir), results_dir, proteins_dir)

    if "qc" in selection:
        etape_qc(genomes, results_dir)
    if "checkm2" in selection:
        if not os.path.exists(results_dir + "/qc.txt"):
            print("    !! WARNING: qc.txt not found in %s" % results_dir)
            print("    !! The 'qc' step must be run at least once before 'checkm2'.")
            print("    !! Re-run with --steps qc,checkm2")
        else:
            etape_checkm2(results_dir, fichiers, args.threads)
    if "proteines" in selection:
        etape_proteines(fichiers, proteins_dir)
    if "ani" in selection:
        etape_ani(fichiers, results_dir, args.threads)
    if "dddh" in selection:
        etape_dddh(genomes_dir, results_dir, args.threads)
    if "aai" in selection:
        etape_aai(fichiers, results_dir, args.threads)
    if "pocp" in selection:
        if not os.path.isdir(proteins_dir) or not glob.glob(proteins_dir + "/*.faa"):
            print("    !! WARNING: no protein files found in %s" % proteins_dir)
            print("    !! The 'proteines' step must be run at least once before POCP.")
            print("    !! Re-run with --steps proteines,pocp")
        else:
            etape_pocp(proteins_dir, results_dir, args.threads)
    if "arbre" in selection:
        etape_arbre(fichiers, results_dir, job_id, args.threads, args.outgroup)
    if "figures" in selection:
        etape_figures(results_dir)

    print("===== PIPELINE DONE =====")


if __name__ == "__main__":
    main()
