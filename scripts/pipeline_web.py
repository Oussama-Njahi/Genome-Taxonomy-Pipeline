#!/usr/bin/env python3
# pipeline_web.py - The pipeline, for the web interface and the command line (pipeline.py calls it).
#
# Genome folder, protein folder and results folder are passed as command-line
# arguments, so several independent jobs (one per web upload) can run without
# touching each other's files. This script checks the genomes, prepares the
# folders and writes versions.txt; the 9 steps themselves are a Snakemake
# workflow (workflow/Snakefile), run here on the selected steps; the pairwise
# summary is written at the end.
#
# Usage:
#   python pipeline_web.py --genomes-dir <dir> --output-dir <dir> [--job-id <id>] [--threads 4] [--steps qc,ani,aai]
#                          [--outgroup <genome>]
#
# Progress is reported on stdout so a caller (e.g. the FastAPI backend) can
# drive a progress bar: first "STEPS <key>,<key>,..." (the steps that will
# run), then "PROGRESS start <key> <label>" and "PROGRESS done <key> <label>"
# for each of them. Several steps can be in progress at the same time.

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import signal
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
SNAKEFILE = BASE + "/workflow/Snakefile"

# Regles du Snakefile de chaque etape, relancees (--forcerun) quand l'etape est demandee
REGLES = {"qc": ["qc"], "checkm2": ["checkm2"], "proteines": ["prodigal"], "ani": ["ani"], "dddh": ["dddh"],
          "aai": ["ezaai_extract", "aai"], "pocp": ["pocp"], "arbre": ["arbre"], "figures": ["figures"]}

# Lignes de Snakemake sans interet dans le journal (une paire par tache lancee)
BRUIT = ("Assuming unrestricted shared filesystem usage", "Using shell:", "Select jobs to execute...",
         "Execute ", "Building DAG of jobs...")


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
    (longueur totale, nombre de G+C, nombre de A/C/G/T, nombre de sequences, empreinte SHA-1, N50)."""
    nom = os.path.basename(chemin)
    empreinte = hashlib.sha1()
    total = gc = acgt = sequences = 0
    longueurs = []
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
                longueurs.append(0)
                continue
            if not vu_entete:
                fail("%s is not a FASTA file: its first line should start with '>'" % nom)
            ligne = ligne.upper()
            inconnues = set(ligne) - IUPAC
            if inconnues:
                fail("%s, line %d: characters that are not DNA (%s) - is it really a genome assembly?"
                     % (nom, numero, "".join(sorted(inconnues))[:10]))
            total += len(ligne)
            longueurs[-1] += len(ligne)
            gc += ligne.count("G") + ligne.count("C")
            acgt += sum(ligne.count(b) for b in "ACGT")
    if total == 0:
        fail("%s contains no sequence (empty file or headers only)" % nom)
    # N50 : longueur du contig qui, avec tous les plus longs, couvre la moitie de l'assemblage
    cumul = 0
    for n50 in sorted(longueurs, reverse=True):
        cumul += n50
        if 2 * cumul >= total:
            break
    return total, gc, acgt, sequences, empreinte.hexdigest(), n50


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
        # le lien prend la date du genome : Snakemake compare aussi celle du lien aux resultats, et un lien
        # recree a chaque run ferait refaire toutes les etapes gardees du run d'avant
        date = os.stat(chemin)
        os.utime(lien, ns=(date.st_atime_ns, date.st_mtime_ns), follow_symlinks=False)
        fichiers.append(lien)
    return entree, fichiers


# Commandes qui affichent la version de chaque outil, par etape (une ligne suffit, voir ecrire_versions)
VERSIONS = [
    ("python", None, [sys.executable, "--version"]),
    ("Snakemake", None, [sys.executable, "-m", "snakemake", "--version"]),
    ("CheckM2", "checkm2", ["conda", "run", "-n", "checkm2env", "checkm2", "--version"]),
    ("Prodigal", "proteines", ["prodigal", "-v"]),
    ("FastANI", "ani", ["fastANI", "--version"]),
    ("BLAST+ (blastn)", "dddh", ["blastn", "-version"]),
    ("EzAAI", "aai", ["EzAAI", "-h"]),
    ("DIAMOND", "pocp", ["diamond", "version"]),
    ("EasyCGTree", "arbre", None),
    ("R", "figures", ["Rscript", "--version"]),
]


def ecrire_versions(results_dir, selection, args):
    """Trace du run dans results/versions.txt : date, commande, version du pipeline (commit git) et version de
    chaque outil utilise. Ne fait jamais echouer le run : un outil introuvable est note comme tel."""
    def premiere_ligne(cmd):
        try:
            sortie = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired):
            return "not found"
        texte = re.sub(r"\x1b\[[0-9;]*m", "", sortie.stdout + sortie.stderr)
        lignes = [l.strip() for l in texte.splitlines() if re.search(r"\d+\.\d+", l)]
        return lignes[0] if lignes else "unknown"

    try:
        commit = subprocess.run(["git", "-C", SCRIPTS, "log", "-1", "--format=%h (%cd)", "--date=short"],
                                capture_output=True, text=True, timeout=30).stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        commit = "unknown"
    with open(results_dir + "/versions.txt", "w") as fh:
        fh.write("date\t%s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
        fh.write("command\t%s\n" % " ".join(sys.argv))
        fh.write("steps\t%s\n" % ",".join(selection))
        fh.write("threads\t%d\n" % args.threads)
        fh.write("pipeline commit\t%s\n" % commit)
        for outil, etape, cmd in VERSIONS:
            if etape is not None and etape not in selection:
                continue
            version = os.path.basename(EASYCG) if cmd is None else premiere_ligne(cmd)
            fh.write("%s\t%s\n" % (outil, version))


# Seuils, identiques a ceux de l'interface (webapp/frontend/js/plates.js, METRICS) et de dddh.py
SEUILS = {"ani": 95.0, "dddh": 0.0412, "aai": 65.0, "pocp": 50.0, "fraction": 0.5}


def ecrire_synthese(results_dir):
    """results/synthese_paires.tsv : une ligne par paire de genomes avec toutes les mesures calculees et les
    verdicts. Ecrit a la fin de chaque run, a partir des resultats presents (une colonne reste a NA si son etape
    n'a jamais ete lancee). Verdict d'espece : ANI et dDDH ; verdict de genre : POCP (seuil publie) et AAI
    (seuil indicatif) ; 'conflict' quand les deux criteres d'un rang se contredisent."""
    def matrice(chemin):
        with open(chemin) as fh:
            lignes = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        noms = lignes[0][1:]
        return {(l[0], n): (None if v in ("", "NA") else float(v)) for l in lignes[1:] for n, v in zip(noms, l[1:])}

    nom = lambda chemin: os.path.splitext(os.path.basename(chemin))[0]
    ani, aai, dddh, estim, pocp = {}, {}, {}, {}, {}
    if os.path.exists(results_dir + "/ani_resultats.txt"):
        for ligne in open(results_dir + "/ani_resultats.txt"):
            c = ligne.split("\t")
            if len(c) >= 5:
                ani[(nom(c[0]), nom(c[1]))] = (float(c[2]), float(c[3]) / float(c[4]))
    if os.path.exists(results_dir + "/aai_resultats.tsv"):
        lignes = [l.rstrip("\n").split("\t") for l in open(results_dir + "/aai_resultats.tsv") if l.strip()]
        tete = [h.strip().lower() for h in lignes[0]]
        ia, ib, iv = tete.index("label 1"), tete.index("label 2"), tete.index("aai")
        for c in lignes[1:]:
            aai[(c[ia], c[ib])] = float(c[iv])
    if os.path.exists(results_dir + "/dddh_out/dddh_distance.tsv"):
        dddh = matrice(results_dir + "/dddh_out/dddh_distance.tsv")
        estim = matrice(results_dir + "/dddh_out/dddh_similarite.tsv")
    if os.path.exists(results_dir + "/pocp_out/pocp_matrice.tsv"):
        pocp = matrice(results_dir + "/pocp_out/pocp_matrice.tsv")
    noms = sorted({n for d in (ani, aai, dddh, pocp) for paire in d for n in paire})
    if len(noms) < 2:
        return

    def verdict(a, b):
        if a is None or b is None:
            return a or b or "NA"
        return a if a == b else "conflict"

    lignes = []
    for a, b in ((a, b) for i, a in enumerate(noms) for b in noms[i + 1:]):
        sens = [ani.get(p) for p in ((a, b), (b, a)) if ani.get(p)]
        ani_ab = sum(v for v, _ in sens) / len(sens) if sens else None
        frac = min(f for _, f in sens) if sens else None
        d = dddh.get((a, b))
        aai_ab = aai.get((a, b), aai.get((b, a)))
        pocp_ab = pocp.get((a, b))
        # une paire absente de FastANI (pas d'alignement, ANI < ~80 %) est d'especes differentes
        esp_ani = None if not ani else ("different" if ani_ab is None else
                                        ("same" if ani_ab >= SEUILS["ani"] else "different"))
        esp_dddh = None if d is None else ("same" if d <= SEUILS["dddh"] else "different")
        espece = verdict(esp_ani, esp_dddh)
        if espece == "same" and frac is not None and frac < SEUILS["fraction"]:
            espece = "same (low aligned fraction)"
        gen_pocp = None if pocp_ab is None else ("same" if pocp_ab >= SEUILS["pocp"] else "different")
        gen_aai = None if aai_ab is None else ("same" if aai_ab >= SEUILS["aai"] else "different")
        fmt = lambda v, f: "NA" if v is None else f % v
        lignes.append([a, b, fmt(ani_ab, "%.2f"), fmt(frac, "%.2f"), fmt(d, "%.4f"), fmt(estim.get((a, b)), "%.1f"),
                       fmt(aai_ab, "%.2f"), fmt(pocp_ab, "%.2f"), espece, verdict(gen_pocp, gen_aai)])
    with open(results_dir + "/synthese_paires.tsv", "w") as fh:
        fh.write("\t".join(["genome_a", "genome_b", "ani", "ani_aligned_fraction", "dddh_distance", "dddh_estimate",
                            "aai", "pocp", "same_species", "same_genus"]) + "\n")
        for l in lignes:
            fh.write("\t".join(l) + "\n")
    print("    -> synthese_paires.tsv: %d pairs, every measure and the species / genus calls" % len(lignes))


def etapes_prevues(selection, noms, proteins_dir):
    """Etapes qui vont tourner : celles demandees, plus celles dont elles ont besoin et qui n'ont pas encore leurs
    fichiers. C'est Snakemake qui les ajoute (voir workflow/Snakefile) ; les connaitre ici permet de les annoncer
    des le depart (ligne STEPS) et de les viser, pour que chacune annonce aussi sa fin."""
    prevues = set(selection)
    if "checkm2" in prevues:
        prevues.add("qc")                    # qc.txt est reecrit avec les valeurs de CheckM2
    if "pocp" in prevues and not all(os.path.exists("%s/%s.faa" % (proteins_dir, n)) for n in noms):
        prevues.add("proteines")
    return [k for k in STEPS_KEYS if k in prevues]


def commande_snakemake(prevues, selection, output_dir, configfile, threads):
    """Snakemake vise les etapes prevues et relance celles qui ont ete demandees ; les autres ne tournent que si
    leurs fichiers manquent. Seule la date des fichiers decide si un resultat est perime (--rerun-triggers mtime),
    pas une modification du Snakefile ou de l'environnement : un resultat garde reste celui du run d'avant."""
    return ([sys.executable, "-m", "snakemake"] + prevues +
            ["--snakefile", SNAKEFILE, "--directory", output_dir, "--configfile", configfile,
             "--cores", str(threads), "--forcerun"] + [r for k in selection for r in REGLES[k]] +
            ["--rerun-triggers", "mtime", "--rerun-incomplete", "--nolock", "--sdm", "conda",
             "--quiet", "rules", "host", "progress"])


def lancer_snakemake(cmd, prevues):
    """Lance Snakemake en recopiant sa sortie (sauf les lignes sans interet, BRUIT) et arrete le pipeline s'il
    echoue, en nommant les etapes qui etaient en cours si aucune n'a deja dit pourquoi elle echouait."""
    commencees, finies, explique, annule = set(), set(), False, []
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

    # Annulation (SIGTERM, envoye par l'interface a tout le groupe de processus) : on la transmet a Snakemake et on
    # continue de lire sa sortie jusqu'a ce qu'il ait arrete ses taches. Sans lecteur, Snakemake resterait bloque.
    def arreter(signum, frame):
        annule.append(signum)
        proc.terminate()
    signal.signal(signal.SIGTERM, arreter)

    for ligne in proc.stdout:
        if ligne.startswith(BRUIT):
            continue
        print(ligne, end="")
        sys.stdout.flush()
        mots = ligne.split()
        if mots[:1] == ["PROGRESS"] and len(mots) >= 3:
            (commencees if mots[1] == "start" else finies).add(mots[2])
        elif ligne.startswith("STEP FAILED"):
            explique = True
    proc.wait()
    if annule:
        print("===== PIPELINE CANCELLED =====")
        sys.exit(1)
    if proc.returncode != 0:
        if explique:
            sys.exit(1)
        en_cours = [STEPS[STEPS_KEYS.index(k)] for k in prevues if k in commencees and k not in finies]
        fail("%s failed, see the log above" % (" / ".join(en_cours) if en_cours else "Snakemake"))


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

    selection = [k for k in STEPS_KEYS if k in etapes_demandees(args.steps)]

    job_id = args.job_id or uuid.uuid4().hex[:12]
    output_dir = os.path.abspath(args.output_dir)
    proteins_dir = output_dir + "/proteins"
    results_dir = output_dir + "/results"

    print("===== PIPELINE START (job %s) =====" % job_id)
    print("Steps selected: %s" % ", ".join(selection))

    genomes = lister_genomes(os.path.abspath(args.genomes_dir), selection)   # toujours necessaire, rapide
    genomes_dir, _ = preparer_dossiers(genomes, output_dir, results_dir, proteins_dir)
    ecrire_versions(results_dir, selection, args)

    prevues = etapes_prevues(selection, list(genomes), proteins_dir)
    ajoutees = [k for k in prevues if k not in selection]
    if ajoutees:
        print("    (also running %s, needed by the selected steps)" % ", ".join(ajoutees))
    print("STEPS %s" % ",".join(prevues))
    sys.stdout.flush()
    # la configuration du Snakefile, gardee avec les resultats : snakemake peut etre relance a la main avec elle
    configfile = results_dir + "/snakemake_config.json"
    with open(configfile, "w") as fh:
        json.dump({"genomes": genomes_dir, "results": results_dir, "proteins": proteins_dir, "steps": prevues,
                   "job_id": job_id, "outgroup": args.outgroup}, fh, indent=2)
    lancer_snakemake(commande_snakemake(prevues, selection, output_dir, configfile, args.threads), prevues)

    ecrire_synthese(results_dir)
    print("===== PIPELINE DONE =====")


if __name__ == "__main__":
    main()
