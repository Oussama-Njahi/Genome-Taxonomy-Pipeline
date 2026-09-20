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
#
# Progress is reported on stdout as lines "PROGRESS <n>/9 <label>" so a caller
# (e.g. the FastAPI backend) can parse them to drive a progress bar.

import argparse
import glob
import os
import shutil
import subprocess
import sys
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

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPTS)
EASYCG = BASE + "/EasyCGTree4/EasyCGTree.v4.2-Linux"


def announce(step_index, extra=""):
    label = STEPS[step_index - 1]
    print("PROGRESS %d/%d %s%s" % (step_index, len(STEPS), label, (" - " + extra) if extra else ""))
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


def lister_genomes(genomes_dir):
    """Toujours necessaire : trouve les genomes et verifie qu'il y en a assez.
    Rapide (juste un glob + un comptage) -> peut rester obligatoire sans cout."""
    fichiers = sorted(glob.glob(genomes_dir + "/*.fas") +
                       glob.glob(genomes_dir + "/*.fasta") +
                       glob.glob(genomes_dir + "/*.fna"))
    if not fichiers:
        fail("no .fas/.fasta/.fna genome files found in %s" % genomes_dir)
    if len(fichiers) < 5:
        fail("at least 5 genomes are required (EasyCGTree cannot build a tree with fewer); got %d" % len(fichiers))
    return fichiers


def lire_checkm2(genomes_dir, results_dir, fichiers):
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
             "-x", ext, "--threads", "4", "--force"])

        rapport = sortie + "/quality_report.tsv"
        with open(rapport) as f:
            next(f)
            for ligne in f:
                cols = ligne.strip().split("\t")
                nom, completude, contamination = cols[0], cols[1], cols[2]
                qualite[nom] = (completude, contamination)

    return qualite


def etape_qc(fichiers, results_dir):
    """Stats de base uniquement (taille, GC%, contigs) - rapide, pas de CheckM2.
    Completeness/Contamination sont ecrits a "NA" ; lancer aussi l'etape
    'checkm2' pour les remplir (voir etape_checkm2 ci-dessous)."""
    announce(1)
    with open(results_dir + "/qc.txt", "w") as rapport:
        rapport.write("Genome\tTaille(pb)\tGC(%)\tContigs\tCompleteness(%)\tContamination(%)\n")
        for chemin in fichiers:
            seqs = []
            seq = ""
            for ligne in open(chemin):
                if ligne.startswith(">"):
                    if seq:
                        seqs.append(seq)
                    seq = ""
                else:
                    seq += ligne.strip()
            if seq:
                seqs.append(seq)
            total = sum(len(s) for s in seqs)
            gc = sum(s.upper().count("G") + s.upper().count("C") for s in seqs)
            nom = os.path.splitext(os.path.basename(chemin))[0]
            rapport.write("%s\t%d\t%.1f\t%d\tNA\tNA\n" % (nom, total, 100 * gc / total, len(seqs)))
    print("    -> qc.txt created (Completeness/Contamination = NA; also run 'checkm2' step to fill them in)")


def etape_checkm2(genomes_dir, results_dir, fichiers):
    """Lance CheckM2 et met a jour les colonnes Completeness/Contamination
    de qc.txt (qui doit deja exister, ecrit par l'etape 'qc' - verifie par
    l'appelant dans main(), comme pour la dependance proteines -> pocp)."""
    announce(2)
    qc_path = results_dir + "/qc.txt"
    qualite = lire_checkm2(genomes_dir, results_dir, fichiers)
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


def etape_dddh(genomes_dir, results_dir):
    announce(5)
    run(["python", SCRIPTS + "/dddh.py", genomes_dir, results_dir + "/dddh_out"])
    print("    -> dddh_out/dddh_distance.tsv and dddh_similarite.tsv created")


def etape_aai(fichiers, results_dir):
    announce(6)
    aai_db = results_dir + "/aai_db"
    os.makedirs(aai_db, exist_ok=True)
    for chemin in fichiers:
        nom = os.path.splitext(os.path.basename(chemin))[0]
        run(["EzAAI", "extract", "-i", chemin, "-o", aai_db + "/" + nom + ".db", "-l", nom])
    run(["EzAAI", "calculate", "-i", aai_db, "-j", aai_db,
         "-o", results_dir + "/aai_resultats.tsv"])
    print("    -> aai_resultats.tsv created")


def etape_pocp(proteins_dir, results_dir):
    announce(7)
    run(["python", SCRIPTS + "/pocp.py", proteins_dir, results_dir + "/pocp_out"])
    print("    -> pocp_out/pocp_matrice.tsv created")


def etape_arbre(fichiers, results_dir, job_id, threads):
    announce(8)
    input_name = "job_" + job_id
    entree = EASYCG + "/" + input_name
    os.makedirs(entree, exist_ok=True)
    for g in fichiers:
        shutil.copy(g, entree)
    run(["perl", "EasyCGTree.pl", "-input", input_name,
         "-hmm", "bac120", "-tree", "sm", "-tree_app", "fasttree",
         "-thread", str(threads)], cwd=EASYCG)
    arbre = EASYCG + "/" + input_name + ".bac120.120.supermatrix.fasttree.tree"
    if os.path.exists(arbre):
        shutil.copy(arbre, results_dir + "/arbre.tree")
        print("    -> arbre.tree copied to results")
    else:
        fail("tree file not found (%s) - EasyCGTree may have failed" % arbre)
    cleanup_easycgtree(input_name)


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
    parser.add_argument("--threads", type=int, default=4, help="Threads for fastANI/EasyCGTree")
    parser.add_argument("--steps", type=str, default="all",
                         help="Comma-separated steps to run (" + ", ".join(STEPS_KEYS) + "). Default: all")
    args = parser.parse_args()

    selection = etapes_demandees(args.steps)

    job_id = args.job_id or uuid.uuid4().hex[:12]
    genomes_dir = os.path.abspath(args.genomes_dir)
    proteins_dir = os.path.abspath(args.output_dir) + "/proteins"
    results_dir = os.path.abspath(args.output_dir) + "/results"
    os.makedirs(proteins_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print("===== PIPELINE START (job %s) =====" % job_id)
    print("Steps selected: %s" % ", ".join(selection))

    fichiers = lister_genomes(genomes_dir)   # toujours necessaire, rapide

    if "qc" in selection:
        etape_qc(fichiers, results_dir)
    if "checkm2" in selection:
        if not os.path.exists(results_dir + "/qc.txt"):
            print("    !! WARNING: qc.txt not found in %s" % results_dir)
            print("    !! The 'qc' step must be run at least once before 'checkm2'.")
            print("    !! Re-run with --steps qc,checkm2")
        else:
            etape_checkm2(genomes_dir, results_dir, fichiers)
    if "proteines" in selection:
        etape_proteines(fichiers, proteins_dir)
    if "ani" in selection:
        etape_ani(fichiers, results_dir, args.threads)
    if "dddh" in selection:
        etape_dddh(genomes_dir, results_dir)
    if "aai" in selection:
        etape_aai(fichiers, results_dir)
    if "pocp" in selection:
        if not os.path.isdir(proteins_dir) or not glob.glob(proteins_dir + "/*.faa"):
            print("    !! WARNING: no protein files found in %s" % proteins_dir)
            print("    !! The 'proteines' step must be run at least once before POCP.")
            print("    !! Re-run with --steps proteines,pocp")
        else:
            etape_pocp(proteins_dir, results_dir)
    if "arbre" in selection:
        etape_arbre(fichiers, results_dir, job_id, args.threads)
    if "figures" in selection:
        etape_figures(results_dir)

    print("===== PIPELINE DONE =====")


if __name__ == "__main__":
    main()
