# pipeline.py - Pipeline complet avec QC enrichi (CheckM2)
import subprocess
import glob
import os
import shutil

BASE     = "/home/oussama/stage"
GENOMES  = BASE + "/genomes"
PROTEINS = BASE + "/proteins"
RESULTS  = BASE + "/results"
SCRIPTS  = BASE + "/scripts"
EASYCG   = BASE + "/EasyCGTree4/EasyCGTree.v4.2-Linux"


def lire_checkm2():
    # Lance CheckM2 (dans son environnement) et renvoie un dictionnaire
    # {nom_genome: (completude, contamination)}
    print("    (CheckM2 en cours, cela prend quelques minutes...)")
    sortie = RESULTS + "/checkm2_out"
    subprocess.run(["conda", "run", "-n", "checkm2env",
                    "checkm2", "predict",
                    "--input", GENOMES,
                    "--output-directory", sortie,
                    "-x", "fas", "--threads", "4", "--force"])
    # lire le rapport produit par CheckM2
    qualite = {}
    rapport = sortie + "/quality_report.tsv"
    with open(rapport) as f:
        next(f)  # sauter la ligne d'en-tete
        for ligne in f:
            cols = ligne.strip().split("\t")
            nom = cols[0]
            completude = cols[1]
            contamination = cols[2]
            qualite[nom] = (completude, contamination)
    return qualite


def etape_qc():
    print(">>> Etape 1 : controle qualite (+ CheckM2)")
    # 1) recuperer completude/contamination via CheckM2
    qualite = lire_checkm2()
    # 2) calculer taille/GC/contigs et ecrire le tableau
    fichiers = sorted(glob.glob(GENOMES + "/*.fas"))
    rapport = open(RESULTS + "/qc.txt", "w")
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
        nom = os.path.basename(chemin).replace(".fas", "")
        comp, cont = qualite.get(nom, ("NA", "NA"))
        rapport.write("%s\t%d\t%.1f\t%d\t%s\t%s\n"
                      % (nom, total, 100 * gc / total, len(seqs), comp, cont))
    rapport.close()
    print("    -> qc.txt cree avec Completeness et Contamination")


def etape_proteines():
    print(">>> Etape 2 : prediction des proteines (Prodigal)")
    for chemin in glob.glob(GENOMES + "/*.fas"):
        nom = os.path.basename(chemin).replace(".fas", "")
        subprocess.run(["prodigal", "-i", chemin, "-a", PROTEINS + "/" + nom + ".faa", "-q"])
    print("    -> proteines ecrites dans", PROTEINS)


def etape_ani():
    print(">>> Etape 3 : calcul de l'ANI (FastANI)")
    liste = RESULTS + "/liste.txt"
    with open(liste, "w") as f:
        for g in glob.glob(GENOMES + "/*.fas"):
            f.write(g + "\n")
    subprocess.run(["fastANI", "--ql", liste, "--rl", liste,
                    "-o", RESULTS + "/ani_resultats.txt", "-t", "4"])
    print("    -> ani_resultats.txt cree")


def etape_aai():
    print(">>> Etape 4 : calcul de l'AAI (EzAAI)")
    aai_db = RESULTS + "/aai_db"
    os.makedirs(aai_db, exist_ok=True)
    for chemin in glob.glob(GENOMES + "/*.fas"):
        nom = os.path.basename(chemin).replace(".fas", "")
        subprocess.run(["EzAAI", "extract", "-i", chemin, "-o", aai_db + "/" + nom + ".db", "-l", nom])
    subprocess.run(["EzAAI", "calculate", "-i", aai_db, "-j", aai_db,
                    "-o", RESULTS + "/aai_resultats.tsv"])
    print("    -> aai_resultats.tsv cree")


def etape_pocp():
    print(">>> Etape 5 : calcul du POCP (script pocp.py + DIAMOND)")
    subprocess.run(["python", SCRIPTS + "/pocp.py", PROTEINS, RESULTS + "/pocp_out"])
    print("    -> pocp_out/pocp_matrice.tsv cree")


def etape_arbre():
    print(">>> Etape 6 : arbre phylogenomique (EasyCGTree, Perl)")
    entree = EASYCG + "/input_genomes"
    os.makedirs(entree, exist_ok=True)
    for g in glob.glob(GENOMES + "/*.fas"):
        shutil.copy(g, entree)
    subprocess.run(["perl", "EasyCGTree.pl", "-input", "input_genomes",
                    "-hmm", "bac120", "-tree", "sm", "-tree_app", "fasttree",
                    "-thread", "4"], cwd=EASYCG)
    arbre = EASYCG + "/input_genomes.bac120.120.supermatrix.fasttree.tree"
    if os.path.exists(arbre):
        shutil.copy(arbre, RESULTS + "/arbre.tree")
        print("    -> arbre.tree copie dans results")
    else:
        print("    !! arbre non trouve")


def etape_figures():
    print(">>> Etape 7 : figures (R : ggplot2 + ggtree)")
    subprocess.run(["Rscript", SCRIPTS + "/plot_all.R"])
    print("    -> figures PNG creees dans results")


print("===== DEBUT DU PIPELINE =====")
etape_qc()
etape_proteines()
etape_ani()
etape_aai()
etape_pocp()
etape_arbre()
etape_figures()
print("===== PIPELINE TERMINE =====")
