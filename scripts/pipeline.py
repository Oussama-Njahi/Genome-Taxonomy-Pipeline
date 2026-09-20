# pipeline.py - Pipeline avec choix des etapes a executer
import subprocess
import glob
import os
import shutil
import argparse
import sys

BASE     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMES  = BASE + "/genomes"
PROTEINS = BASE + "/proteins"
RESULTS  = BASE + "/results"
SCRIPTS  = BASE + "/scripts"
EASYCG   = BASE + "/EasyCGTree4/EasyCGTree.v4.2-Linux"

# La liste officielle des noms d'etapes valides, dans l'ordre logique
ETAPES_VALIDES = ["qc", "checkm2", "proteines", "ani", "dddh", "aai", "pocp", "arbre", "figures"]


def lire_checkm2():
    print("    (CheckM2 en cours, cela prend quelques minutes...)")
    sortie = RESULTS + "/checkm2_out"
    subprocess.run(["conda", "run", "-n", "checkm2env",
                    "checkm2", "predict",
                    "--input", GENOMES,
                    "--output-directory", sortie,
                    "-x", "fas", "--threads", "4", "--force"])
    qualite = {}
    rapport = sortie + "/quality_report.tsv"
    with open(rapport) as f:
        next(f)
        for ligne in f:
            cols = ligne.strip().split("\t")
            qualite[cols[0]] = (cols[1], cols[2])
    return qualite


def etape_qc():
    print(">>> Etape : controle qualite basique (taille, GC%, contigs)")
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
        rapport.write("%s\t%d\t%.1f\t%d\tNA\tNA\n"
                      % (nom, total, 100 * gc / total, len(seqs)))
    rapport.close()
    print("    -> qc.txt cree (Completeness/Contamination = NA ; lancer aussi 'checkm2' pour les avoir)")


def etape_checkm2():
    print(">>> Etape : CheckM2 (completude et contamination)")
    qc_path = RESULTS + "/qc.txt"
    if not os.path.exists(qc_path):
        print("    !! ATTENTION : qc.txt introuvable.")
        print("    !! L'etape 'qc' doit avoir ete executee au moins une fois avant 'checkm2'.")
        print("    !! Relance avec : python pipeline.py --steps qc,checkm2")
        return
    qualite = lire_checkm2()
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
    print("    -> qc.txt mis a jour avec Completeness et Contamination")


def etape_proteines():
    print(">>> Etape : prediction des proteines (Prodigal)")
    for chemin in glob.glob(GENOMES + "/*.fas"):
        nom = os.path.basename(chemin).replace(".fas", "")
        subprocess.run(["prodigal", "-i", chemin, "-a", PROTEINS + "/" + nom + ".faa", "-q"])
    print("    -> proteines ecrites dans", PROTEINS)


def etape_ani():
    print(">>> Etape : calcul de l'ANI (FastANI)")
    liste = RESULTS + "/liste.txt"
    with open(liste, "w") as f:
        for g in glob.glob(GENOMES + "/*.fas"):
            f.write(g + "\n")
    subprocess.run(["fastANI", "--ql", liste, "--rl", liste,
                    "-o", RESULTS + "/ani_resultats.txt", "-t", "4"])
    print("    -> ani_resultats.txt cree")


def etape_dddh():
    print(">>> Etape : calcul du dDDH (script dddh.py + BLASTN)")
    subprocess.run(["python", SCRIPTS + "/dddh.py", GENOMES, RESULTS + "/dddh_out"])
    print("    -> dddh_out/dddh_distance.tsv et dddh_similarite.tsv crees")


def etape_aai():
    print(">>> Etape : calcul de l'AAI (EzAAI)")
    aai_db = RESULTS + "/aai_db"
    os.makedirs(aai_db, exist_ok=True)
    for chemin in glob.glob(GENOMES + "/*.fas"):
        nom = os.path.basename(chemin).replace(".fas", "")
        subprocess.run(["EzAAI", "extract", "-i", chemin, "-o", aai_db + "/" + nom + ".db", "-l", nom])
    subprocess.run(["EzAAI", "calculate", "-i", aai_db, "-j", aai_db,
                    "-o", RESULTS + "/aai_resultats.tsv"])
    print("    -> aai_resultats.tsv cree")


def etape_pocp():
    print(">>> Etape : calcul du POCP (script pocp.py + DIAMOND)")
    if not os.path.isdir(PROTEINS) or not glob.glob(PROTEINS + "/*.faa"):
        print("    !! ATTENTION : aucun fichier de proteines trouve dans", PROTEINS)
        print("    !! L'etape 'proteines' doit avoir ete executee au moins une fois avant POCP.")
        print("    !! Relance avec : python pipeline.py --steps proteines,pocp")
        return
    subprocess.run(["python", SCRIPTS + "/pocp.py", PROTEINS, RESULTS + "/pocp_out"])
    print("    -> pocp_out/pocp_matrice.tsv cree")


def etape_arbre():
    print(">>> Etape : arbre phylogenomique (EasyCGTree, Perl)")
    if not os.path.exists(EASYCG + "/EasyCGTree.pl"):
        print("    !! ATTENTION : EasyCGTree4 introuvable dans", EASYCG)
        print("    !! Cet outil ne s'installe pas via conda, il faut l'installer manuellement")
        print("    !! (voir la section Installation du README).")
        return
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
    print(">>> Etape : figures (R : ggplot2 + ggtree)")
    subprocess.run(["Rscript", SCRIPTS + "/plot_all.R", RESULTS])
    print("    -> figures PNG creees dans results")


# ----- Fonctions liees aux arguments de ligne de commande -----

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pipeline de taxonomie genomique (QC, CheckM2, ANI, dDDH, AAI, POCP, arbre, figures)."
    )
    parser.add_argument(
        "--steps",
        type=str,
        default="all",
        help=("Etapes a executer, separees par des virgules. "
              "Choix possibles : " + ", ".join(ETAPES_VALIDES) + ". "
              "Utilise 'all' pour tout executer (comportement par defaut).")
    )
    return parser.parse_args()


def etapes_demandees(args):
    if args.steps == "all":
        return ETAPES_VALIDES
    choisies = [s.strip() for s in args.steps.split(",")]
    for s in choisies:
        if s not in ETAPES_VALIDES:
            print("Etape inconnue :", s)
            print("Etapes valides :", ", ".join(ETAPES_VALIDES))
            sys.exit(1)
    return choisies


FONCTIONS = {
    "qc": etape_qc,
    "checkm2": etape_checkm2,
    "proteines": etape_proteines,
    "ani": etape_ani,
    "dddh": etape_dddh,
    "aai": etape_aai,
    "pocp": etape_pocp,
    "arbre": etape_arbre,
    "figures": etape_figures,
}


if __name__ == "__main__":
    args = parse_args()
    selection = etapes_demandees(args)

    print("===== DEBUT DU PIPELINE =====")
    print("Etapes selectionnees :", ", ".join(selection))

    for nom_etape in ETAPES_VALIDES:       # on respecte toujours l'ordre logique
        if nom_etape in selection:
            FONCTIONS[nom_etape]()

    print("===== PIPELINE TERMINE =====")
