import subprocess, glob, os, shutil
BASE   = "/home/oussama/stage"
GENOMES= BASE + "/genomes"
RESULTS= BASE + "/results"
EASYCG = BASE + "/EasyCGTree4/EasyCGTree.v4.2-Linux"

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

etape_arbre()
