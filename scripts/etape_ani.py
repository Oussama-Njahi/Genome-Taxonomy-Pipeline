# etape_ani.py - Mon premier script d'orchestration (Seance 2)
import subprocess   # module qui permet de lancer d'autres programmes
import glob         # module pour lister les fichiers d'un dossier

# 1) Trouver tous les genomes .fas et ecrire leur liste dans un fichier
genomes = glob.glob("/home/oussama/stage/genomes/*.fas")
with open("/home/oussama/stage/results/liste.txt", "w") as f:
    for g in genomes:
        f.write(g + "\n")
print("Liste creee :", len(genomes), "genomes")

# 2) Lancer FastANI via subprocess (Python "tape" la commande a ta place)
print("Lancement de FastANI...")
subprocess.run([
    "fastANI",
    "--ql", "/home/oussama/stage/results/liste.txt",
    "--rl", "/home/oussama/stage/results/liste.txt",
    "-o",  "/home/oussama/stage/results/ani_resultats.txt",
    "-t",  "4"
])
print("Termine ! Resultat dans results/ani_resultats.txt")

