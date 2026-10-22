"""pipeline.py - Point d'entree en ligne de commande du pipeline.

Lance le meme pipeline que l'interface web (pipeline_web.py, qui execute les
etapes avec Snakemake : workflow/Snakefile) sur les dossiers
habituels du projet : genomes/ en entree (.fas, .fasta ou .fna, au moins 5),
proteins/ et results/ en sortie. Les 9 etapes n'existent donc qu'en un seul
exemplaire, avec les memes verifications d'erreurs et les memes figures en
ligne de commande et sur le site.

Usage :
  python scripts/pipeline.py                      # les 9 etapes
  python scripts/pipeline.py --steps qc,ani,dddh  # seulement certaines
  python scripts/pipeline.py --threads 8          # threads pour FastANI/EasyCGTree
  python scripts/pipeline.py --outgroup Staphylococcus-aureus-NCTC8325
                                                  # racine de l'arbre (defaut : point median)

Etapes (--steps, separees par des virgules) :
  qc, checkm2, proteines, ani, dddh, aai, pocp, arbre, figures
  (pocp a besoin des proteines : l'etape proteines tourne aussi si elles manquent ;
   qc.txt est refait apres checkm2 pour prendre ses valeurs)
"""
import os
import subprocess
import sys

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPTS)

if __name__ == "__main__":
    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print(__doc__)
        sys.exit(0)
    # Les options (--steps, --threads, --outgroup) sont transmises telles quelles.
    cmd = [sys.executable, os.path.join(SCRIPTS, "pipeline_web.py"),
           "--genomes-dir", os.path.join(BASE, "genomes"),
           "--output-dir", BASE] + sys.argv[1:]
    sys.exit(subprocess.call(cmd))
