#!/usr/bin/env python3
# qc.py - Rapport qualite des genomes (results/qc.txt), etape 'qc' du pipeline (voir workflow/Snakefile).
#
# Taille, %GC, nombre de contigs et N50 viennent de lire_fasta (pipeline_web.py). Le %GC est calcule sur les seules
# bases A/C/G/T, pour que les N d'un assemblage ne le fassent pas baisser. Completeness/Contamination viennent du
# rapport de CheckM2 (checkm2_out_fna/quality_report.tsv) quand il existe, et valent NA sinon : le Snakefile
# relance cette etape apres CheckM2, si bien que CheckM2 n'a plus a reecrire qc.txt lui-meme.
#
# Usage : python qc.py <dossier_genomes> <dossier_resultats>
import glob
import os
import sys

import pipeline_web as pw


def lire_checkm2(rapport):
    """{genome: (completude, contamination)} lus dans le quality_report.tsv de CheckM2, vide s'il n'existe pas.
    Les colonnes sont reperees par leur nom, pas par leur position, au cas ou une version de CheckM2 les deplace."""
    qualite = {}
    if not os.path.exists(rapport):
        return qualite
    with open(rapport) as fh:
        entete = fh.readline().rstrip("\n").split("\t")
        i_nom, i_comp, i_cont = (entete.index(c) for c in ("Name", "Completeness", "Contamination"))
        for ligne in fh:
            cols = ligne.rstrip("\n").split("\t")
            qualite[cols[i_nom]] = (cols[i_comp], cols[i_cont])
    return qualite


def ecrire_qc(genomes_dir, results_dir):
    qualite = lire_checkm2(results_dir + "/checkm2_out_fna/quality_report.tsv")
    with open(results_dir + "/qc.txt", "w") as rapport:
        rapport.write("Genome\tTaille(pb)\tGC(%)\tContigs\tN50\tCompleteness(%)\tContamination(%)\n")
        for chemin in sorted(glob.glob(genomes_dir + "/*.fna")):
            nom = os.path.splitext(os.path.basename(chemin))[0]
            total, gc, acgt, sequences, _, n50 = pw.lire_fasta(chemin)
            taux = "%.1f" % (100.0 * gc / acgt) if acgt else "NA"
            completude, contamination = qualite.get(nom, ("NA", "NA"))
            rapport.write("%s\t%d\t%s\t%d\t%d\t%s\t%s\n" % (nom, total, taux, sequences, n50, completude, contamination))
    if qualite:
        print("    -> qc.txt created, with CheckM2's completeness and contamination")
    else:
        print("    -> qc.txt created (Completeness/Contamination = NA; also run the 'checkm2' step to fill them in)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python qc.py <genomes_dir> <results_dir>")
    ecrire_qc(sys.argv[1], sys.argv[2])
