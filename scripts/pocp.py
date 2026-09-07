#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pocp.py  --  Calcul du POCP (Percentage of Conserved Proteins)
============================================================================
Le POCP sert a delimiter les GENRES bacteriens (Qin et al., 2014).
Regle : deux genomes avec POCP > 50 % appartiennent au meme genre.

Principe (pour chaque paire de genomes A et B) :
  POCP = (proteines_conservees_A->B + proteines_conservees_B->A)
         / (proteines_totales_A + proteines_totales_B)  x 100

Une proteine est "conservee" si elle a, dans l'autre genome, un equivalent
trouve par DIAMOND avec : identite >= 40 %, couverture >= 50 %, e-value <= 1e-5.

USAGE :
    python pocp.py <dossier_proteines> <dossier_sortie>
        <dossier_proteines> : contient un fichier .faa par genome (proteines).
        <dossier_sortie>    : sera cree ; contiendra la matrice POCP.
============================================================================
"""

import os, sys, glob, subprocess, itertools, csv

# ----- 0) Lire les arguments de la ligne de commande -----
if len(sys.argv) != 3:
    print("Usage : python pocp.py <dossier_proteines> <dossier_sortie>")
    sys.exit(1)

prot_dir = sys.argv[1]                 # ou sont les .faa
out_dir  = sys.argv[2]                 # ou ecrire les resultats
db_dir   = os.path.join(out_dir, "dbs")         # bases DIAMOND
aln_dir  = os.path.join(out_dir, "alignements")  # resultats bruts diamond
for d in (out_dir, db_dir, aln_dir):
    os.makedirs(d, exist_ok=True)

# ----- 1) Lister les proteomes (.faa) -----
faa_files = sorted(glob.glob(os.path.join(prot_dir, "*.faa")))
if not faa_files:
    print("ERREUR : aucun fichier .faa trouve dans", prot_dir)
    sys.exit(1)
names      = [os.path.basename(f)[:-4] for f in faa_files]  # nom sans ".faa"
faa_byname = dict(zip(names, faa_files))
print(len(faa_files), "proteomes detectes.")

# ----- 2) Compter les proteines + creer une base DIAMOND par genome -----
total_prot = {}
for name, f in faa_byname.items():
    total_prot[name] = sum(1 for line in open(f) if line.startswith(">"))
    subprocess.run(["diamond", "makedb", "--in", f,
                    "-d", os.path.join(db_dir, name), "--quiet"], check=True)
print("Bases DIAMOND creees.")

# ----- 3) Fonction : nombre de proteines de A conservees dans B -----
def conserved(query_name, subject_name):
    """Aligne les proteines de A (query) contre B (subject) et compte
    celles qui passent les seuils (>=40% identite, >=50% couverture)."""
    out_tsv = os.path.join(aln_dir, query_name + "__vs__" + subject_name + ".tsv")
    subprocess.run([
        "diamond", "blastp",
        "-q", faa_byname[query_name],
        "-d", os.path.join(db_dir, subject_name),
        "-o", out_tsv,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen",
        "--evalue", "1e-5", "--max-target-seqs", "1", "--quiet"
    ], check=True)
    conserved_ids = set()
    with open(out_tsv) as fh:
        for line in fh:
            qid, sid, pident, length, qlen = line.rstrip("\n").split("\t")
            identite   = float(pident)
            couverture = float(length) / float(qlen)     # part de la proteine alignee
            if identite >= 40.0 and couverture >= 0.5:
                conserved_ids.add(qid)                   # proteine conservee
    return len(conserved_ids)

# ----- 4) Calculer le POCP pour toutes les paires -----
matrix = {a: {b: 100.0 for b in names} for a in names}   # diagonale = 100 %
pairs  = list(itertools.combinations(names, 2))
for i, (a, b) in enumerate(pairs, 1):
    print("  paire %d/%d : %s  vs  %s" % (i, len(pairs), a, b))
    c_ab = conserved(a, b)                # proteines de A conservees dans B
    c_ba = conserved(b, a)                # proteines de B conservees dans A
    pocp = (c_ab + c_ba) / (total_prot[a] + total_prot[b]) * 100.0
    matrix[a][b] = matrix[b][a] = round(pocp, 2)

# ----- 5) Ecrire la matrice POCP -----
out_matrix = os.path.join(out_dir, "pocp_matrice.tsv")
with open(out_matrix, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([""] + names)
    for a in names:
        w.writerow([a] + [matrix[a][b] for b in names])

print("\nTERMINE. Matrice POCP ecrite dans :", out_matrix)
