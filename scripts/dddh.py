#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dddh.py  --  Calcul du dDDH (digital DNA-DNA Hybridization), formule 2 (GBDP)
============================================================================
Le dDDH sert a delimiter les ESPECES bacteriennes (seuil ~70 %), en
complement de l'ANI -- calcule ici de facon INDEPENDANTE de l'ANI (pas une
conversion depuis le chiffre d'ANI), via alignement BLASTN direct entre
chaque paire de genomes.

Methode et coefficients publies dans :
  Auch AF, Klenk HP, Goker M (2010). "Standard operating procedure for
  calculating genome-to-genome distances based on high-scoring segment
  pairs." Stand Genomic Sci 2(1):142-148.  -- formule (2), recommandee
  (la plus robuste pour des genomes incomplets / MAG fragmentes).

Principe (pour chaque paire de genomes A et B) :
  1. BLASTN A -> B, et B -> A (deux sens, la recherche est asymetrique).
  2. Ne garder que les alignements (HSP) avec e-value <= 1e-2.
  3. Corriger les chevauchements entre HSP sur une meme sequence
     ("greedy-with-trimming") : une region dupliquee ne doit etre
     comptee qu'une fois -- on garde en priorite les HSP les plus
     longues, et on ne compte des HSP suivantes que la portion pas
     deja couverte.
  4. d = 1 - (identites dans les HSP, 2 sens) / (longueur totale des
           HSP apres trimming, 2 sens)
  5. similarite estimee (~%dDDH) = 90.3998 - 438.3134 * d
     (coefficients publies, NCBI BLAST, formule 2)
  6. Seuil espece : d > 0.0412  =>  especes differentes (analogue a
     <70 % DDH). d <= 0.0412 => meme espece.

Limite connue : le "greedy-with-trimming" implemente ici suppose une
identite uniforme le long de chaque HSP tronquee (on n'a pas acces au
detail base-par-base de l'alignement, seulement au %identite global de
chaque HSP). C'est une approximation raisonnable et couramment admise,
mais ce n'est pas un reimplementation bit-a-bit du code interne de GGDC.

USAGE :
    python dddh.py <dossier_genomes> <dossier_sortie>
        <dossier_genomes> : contient un fichier .fas/.fasta/.fna par genome.
        <dossier_sortie>  : sera cree ; contiendra les matrices dDDH.
============================================================================
"""

import os, sys, glob, subprocess, itertools, csv

# ----- 0) Lire les arguments de la ligne de commande -----
if len(sys.argv) != 3:
    print("Usage : python dddh.py <dossier_genomes> <dossier_sortie>")
    sys.exit(1)

genomes_dir = sys.argv[1]
out_dir = sys.argv[2]
db_dir = os.path.join(out_dir, "dbs")
aln_dir = os.path.join(out_dir, "alignements")
for d in (out_dir, db_dir, aln_dir):
    os.makedirs(d, exist_ok=True)

EVALUE_MAX = 1e-2
# Coefficients publies (Auch, Klenk & Goker 2010, Table 4 -- NCBI BLAST, formule 2)
INTERCEPT = 90.3998
SLOPE = -438.3134
SEUIL_ESPECE = 0.0412  # d au-dessus de ce seuil => especes differentes (~<70% DDH)

# ----- 1) Lister les genomes -----
genome_files = sorted(glob.glob(os.path.join(genomes_dir, "*.fas")) +
                       glob.glob(os.path.join(genomes_dir, "*.fasta")) +
                       glob.glob(os.path.join(genomes_dir, "*.fna")))
if not genome_files:
    print("ERREUR : aucun fichier de genome trouve dans", genomes_dir)
    sys.exit(1)
names = [os.path.splitext(os.path.basename(f))[0] for f in genome_files]
fasta_byname = dict(zip(names, genome_files))
print(len(genome_files), "genomes detectes.")

# ----- 2) Construire une base BLAST par genome -----
for name, fasta in fasta_byname.items():
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl",
                     "-out", os.path.join(db_dir, name)],
                    check=True, stdout=subprocess.DEVNULL)
print("Bases BLAST creees.")


def hsps_diriges(query_name, subject_name):
    """BLASTN query->subject, renvoie la liste des HSP filtrees (e-value
    <= 1e-2) : (sequence_query, debut, fin, longueur, %identite)."""
    out_tsv = os.path.join(aln_dir, query_name + "__vs__" + subject_name + ".tsv")
    subprocess.run([
        "blastn",
        "-query", fasta_byname[query_name],
        "-db", os.path.join(db_dir, subject_name),
        "-out", out_tsv,
        "-outfmt", "6 qseqid qstart qend length pident",
        "-evalue", str(EVALUE_MAX),
        "-max_target_seqs", "100000",
        "-num_threads", "4",
    ], check=True)
    hsps = []
    with open(out_tsv) as fh:
        for ligne in fh:
            qseqid, qstart, qend, longueur, pident = ligne.rstrip("\n").split("\t")
            qstart, qend = int(qstart), int(qend)
            if qstart > qend:                      # brin -, coordonnees inversees
                qstart, qend = qend, qstart
            hsps.append((qseqid, qstart, qend, int(longueur), float(pident)))
    return hsps


def portions_non_couvertes(debut, fin, intervalles_couverts):
    """Renvoie la liste des sous-intervalles de [debut,fin] qui ne
    chevauchent aucun intervalle deja couvert."""
    libres = [(debut, fin)]
    for c_debut, c_fin in intervalles_couverts:
        nouveaux = []
        for l_debut, l_fin in libres:
            if c_fin < l_debut or c_debut > l_fin:
                nouveaux.append((l_debut, l_fin))          # pas de chevauchement
                continue
            if c_debut > l_debut:
                nouveaux.append((l_debut, c_debut - 1))    # portion avant le chevauchement
            if c_fin < l_fin:
                nouveaux.append((c_fin + 1, l_fin))        # portion apres le chevauchement
        libres = [seg for seg in nouveaux if seg[0] <= seg[1]]
    return libres


def longueur_et_identites_apres_trimming(hsps):
    """Corrige les chevauchements entre HSP sur une meme sequence query
    (greedy-with-trimming, Auch et al. 2010) : on traite d'abord les HSP
    les plus longues, et pour chaque HSP suivante on ne compte que la
    portion de son intervalle pas deja couverte par une HSP precedente
    (donc plus longue ou egale). Renvoie (longueur_totale, identites_totales)."""
    hsps_triees = sorted(hsps, key=lambda h: h[3], reverse=True)
    couverture = {}   # qseqid -> liste d'intervalles (debut, fin) deja couverts
    longueur_totale = 0
    identites_totales = 0.0
    for qseqid, debut, fin, longueur, pident in hsps_triees:
        intervalles = couverture.setdefault(qseqid, [])
        for p_debut, p_fin in portions_non_couvertes(debut, fin, intervalles):
            p_longueur = p_fin - p_debut + 1
            longueur_totale += p_longueur
            # identite supposee uniforme le long de la HSP (voir limite en en-tete)
            identites_totales += pident / 100.0 * p_longueur
            intervalles.append((p_debut, p_fin))
    return longueur_totale, identites_totales


def dddh_paire(a, b):
    """Calcule (d, similarite) pour la paire de genomes (a, b), ou None
    si aucun alignement significatif n'a ete trouve entre les deux."""
    long_ab, ident_ab = longueur_et_identites_apres_trimming(hsps_diriges(a, b))
    long_ba, ident_ba = longueur_et_identites_apres_trimming(hsps_diriges(b, a))
    longueur_totale = long_ab + long_ba
    identites_totales = ident_ab + ident_ba
    if longueur_totale == 0:
        return None
    d = 1 - identites_totales / longueur_totale
    similarite = INTERCEPT + SLOPE * d
    return d, similarite


# ----- 3) Calculer le dDDH pour toutes les paires -----
matrix_d = {a: {b: 0.0 for b in names} for a in names}          # diagonale = 0 (identique)
matrix_s = {a: {b: 100.0 for b in names} for a in names}        # diagonale = 100 %
pairs = list(itertools.combinations(names, 2))
for i, (a, b) in enumerate(pairs, 1):
    print("  paire %d/%d : %s  vs  %s" % (i, len(pairs), a, b))
    resultat = dddh_paire(a, b)
    if resultat is None:
        print("    !! aucun alignement significatif trouve, paire ignoree")
        continue
    d, s = resultat
    matrix_d[a][b] = matrix_d[b][a] = round(d, 4)
    matrix_s[a][b] = matrix_s[b][a] = round(s, 2)

# ----- 4) Ecrire les matrices -----
out_matrix_d = os.path.join(out_dir, "dddh_distance.tsv")
out_matrix_s = os.path.join(out_dir, "dddh_similarite.tsv")
for out_matrix, matrix in ((out_matrix_d, matrix_d), (out_matrix_s, matrix_s)):
    with open(out_matrix, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([""] + names)
        for a in names:
            w.writerow([a] + [matrix[a][b] for b in names])

print("\nTERMINE.")
print("  -> distance GBDP (formule 2)     :", out_matrix_d)
print("  -> similarite estimee (~%%dDDH)   :", out_matrix_s)
print("  Seuil espece : d > %.4f (similarite < ~70%%) => especes differentes" % SEUIL_ESPECE)
