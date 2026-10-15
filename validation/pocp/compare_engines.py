#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_engines.py  --  POCP : DIAMOND (rapide, ultra-sensible) contre BLASTP
============================================================================
Qin et al. (2014) ont defini le POCP avec BLASTP. Le pipeline utilise DIAMOND,
beaucoup plus rapide. Ce script calcule le POCP d'une paire de proteomes avec
trois moteurs, EXACTEMENT les memes criteres que scripts/pocp.py (meilleur
hit, e-value <= 1e-5, identite >= 40 %, alignement >= 50 % de la requete),
pour mesurer l'ecart du a l'outil d'alignement seul.

USAGE (dans l'environnement conda du pipeline) :
    python validation/pocp/compare_engines.py <A.faa> <B.faa> [threads]
============================================================================
"""
import os
import subprocess
import sys
import tempfile
import time

COLONNES = ["qseqid", "sseqid", "pident", "length", "qlen"]


def nombre(faa):
    return sum(1 for ligne in open(faa) if ligne.startswith(">"))


def conservees(tsv):
    """Proteines requete dont le meilleur hit passe les seuils (une seule ligne lue par requete)."""
    vues, ok = set(), set()
    with open(tsv) as fh:
        for ligne in fh:
            qid, _, pident, longueur, qlen = ligne.rstrip("\n").split("\t")
            if qid in vues:
                continue
            vues.add(qid)
            if float(pident) >= 40.0 and float(longueur) / float(qlen) >= 0.5:
                ok.add(qid)
    return len(ok)


def diamond(requete, sujet, sortie, threads, mode, tmp):
    base = os.path.join(tmp, "dmnd_" + os.path.basename(sujet))
    subprocess.run(["diamond", "makedb", "--in", sujet, "-d", base, "--quiet"], check=True)
    cmd = ["diamond", "blastp", "-q", requete, "-d", base, "-o", sortie, "--outfmt", "6"] + COLONNES + \
          ["--evalue", "1e-5", "--max-target-seqs", "1", "--threads", str(threads), "--quiet"]
    subprocess.run(cmd + ([mode] if mode else []), check=True)


def blastp(requete, sujet, sortie, threads, mode, tmp):
    base = os.path.join(tmp, "blast_" + os.path.basename(sujet))
    subprocess.run(["makeblastdb", "-in", sujet, "-dbtype", "prot", "-out", base],
                   check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["blastp", "-query", requete, "-db", base, "-out", sortie,
                    "-outfmt", "6 " + " ".join(COLONNES), "-evalue", "1e-5",
                    "-max_target_seqs", "1", "-max_hsps", "1", "-num_threads", str(threads)],
                   check=True, stderr=subprocess.DEVNULL)


MOTEURS = [("DIAMOND default (fast)", diamond, None),
           ("DIAMOND --ultra-sensitive", diamond, "--ultra-sensitive"),
           ("BLASTP (Qin et al. 2014)", blastp, None)]

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage : python compare_engines.py <A.faa> <B.faa> [threads]")
        sys.exit(1)
    a, b = sys.argv[1], sys.argv[2]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    total = nombre(a) + nombre(b)
    print("%s (%d proteins)  vs  %s (%d proteins)" % (os.path.basename(a), nombre(a), os.path.basename(b), nombre(b)))
    with tempfile.TemporaryDirectory() as tmp:
        for nom, moteur, mode in MOTEURS:
            debut = time.time()
            moteur(a, b, os.path.join(tmp, "ab.tsv"), threads, mode, tmp)
            moteur(b, a, os.path.join(tmp, "ba.tsv"), threads, mode, tmp)
            c = conservees(os.path.join(tmp, "ab.tsv")) + conservees(os.path.join(tmp, "ba.tsv"))
            print("  %-28s POCP = %6.2f %%   (%d conserved, %.0f s)" % (nom, 100.0 * c / total, c, time.time() - debut))
