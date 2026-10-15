#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_dddh.py  --  Validation de scripts/dddh.py sur un jeu de reference
============================================================================
Le dDDH du pipeline est une implementation maison (Auch et al. 2010,
formule 2). Ce script mesure a quel point ses DECISIONS d'espece sont
fiables, sur des genomes de reference choisis pour encadrer le seuil
(memes especes, especes tres proches, especes eloignees) : benchmark.tsv.

Pour chaque groupe de benchmark.tsv :
  1. telecharge les genomes depuis NCBI (une seule fois) ;
  2. lance scripts/dddh.py tel quel (le script valide est celui du pipeline) ;
  3. lance FastANI sur les memes genomes.

Puis compare, paire par paire, la decision du dDDH (d <= 0.0412) :
  - a celle de l'ANI (>= 95 %), critere d'espece independant ;
  - a GGDC (dDDH formule 2 >= 70 %), SI ggdc_reference.tsv a ete rempli avec
    les valeurs obtenues sur https://ggdc.dsmz.de (voir README.md).

Sorties : pairs.tsv (toutes les paires), RESULTS.md (resume), et
ggdc_reference.tsv (modele a remplir) s'il n'existe pas encore.

USAGE (dans l'environnement conda du pipeline) :
    python validation/dddh/validate_dddh.py [--force]
        --force : relance dDDH et FastANI meme si leurs resultats existent
============================================================================
"""

import csv
import io
import itertools
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile

ICI = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ICI, "..", "..", "scripts")
GENOMES = os.path.join(ICI, "genomes")     # telecharges, non suivis par git
TRAVAIL = os.path.join(ICI, "work")        # sorties brutes de dddh.py / fastANI, non suivies par git

SEUIL_D = 0.0412      # meme seuil que scripts/dddh.py
SEUIL_ANI = 95.0
SEUIL_GGDC = 70.0
URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/download?include_annotation_type=GENOME_FASTA"


def lire_benchmark():
    groupes = {}
    with open(os.path.join(ICI, "benchmark.tsv")) as fh:
        for ligne in csv.DictReader(fh, delimiter="\t"):
            groupes.setdefault(ligne["group"], []).append(ligne)
    return groupes


def telecharger(groupes):
    for groupe, genomes in groupes.items():
        dossier = os.path.join(GENOMES, groupe)
        os.makedirs(dossier, exist_ok=True)
        for g in genomes:
            cible = os.path.join(dossier, g["name"] + ".fna")
            if os.path.exists(cible) and os.path.getsize(cible) > 0:
                continue
            print("  downloading %s (%s)" % (g["name"], g["accession"]))
            with urllib.request.urlopen(URL % g["accession"]) as reponse:
                archive = zipfile.ZipFile(io.BytesIO(reponse.read()))
            fastas = [n for n in archive.namelist() if n.endswith(".fna")]
            if not fastas:
                sys.exit("ERREUR : aucune sequence dans l'archive NCBI de %s" % g["accession"])
            with open(cible + ".part", "wb") as sortie:
                for nom in fastas:
                    sortie.write(archive.read(nom))
            os.replace(cible + ".part", cible)


def lancer(cmd, journal):
    with open(journal, "w") as fh:
        code = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if code != 0:
        sys.exit("ERREUR : la commande a echoue (voir %s) : %s" % (journal, " ".join(cmd)))


def calculer(groupes, force):
    for groupe in groupes:
        genomes = os.path.join(GENOMES, groupe)
        sortie = os.path.join(TRAVAIL, groupe)
        os.makedirs(sortie, exist_ok=True)
        if force or not os.path.exists(os.path.join(sortie, "dddh", "dddh_distance.tsv")):
            print("  dDDH   %s" % groupe)
            shutil.rmtree(os.path.join(sortie, "dddh"), ignore_errors=True)
            lancer([sys.executable, os.path.join(SCRIPTS, "dddh.py"), genomes, os.path.join(sortie, "dddh")],
                   os.path.join(sortie, "dddh.log"))
        if force or not os.path.exists(os.path.join(sortie, "ani.txt")):
            print("  FastANI %s" % groupe)
            liste = os.path.join(sortie, "liste.txt")
            with open(liste, "w") as fh:
                for f in sorted(os.listdir(genomes)):
                    if f.endswith(".fna"):
                        fh.write(os.path.join(genomes, f) + "\n")
            lancer(["fastANI", "--ql", liste, "--rl", liste, "-o", os.path.join(sortie, "ani.txt"), "-t", "4"],
                   os.path.join(sortie, "ani.log"))


def lire_matrice(chemin):
    with open(chemin) as fh:
        lignes = list(csv.reader(fh, delimiter="\t"))
    noms = lignes[0][1:]
    return {(l[0], n): (None if v in ("", "NA") else float(v)) for l in lignes[1:] for n, v in zip(noms, l[1:])}


def lire_ani(chemin):
    """ANI moyen des deux sens et plus petite fraction alignee, par paire non ordonnee."""
    sens = {}
    with open(chemin) as fh:
        for ligne in fh:
            c = ligne.rstrip("\n").split("\t")
            a, b = (os.path.splitext(os.path.basename(x))[0] for x in c[:2])
            sens[(a, b)] = (float(c[2]), float(c[3]) / float(c[4]))
    paires = {}
    for (a, b), (ani, fraction) in sens.items():
        if a < b and (b, a) in sens:
            paires[(a, b)] = ((ani + sens[(b, a)][0]) / 2, min(fraction, sens[(b, a)][1]))
    return paires


def espece(organisme):
    """Espece selon le nom NCBI : les deux premiers mots, sans les crochets des noms contestes."""
    return " ".join(organisme.replace("[", "").replace("]", "").split()[:2])


def lire_ggdc():
    """Valeurs GGDC saisies par l'utilisateur ; cree le modele a remplir s'il n'existe pas."""
    chemin = os.path.join(ICI, "ggdc_reference.tsv")
    valeurs = {}
    if os.path.exists(chemin):
        with open(chemin) as fh:
            for ligne in csv.DictReader(fh, delimiter="\t"):
                v = (ligne.get("ggdc_ddh_formula2") or "").strip().replace(",", ".").rstrip("%")
                if v:
                    valeurs[tuple(sorted((ligne["genome_a"], ligne["genome_b"])))] = float(v)
    return chemin, valeurs


def accord(lignes, cle_ref):
    """Tableau de contingence dDDH (pipeline) contre un critere de reference."""
    avec = [l for l in lignes if l["dddh_same"] is not None and l[cle_ref] is not None]
    vp = sum(1 for l in avec if l["dddh_same"] and l[cle_ref])
    vn = sum(1 for l in avec if not l["dddh_same"] and not l[cle_ref])
    fp = [l for l in avec if l["dddh_same"] and not l[cle_ref]]       # dit "meme espece" a tort
    fn = [l for l in avec if not l["dddh_same"] and l[cle_ref]]       # dit "especes differentes" a tort
    return len(avec), vp, vn, fp, fn


def spearman(xs, ys):
    def rangs(v):
        ordre = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(ordre):
            j = i
            while j + 1 < len(ordre) and v[ordre[j + 1]] == v[ordre[i]]:
                j += 1
            for k in range(i, j + 1):
                r[ordre[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = rangs(xs), rangs(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    force = "--force" in sys.argv
    groupes = lire_benchmark()
    print("1/3 genomes")
    telecharger(groupes)
    print("2/3 dDDH and ANI")
    calculer(groupes, force)
    print("3/3 comparison")
    chemin_ggdc, ggdc = lire_ggdc()

    lignes = []
    for groupe, genomes in groupes.items():
        d = lire_matrice(os.path.join(TRAVAIL, groupe, "dddh", "dddh_distance.tsv"))
        s = lire_matrice(os.path.join(TRAVAIL, groupe, "dddh", "dddh_similarite.tsv"))
        ani = lire_ani(os.path.join(TRAVAIL, groupe, "ani.txt"))
        organisme = {g["name"]: g["organism"] for g in genomes}
        for a, b in itertools.combinations(sorted(organisme), 2):
            dist = d.get((a, b))
            ani_ab = ani.get((a, b))                 # absent si FastANI n'aligne pas (< ~80 %)
            ref = ggdc.get((a, b))
            lignes.append({
                "group": groupe, "genome_a": a, "genome_b": b,
                "same_ncbi_species": espece(organisme[a]) == espece(organisme[b]),
                "dddh_distance": dist, "dddh_estimate": s.get((a, b)),
                "dddh_same": None if dist is None else dist <= SEUIL_D,
                "ani": ani_ab[0] if ani_ab else None,
                "aligned_fraction": ani_ab[1] if ani_ab else None,
                "ani_same": (ani_ab[0] >= SEUIL_ANI) if ani_ab else False,
                "ggdc_ddh_formula2": ref,
                "ggdc_same": None if ref is None else ref >= SEUIL_GGDC,
            })

    def fmt(v, decimales=None):
        if v is None:
            return "NA"
        if isinstance(v, bool):
            return "yes" if v else "no"
        return ("%." + str(decimales) + "f") % v if decimales is not None else str(v)

    colonnes = [("group", None), ("genome_a", None), ("genome_b", None), ("same_ncbi_species", None),
                ("dddh_distance", 4), ("dddh_estimate", 1), ("dddh_same", None), ("ani", 2),
                ("aligned_fraction", 2), ("ani_same", None), ("ggdc_ddh_formula2", 1), ("ggdc_same", None)]
    with open(os.path.join(ICI, "pairs.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([c for c, _ in colonnes])
        for l in sorted(lignes, key=lambda l: (l["dddh_distance"] is None, l["dddh_distance"])):
            w.writerow([fmt(l[c], dec) for c, dec in colonnes])

    if not os.path.exists(chemin_ggdc):
        with open(chemin_ggdc, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["genome_a", "genome_b", "ggdc_ddh_formula2"])
            for l in lignes:
                w.writerow([l["genome_a"], l["genome_b"], ""])

    # ---- resume ----
    out = ["# dDDH validation results", "",
           "Generated by `validate_dddh.py` from `benchmark.tsv` (%d genomes, %d pairs). Do not edit by hand." %
           (sum(len(g) for g in groupes.values()), len(lignes)), ""]

    def section(titre, cle_ref, note):
        n, vp, vn, fp, fn = accord(lignes, cle_ref)
        out.extend(["## " + titre, ""])
        if n == 0:
            out.extend([note, ""])
            return
        out.extend(["%d pairs compared, **%d agree (%.0f %%)**." % (n, vp + vn, 100.0 * (vp + vn) / n), "",
                    "| | reference: same species | reference: different species |", "|---|---|---|",
                    "| dDDH: same species (d <= %.4f) | %d | %d |" % (SEUIL_D, vp, len(fp)),
                    "| dDDH: different species | %d | %d |" % (len(fn), vn), ""])
        for nom, cas in (("dDDH says same species, reference says different", fp),
                         ("dDDH says different species, reference says same", fn)):
            if cas:
                out.append("%s:" % nom)
                out.append("")
                for l in cas:
                    out.append("- %s / %s: d = %s, ANI = %s%s" % (
                        l["genome_a"], l["genome_b"], fmt(l["dddh_distance"], 4), fmt(l["ani"], 2),
                        "" if l["ggdc_ddh_formula2"] is None else ", GGDC = %s" % fmt(l["ggdc_ddh_formula2"], 1)))
                out.append("")

    section("Against ANI (FastANI, >= 95 %)", "ani_same", "")
    section("Against GGDC (formula 2, >= 70 %)", "ggdc_same",
            "Not done yet: fill `ggdc_reference.tsv` with the values from https://ggdc.dsmz.de (see README.md) "
            "and run the script again.")
    avec = [l for l in lignes if l["dddh_distance"] is not None and l["ani"] is not None]
    if len(avec) > 2:
        out.extend(["## Rank correlation", "",
                    "Spearman correlation between the dDDH distance and ANI: %.3f over %d pairs "
                    "(expected close to -1)." % (spearman([l["dddh_distance"] for l in avec], [l["ani"] for l in avec]), len(avec))])
        avec_g = [l for l in lignes if l["dddh_distance"] is not None and l["ggdc_ddh_formula2"] is not None]
        if len(avec_g) > 2:
            out.append("Spearman correlation between the dDDH distance and GGDC's dDDH: %.3f over %d pairs." % (
                spearman([l["dddh_distance"] for l in avec_g], [l["ggdc_ddh_formula2"] for l in avec_g]), len(avec_g)))
        out.append("")
    proches = [l for l in lignes if l["dddh_distance"] is not None and SEUIL_D / 2 <= l["dddh_distance"] <= SEUIL_D * 2]
    out.extend(["## Pairs near the threshold", "",
                "%d pair(s) with a distance between %.4f and %.4f, where the decision is the most fragile:" %
                (len(proches), SEUIL_D / 2, SEUIL_D * 2), ""])
    for l in sorted(proches, key=lambda l: l["dddh_distance"]):
        out.append("- %s / %s: d = %s (%s), ANI = %s" % (l["genome_a"], l["genome_b"], fmt(l["dddh_distance"], 4),
                                                     "same" if l["dddh_same"] else "different", fmt(l["ani"], 2)))
    with open(os.path.join(ICI, "RESULTS.md"), "w") as fh:
        fh.write("\n".join(out).rstrip() + "\n")
    print("\n".join(out))
    print("\n-> pairs.tsv, RESULTS.md written in %s" % ICI)


if __name__ == "__main__":
    main()
