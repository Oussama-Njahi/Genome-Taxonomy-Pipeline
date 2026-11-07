#!/usr/bin/env python3
# arbre.py - Arbre phylogenomique, etape 'arbre' du pipeline (voir workflow/Snakefile).
#
# EasyCGTree4 (Perl) aligne les 120 genes marqueurs bacteriens (bac120 de GTDB) et FastTree construit l'arbre sur
# la supermatrice. Les feuilles reprennent ensuite les noms des genomes, et l'arbre non enracine de FastTree est
# enracine par root_tree.py, au point median ou sur le groupe externe donne avec --outgroup.
#
# Usage : python arbre.py <dossier_genomes> <dossier_resultats> <job_id> <threads> [--outgroup <genome>[,...]]
import argparse
import glob
import os
import re
import shutil
import sys

import pipeline_web as pw


def nom_easycgtree(nom):
    """EasyCGTree renomme ses fichiers d'entree avant de les lire (sub rename de EasyCGTree.pl) : espaces et '_'
    deviennent '-', et ce nom-la sert d'etiquette dans l'arbre. Les autres caracteres permis (NOM_PERMIS) restent."""
    return re.sub(r"[\s_]", "-", nom)


def cleanup_easycgtree(input_name):
    """Remove this job's working files from the shared EasyCGTree tool directory."""
    for path in glob.glob(pw.EASYCG + "/" + input_name + "*"):
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                os.remove(path)
            except OSError:
                pass


def fichier_arbre(input_name):
    """L'arbre d'EasyCGTree, <input>.bac120.<n>.supermatrix.fasttree.tree : n est le nombre de genes gardes dans la
    supermatrice, 120 seulement si chaque gene est present dans assez de genomes (-gene_cutoff 0.8, soit 4 genomes
    sur 5). Il en garde moins quand des genomes n'ont pas tous les genes (118 avec deux Klebsiella qui en ont 118).
    Renvoie (chemin, n), ou (None, None) s'il n'y a pas d'arbre."""
    for chemin in glob.glob(pw.EASYCG + "/" + glob.escape(input_name) + ".bac120.*.supermatrix.fasttree.tree"):
        n = chemin[len(pw.EASYCG + "/" + input_name + ".bac120."):].split(".")[0]
        if n.isdigit():
            return chemin, int(n)
    return None, None


def construire_arbre(fichiers, results_dir, job_id, threads, outgroup=None):
    noms = [os.path.splitext(os.path.basename(f))[0] for f in fichiers]
    vrai_nom = {}
    for nom in noms:
        autre = vrai_nom.setdefault(nom_easycgtree(nom), nom)
        if autre != nom:
            pw.fail("EasyCGTree turns '_' into '-', so '%s' and '%s' would get the same name in the tree: "
                    "rename one of them" % (autre, nom))
    if not os.path.exists(pw.EASYCG + "/EasyCGTree.pl"):
        pw.fail("EasyCGTree4 not found at %s - it must be installed manually "
                "(not available via conda), see the README's Installation section" % pw.EASYCG)
    input_name = "job_" + job_id
    entree = pw.EASYCG + "/" + input_name
    os.makedirs(entree, exist_ok=True)
    for g in fichiers:
        shutil.copy(g, entree)
    pw.run(["perl", "EasyCGTree.pl", "-input", input_name,
            "-hmm", "bac120", "-tree", "sm", "-tree_app", "fasttree",
            "-thread", str(threads)], cwd=pw.EASYCG)
    arbre, n_genes = fichier_arbre(input_name)
    if arbre is None:
        pw.fail("tree file not found (%s/%s.bac120.<genes>.supermatrix.fasttree.tree) - EasyCGTree may have failed"
                % (pw.EASYCG, input_name))
    print("    (tree built on %d of the 120 bac120 marker genes)" % n_genes)
    # les feuilles reprennent les noms des genomes (ceux des matrices), pas ceux qu'EasyCGTree leur a donnes
    with open(arbre) as fh:
        newick = fh.read()
    feuilles = re.findall(r"(?<=[(,])([^(),:;]+)(?=[:),;])", newick)
    inconnues = sorted(set(feuilles) - set(vrai_nom))
    if inconnues:
        pw.fail("unexpected leaf names in EasyCGTree's tree: %s" % ", ".join(inconnues))
    newick = re.sub(r"(?<=[(,])([^(),:;]+)(?=[:),;])", lambda m: vrai_nom[m.group(1)], newick)
    # FastTree ecrit un arbre non enracine : on le garde tel quel, et arbre.tree
    # (lu par les figures et l'interface) est sa version enracinee, au point
    # median ou sur le groupe externe donne avec --outgroup.
    brut = results_dir + "/arbre_non_enracine.tree"
    with open(brut, "w") as fh:
        fh.write(newick)
    cleanup_easycgtree(input_name)
    cmd = [sys.executable, pw.SCRIPTS + "/root_tree.py", brut, results_dir + "/arbre.tree"]
    if outgroup:
        cmd += ["--outgroup", outgroup]
    pw.run(cmd)
    print("    -> arbre.tree (rooted) and arbre_non_enracine.tree written to results")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the rooted phylogenomic tree of a genome folder.")
    parser.add_argument("genomes_dir")
    parser.add_argument("results_dir")
    parser.add_argument("job_id")
    parser.add_argument("threads", type=int)
    parser.add_argument("--outgroup", default=None)
    args = parser.parse_args()
    construire_arbre(sorted(glob.glob(args.genomes_dir + "/*.fna")), args.results_dir, args.job_id, args.threads,
                     args.outgroup)
