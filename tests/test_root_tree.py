"""Tests de scripts/root_tree.py. Lancer avec : python tests/test_root_tree.py  (ou pytest)."""
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from root_tree import distances, enraciner, feuilles_du_cote, lire_newick  # noqa: E402

# arbre non enracine tel que FastTree l'ecrit (trifurcation a la racine du fichier)
NON_ENRACINE = "(A:0.03,(B:0.001,C:0.001)0.99:0.02,(D:0.05,OUT:0.36)0.87:0.03);"


def patristiques(newick):
    """Distance entre chaque paire de feuilles : ne doit pas dependre de la racine."""
    arbre, _ = lire_newick(newick)
    par_nom = {nom: u for u, nom in arbre.noms.items()}
    return {(a, b): round(distances(arbre, par_nom[a])[0][par_nom[b]], 6)
            for a, b in itertools.combinations(sorted(par_nom), 2)}


def supports(newick):
    """Support de chaque bipartition (cote sans la feuille de reference, pour etre independant de la racine)."""
    arbre, _ = lire_newick(newick)
    reference = min(arbre.noms.values())
    trouves = {}
    for u in arbre.voisins:
        for v, (_, support) in arbre.voisins[u].items():
            cote = feuilles_du_cote(arbre, u, v)
            if support is not None and reference not in cote:
                trouves[frozenset(cote)] = support
    return trouves


def enfants_de_la_racine(newick):
    arbre, racine = lire_newick(newick)
    return [feuilles_du_cote(arbre, racine, v) for v in arbre.voisins[racine]]


def test_point_median_isole_la_longue_branche():
    resultat, methode = enraciner(NON_ENRACINE)
    assert methode == "midpoint"
    cotes = enfants_de_la_racine(resultat)
    assert len(cotes) == 2 and {"OUT"} in cotes
    # la racine est a egale distance des deux feuilles les plus eloignees (OUT et D exclue : OUT-A = 0.42)
    arbre, racine = lire_newick(resultat)
    d = distances(arbre, racine)[0]
    par_nom = {nom: u for u, nom in arbre.noms.items()}
    assert abs(d[par_nom["OUT"]] - d[par_nom["A"]]) < 1e-6


def test_distances_et_supports_conserves():
    for groupe in (None, ["OUT"], ["B", "C"], ["D", "OUT"]):
        resultat, _ = enraciner(NON_ENRACINE, groupe)
        assert patristiques(resultat) == patristiques(NON_ENRACINE), groupe
        assert supports(resultat) == supports(NON_ENRACINE), groupe


def test_groupe_externe():
    resultat, methode = enraciner(NON_ENRACINE, ["B", "C"])
    assert methode == "outgroup (B, C)"
    assert {"B", "C"} in enfants_de_la_racine(resultat)


def test_arbre_deja_enracine():
    # un arbre deja enracine ailleurs est re-enracine sans garder l'ancienne racine comme noeud de passage
    deja = "((A:0.03,(B:0.001,C:0.001)0.99:0.02):0.01,(D:0.05,OUT:0.36)0.87:0.02);"
    resultat, _ = enraciner(deja)
    assert patristiques(resultat) == patristiques(NON_ENRACINE)
    assert {"OUT"} in enfants_de_la_racine(resultat)


def test_erreurs():
    for groupe, message in ((["X"], "not found"), (["B", "OUT"], "not monophyletic"),
                            (["A", "B", "C", "D", "OUT"], "every genome")):
        try:
            enraciner(NON_ENRACINE, groupe)
        except ValueError as erreur:
            assert message in str(erreur)
        else:
            raise AssertionError("ValueError attendue pour %s" % groupe)


if __name__ == "__main__":
    for nom, fonction in sorted(globals().items()):
        if nom.startswith("test_"):
            fonction()
            print("ok  ", nom)
