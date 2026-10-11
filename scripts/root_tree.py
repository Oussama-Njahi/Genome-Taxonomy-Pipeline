#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
root_tree.py  --  Enracinement d'un arbre Newick
============================================================================
FastTree (via EasyCGTree) produit un arbre NON enracine : la "racine" du
fichier Newick est un noeud arbitraire a trois branches, et la lire comme une
vraie racine fait apparaitre de faux groupes freres (un groupe externe se
retrouve "frere" d'un des genomes etudies).

Deux methodes :
  - point median (par defaut) : la racine est placee au milieu du plus long
    chemin entre deux feuilles. Ne demande aucune connaissance prealable ;
    suppose des vitesses d'evolution a peu pres egales.
  - groupe externe (--outgroup) : la racine est placee au milieu de la branche
    qui separe le groupe externe du reste. C'est la methode a preferer quand
    on connait un groupe externe ; il doit former un groupe monophyletique.

Les valeurs de support (etiquettes des noeuds internes) appartiennent aux
BRANCHES, pas aux noeuds : elles sont stockees sur les branches pour rester
attachees a la bonne bipartition quand la racine change de place.

USAGE :
    python root_tree.py <arbre_entree> <arbre_sortie> [--outgroup NOM1,NOM2]
============================================================================
"""

import sys


class Arbre:
    """Arbre non oriente : voisins[u] = {v: [longueur, support]}, noms[u] pour les feuilles."""

    def __init__(self):
        self.voisins = {}
        self.noms = {}
        self.suivant = 0     # compteur : un numero supprime n'est jamais reutilise

    def noeud(self):
        u = self.suivant
        self.suivant += 1
        self.voisins[u] = {}
        return u

    def relier(self, u, v, longueur, support):
        self.voisins[u][v] = [longueur, support]
        self.voisins[v][u] = [longueur, support]

    def couper(self, u, v):
        del self.voisins[u][v]
        del self.voisins[v][u]

    def feuilles(self):
        return [u for u in self.voisins if u in self.noms]


def lire_newick(texte):
    """Lit un Newick (noms sans guillemets, longueurs et supports optionnels).
    Renvoie (arbre, racine_du_fichier)."""
    s = texte.strip().rstrip(";").strip()
    arbre = Arbre()
    pos = [0]

    def etiquette():
        debut = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ":,()":
            pos[0] += 1
        return s[debut:pos[0]].strip()

    def lire_noeud():
        u = arbre.noeud()
        enfants = []
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            while True:
                enfants.append(lire_noeud())
                if pos[0] < len(s) and s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if pos[0] < len(s) and s[pos[0]] == ")":
                    pos[0] += 1
                    break
                raise ValueError("malformed Newick tree near character %d" % pos[0])
        nom = etiquette()
        longueur = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            debut = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                pos[0] += 1
            longueur = float(s[debut:pos[0]])
        for v, l_v, etiq_v in enfants:
            # l'etiquette d'un noeud interne est le support de la branche au-dessus de lui
            arbre.relier(u, v, l_v, etiq_v if v not in arbre.noms else None)
        if not enfants:
            if not nom:
                raise ValueError("a leaf has no name")
            arbre.noms[u] = nom
        return u, longueur, (nom or None)

    racine, _, _ = lire_noeud()
    if len(arbre.feuilles()) < 2:
        raise ValueError("the tree has fewer than two leaves")
    return arbre, racine


def supprimer_noeud_de_passage(arbre, u):
    """Un noeud interne a 2 voisins n'est pas une vraie bifurcation (ancienne
    racine d'un arbre deja enracine) : on fusionne ses deux branches."""
    if u in arbre.noms or len(arbre.voisins[u]) != 2:
        return
    (a, (l_a, s_a)), (b, (l_b, s_b)) = arbre.voisins[u].items()
    arbre.couper(u, a)
    arbre.couper(u, b)
    del arbre.voisins[u]
    arbre.relier(a, b, l_a + l_b, s_a if s_a is not None else s_b)


def distances(arbre, depart):
    """Distance et predecesseur de chaque noeud depuis `depart`."""
    dist, pred = {depart: 0.0}, {depart: None}
    pile = [depart]
    while pile:
        u = pile.pop()
        for v, (longueur, _) in arbre.voisins[u].items():
            if v not in dist:
                dist[v] = dist[u] + longueur
                pred[v] = u
                pile.append(v)
    return dist, pred


def inserer_racine(arbre, u, v, distance_depuis_u):
    """Coupe la branche u-v et y place un nouveau noeud racine."""
    longueur, support = arbre.voisins[u][v]
    arbre.couper(u, v)
    racine = arbre.noeud()
    arbre.relier(racine, u, distance_depuis_u, support)
    arbre.relier(racine, v, longueur - distance_depuis_u, support)
    return racine


def racine_point_median(arbre):
    """Racine au milieu du plus long chemin feuille-feuille (diametre)."""
    feuilles = arbre.feuilles()
    dist, _ = distances(arbre, feuilles[0])
    a = max(feuilles, key=lambda f: dist[f])
    dist, pred = distances(arbre, a)
    b = max(feuilles, key=lambda f: dist[f])
    moitie = dist[b] / 2.0
    # on remonte de b vers a jusqu'a la branche qui contient le milieu
    v = b
    while pred[v] is not None and dist[pred[v]] > moitie:
        v = pred[v]
    u = pred[v]
    if u is None:                    # toutes les longueurs sont nulles
        u, v = v, next(iter(arbre.voisins[v]))
        return inserer_racine(arbre, u, v, 0.0)
    return inserer_racine(arbre, u, v, moitie - dist[u])


def feuilles_du_cote(arbre, u, v):
    """Noms des feuilles atteintes depuis v sans repasser par u."""
    noms, pile, vus = set(), [v], {u, v}
    while pile:
        x = pile.pop()
        if x in arbre.noms:
            noms.add(arbre.noms[x])
        for y in arbre.voisins[x]:
            if y not in vus:
                vus.add(y)
                pile.append(y)
    return noms


def racine_groupe_externe(arbre, groupe):
    """Racine au milieu de la branche qui isole exactement le groupe externe."""
    connus = set(arbre.noms.values())
    inconnus = sorted(set(groupe) - connus)
    if inconnus:
        raise ValueError("outgroup not found in the tree: %s (leaves: %s)"
                         % (", ".join(inconnus), ", ".join(sorted(connus))))
    if set(groupe) == connus:
        raise ValueError("the outgroup cannot contain every genome")
    for u in list(arbre.voisins):
        for v in list(arbre.voisins[u]):
            if feuilles_du_cote(arbre, u, v) == set(groupe):
                return inserer_racine(arbre, u, v, arbre.voisins[u][v][0] / 2.0)
    raise ValueError("the outgroup is not monophyletic in this tree: %s" % ", ".join(sorted(groupe)))


def ecrire_newick(arbre, racine):
    def nombre(x):
        return ("%.5f" % x).rstrip("0").rstrip(".") if x else "0"

    def ecrire(u, parent):
        if u in arbre.noms:
            return arbre.noms[u]
        enfants = [v for v in arbre.voisins[u] if v != parent]
        parties = []
        for v in enfants:
            longueur, support = arbre.voisins[u][v]
            etiq = support if (support is not None and v not in arbre.noms) else ""
            parties.append(ecrire(v, u) + etiq + ":" + nombre(longueur))
        return "(" + ",".join(parties) + ")"

    return ecrire(racine, None) + ";\n"


def enraciner(texte, groupe_externe=None):
    """Renvoie (newick_enracine, methode)."""
    arbre, racine_fichier = lire_newick(texte)
    supprimer_noeud_de_passage(arbre, racine_fichier)
    if groupe_externe:
        racine = racine_groupe_externe(arbre, groupe_externe)
        methode = "outgroup (%s)" % ", ".join(groupe_externe)
    else:
        racine = racine_point_median(arbre)
        methode = "midpoint"
    return ecrire_newick(arbre, racine), methode


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    groupe = None
    if "--outgroup" in sys.argv:
        valeur = sys.argv[sys.argv.index("--outgroup") + 1]
        args.remove(valeur)
        groupe = [g.strip() for g in valeur.split(",") if g.strip()]
    if len(args) != 2:
        print("Usage : python root_tree.py <arbre_entree> <arbre_sortie> [--outgroup NOM1,NOM2]")
        sys.exit(1)
    try:
        with open(args[0]) as fh:
            resultat, methode = enraciner(fh.read(), groupe)
    except ValueError as erreur:
        print("ERREUR : %s" % erreur)
        sys.exit(1)
    with open(args[1], "w") as fh:
        fh.write(resultat)
    print("Arbre enracine (%s) ecrit dans : %s" % (methode, args[1]))
