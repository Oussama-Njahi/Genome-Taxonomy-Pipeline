# TGA — Taxo Genomic Analyses · Guide d'usage du logo

## 1. Le logo
- **Idée** : les lettres T, G et A sont les feuilles d'un dendrogramme, la figure même que trace le pipeline.
  G et A sont sœurs (réunies en bas), T les rejoint plus haut.
- **Versions** :
  - logotype seul (`tga-wordmark-*`), la version principale ;
  - empilée avec le nom (`tga-stacked-*`) ;
  - horizontale avec le nom sur trois lignes (`tga-horizontal-*`), pour les en-têtes ;
  - symbole (`tga-symbol-color.svg`), de 32 à 64 px ;
  - favicon (`tga-favicon-*`, `web/`), à 32 px et moins.
- **Fichiers** : `final/*.svg` (masters vectoriels), `final/*-1200.png` (écran), `final/web/` (favicon, icônes, manifeste).
  Les lettres sont des contours vectoriels, pas du texte : ne pas les retaper.

## 2. Zone de protection
Laisser autour du logo un espace libre d'au moins **1 × x** de chaque côté. **x** est l'écart vertical entre les deux
barres horizontales de l'arbre. La zone grandit avec le logo : ce n'est jamais une distance fixe.

## 3. Taille minimale
| Version | Écran | Impression |
|---|---|---|
| Logotype, empilée, horizontale | 120 px de large | 30 mm de large |
| Symbole (`tga-symbol-color.svg`) | 32 px | 10 mm |
| Favicon (`tga-favicon-*`) | 16 px | — |

En dessous de 120 px, les déliés du A et du G disparaissent : il faut passer au symbole, puis au favicon.

## 4. Couleurs
| Nom | Usage | HEX | RGB | CMJN (approx.) |
|---|---|---|---|---|
| Encre | lettres | `#2A211B` | 42 33 27 | 0 21 36 84 |
| Bleu de Prusse | branches | `#2F4A63` | 47 74 99 | 53 25 0 61 |
| Papier | fond, tuile du favicon | `#F6EEE3` | 246 238 227 | 0 3 8 4 |

Ce sont les encres de l'interface et des figures R (`plot_all_web.R`, `plates.js`). Les valeurs CMJN sont des
conversions directes : il faut les faire valider par l'imprimeur, et aucun Pantone n'est encore fixé.

**Associations autorisées :**
- couleur sur papier ou sur blanc ;
- version inversée (`tga-wordmark-reversed.svg`) sur encre ;
- noir sur blanc ;
- blanc sur un fond sombre uni ;
- une couleur, en bleu (`-mono-2f4a63`) ou à l'encre (`-mono-2a211b`).

## 5. Typographie
- **TGA** : Bodoni Moda, graisse 600, taille optique 36 (le symbole utilise 800 / 6, plus robuste).
- **Nom « TAXO GENOMIC ANALYSES »** : Spectral SemiBold, en capitales espacées.
- **Licences** : SIL Open Font License pour les deux polices (`webapp/frontend/fonts/OFL-*.txt`), ce qui autorise
  leur usage dans un logo.

## 6. À ne pas faire
- Déformer, faire pivoter ou recolorer le logo hors palette.
- Ajouter une ombre, un contour ou un dégradé.
- Déplacer les branches ou changer la topologie : T reste seul, G et A restent sœurs.
- Retaper « TGA » dans une police.
- Utiliser le logotype sous 120 px.
- Le poser sur un fond chargé sans tuile papier.

## 7. Masters
Tous les fichiers sont régénérés par `python3 final.py`, depuis `branding/tga/`, à partir des polices de
`webapp/frontend/fonts/`. Les exports (noir, blanc, mono, PNG, icônes web) viennent de `export_variants.py` du skill
logo-design. La planche de présentation est dans `board/tga-board.html` et `board/slides/`.

**Point ouvert :** l'interface web affiche encore « Taxo Genomics Analyses » (avec un s). Le logo utilise
« Taxo Genomic Analyses ».
