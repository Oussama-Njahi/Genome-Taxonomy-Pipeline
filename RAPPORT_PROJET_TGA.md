# TGA — Taxo Genomics Analyses

**Pipeline intégré de taxonomie génomique bactérienne**

Stage de recherche · Institut des Régions Arides (IRA), Médenine, Tunisie · Laboratoire de Recherche « Écologie Pastorale »
Auteur : Oussama Njahi · Dépôt : https://github.com/Oussama-Njahi/Genome-Taxonomy-Pipeline
État du document : 14 octobre 2026, version `main` = `b5205b0` (37 commits), plus la validation non encore commitée (§ 10.2).

---

## Sommaire

1. [Résumé](#1-résumé)
2. [Contexte et objectif](#2-contexte-et-objectif)
3. [Chronologie du projet](#3-chronologie-du-projet)
4. [Architecture](#4-architecture)
5. [Les neuf étapes du pipeline](#5-les-neuf-étapes-du-pipeline)
6. [Méthodes implémentées dans le projet](#6-méthodes-implémentées-dans-le-projet)
7. [Interface web](#7-interface-web)
8. [Identité visuelle : le logo TGA](#8-identité-visuelle--le-logo-tga)
9. [Résultats sur le jeu de test](#9-résultats-sur-le-jeu-de-test)
10. [Validation](#10-validation)
11. [Installation et utilisation](#11-installation-et-utilisation)
12. [Organisation du dépôt](#12-organisation-du-dépôt)
13. [Limites connues et travaux restants](#13-limites-connues-et-travaux-restants)
14. [Références](#14-références)

---

## 1. Résumé

TGA prend en entrée un dossier d'assemblages de génomes bactériens (FASTA, au moins 5). En une seule exécution, il
calcule tous les critères génomiques utilisés pour délimiter les espèces et les genres :

- **qualité des génomes** : taille, %GC, contigs ; complétude et contamination avec CheckM2 ;
- **délimitation d'espèce** : ANI (FastANI) et dDDH (implémentation de la formule 2 du GBDP) ;
- **délimitation de genre** : AAI (EzAAI) et POCP (implémentation sur DIAMOND) ;
- **phylogénomique** : arbre sur 120 gènes marqueurs bactériens (bac120, EasyCGTree + FastTree), enraciné ;
- **figures** : heatmaps et arbre, à 300 dpi.

Le pipeline s'utilise en ligne de commande ou par une interface web locale. Dans l'interface, on dépose les génomes, on
choisit les étapes, on suit la progression, puis on lit les résultats sous forme de planches interactives.

Ce qui distingue l'outil : tout tourne **en local et hors ligne**, les quatre critères sont réunis dans un seul run,
et chaque étape est optionnelle.

Sur le jeu de validation, les verdicts d'espèce du dDDH concordent avec ceux de l'ANI pour **37 paires sur 37**
(corrélation de rang −0,986).

---

## 2. Contexte et objectif

Décrire ou reclasser une espèce bactérienne demande aujourd'hui des critères génomiques :

| Critère | Niveau | Seuil usuel | Référence |
|---|---|---|---|
| ANI (Average Nucleotide Identity) | espèce | ≥ 95–96 % | Jain et al. 2018 |
| dDDH (digital DNA–DNA hybridization) | espèce | ≥ 70 % (distance d ≤ 0,0412) | Auch et al. 2010 |
| AAI (Average Amino acid Identity) | genre | ≈ 65 % (indicatif) | Kim et al. 2021 |
| POCP (Percentage of Conserved Proteins) | genre | ≥ 50 % | Qin et al. 2014 |
| Phylogénomie sur gènes conservés | toutes | topologie | Parks et al. 2018 (bac120) |

Ces calculs se font normalement avec des outils séparés : des services web (GGDC, TYGS, JSpeciesWS) et plusieurs
programmes, chacun avec son format. **L'objectif du stage** était d'en faire un pipeline unique, reproductible et
utilisable par un biologiste sans ligne de commande, qui produise des résultats directement lisibles.

---

## 3. Chronologie du projet

| Période | Étape |
|---|---|
| **7 – 14 sept.** | Premier pipeline en ligne de commande : QC, ANI, AAI, POCP, arbre phylogénomique. Option `--steps` pour ne lancer que certaines étapes. |
| **17 – 20 sept.** | Publication de l'interface web avec sélecteur d'étapes. **Ajout du dDDH** (`dddh.py`). CheckM2 devient une étape indépendante et optionnelle. Annulation d'un run depuis l'interface. Corrections pour qu'un clone neuf du dépôt fonctionne (chemins en dur, dépendance BLAST manquante). Message clair si EasyCGTree4 n'est pas installé. |
| **29 sept. – 2 oct.** | Logigramme du pipeline (`tga-pipeline-flowchart.png`). Refonte visuelle de l'interface. Fiabilité : les figures et onglets des étapes non lancées sont ignorés au lieu de faire échouer le run ; une paire dDDH sans alignement est notée NA, jamais « identique » ; `pipeline.py` lance exactement le même pipeline que l'interface ; ajout de la dépendance `r-reshape2`. |
| **5 – 6 oct.** | Documentation de la méthode dDDH implémentée (Auch 2010 contre le GGDC actuel). Progression numérotée parmi les étapes choisies. Refonte de l'interface et des figures R en « planche d'atlas ancien ». Test complet sur 5 génomes NCBI. **Logo TGA**. Colonne « Proteins » dans le registre. Jobs terminés retrouvés après un redémarrage du serveur. |
| **11 oct.** | **Corrections scientifiques** : arbre enraciné ; dDDH avec BLASTN classique ; POCP avec DIAMOND ultra-sensible ; fraction alignée de l'ANI. Interface : vocabulaire « step », nouvelle barre de progression. |
| **11 – 14 oct.** | **Validation** : dDDH sur 37 paires de génomes de référence ; POCP DIAMOND contre BLASTP. *(non commitée)* |

---

## 4. Architecture

```
                 ┌──────────────────────────────┐
  navigateur ──► │ webapp/frontend (HTML/CSS/JS)│
                 └──────────────┬───────────────┘
                                │ HTTP /api/...
                 ┌──────────────▼───────────────┐
                 │ webapp/backend (FastAPI)     │  file d'attente, un job à la fois
                 │ app.py · job_manager.py      │  suit la progression sur stdout
                 └──────────────┬───────────────┘
  ligne de                      │ subprocess
  commande ──► scripts/pipeline.py ──► scripts/pipeline_web.py   (les 9 étapes, un seul exemplaire)
                                │
        ┌───────────┬───────────┼────────────┬──────────────┬─────────────┐
     Python       BLAST+     DIAMOND       FastANI        EzAAI (Java)   Perl : EasyCGTree4 → FastTree
   (QC, dDDH,   (dDDH)       (POCP)        (ANI)          (AAI)          R : ggplot2, ggtree (figures)
    POCP, racine)                                         CheckM2 (env. conda séparé)
```

- **Python** orchestre tout par `subprocess`. Il appelle **Perl** (EasyCGTree) et **R** (figures), sur le même modèle
  multi-langage que des pipelines de référence comme PGCGAP, qui ne calcule ni AAI ni POCP.
- **Un seul pipeline** : la ligne de commande et l'interface web lancent le même `pipeline_web.py`. Un résultat obtenu
  d'une façon est donc reproductible de l'autre.
- **Progression** : chaque étape écrit `PROGRESS n/total <nom>` sur la sortie standard, que le backend lit pour
  piloter la barre de progression.

### Environnement logiciel (versions utilisées)

| Outil | Version | Rôle |
|---|---|---|
| Python | 3.14.7 | orchestration, QC, dDDH, POCP, enracinement |
| FastANI | 1.34 | ANI |
| BLAST+ (`blastn`) | via bioconda | alignements du dDDH |
| DIAMOND | 2.2.8 | alignements du POCP |
| Prodigal | 2.6.3 | prédiction des protéines |
| EzAAI | 1.2.4 | AAI |
| CheckM2 | 1.1.0 (env. `checkm2env`) | complétude et contamination |
| EasyCGTree | 4.2 | extraction et alignement des 120 gènes bac120 |
| FastTree | 2.2.0 | arbre de maximum de vraisemblance |
| R | 4.5.3 (ggplot2, ggtree, ape, reshape2) | figures |
| FastAPI + uvicorn | via conda-forge | serveur web |

Machine de développement : Ubuntu sous WSL2, 12 cœurs, 11 Go de RAM.

---

## 5. Les neuf étapes du pipeline

Clés `--steps` : `qc, checkm2, proteines, ani, dddh, aai, pocp, arbre, figures`. Dépendances : `checkm2` a besoin de
`qc`, et `pocp` a besoin de `proteines`. Si on saute une étape dont une autre dépend, seule celle-ci est ignorée,
avec un avertissement ; le reste tourne.

| # | Étape | Outil et paramètres | Sorties (`results/`) |
|---|---|---|---|
| 1 | **Contrôle qualité** | script Python : taille, %GC, nombre de contigs | `qc.txt` |
| 2 | **Complétude / contamination** | CheckM2 `predict`, 4 threads, environnement `checkm2env` | colonnes ajoutées à `qc.txt` |
| 3 | **Prédiction des protéines** | Prodigal `-q` | `proteins/<génome>.faa` |
| 4 | **ANI** | FastANI tous contre tous (`--ql`/`--rl`), `--threads` | `ani_resultats.txt` (ANI, fragments alignés, fragments totaux) |
| 5 | **dDDH** | `dddh.py` : BLASTN `-task blastn`, e-value ≤ 1e-2, formule 2, trimming (§ 6.1) | `dddh_out/dddh_distance.tsv`, `dddh_similarite.tsv` |
| 6 | **AAI** | EzAAI `extract` puis `calculate` | `aai_resultats.tsv` |
| 7 | **POCP** | `pocp.py` : DIAMOND `blastp --ultra-sensitive` (§ 6.2) | `pocp_out/pocp_matrice.tsv` |
| 8 | **Arbre phylogénomique** | EasyCGTree `-hmm bac120 -tree sm -tree_app fasttree`, puis `root_tree.py` (§ 6.3) | `arbre.tree` (enraciné), `arbre_non_enracine.tree` |
| 9 | **Figures** | `plot_all_web.R` | `heatmap_ani.png`, `heatmap_dddh.png`, `heatmap_aai.png`, `heatmap_pocp.png`, `tree.png` |

**Figures (étape 9).** Les heatmaps utilisent une échelle divergente coupée exactement au seuil : bleu de Prusse du
côté « même taxon », terre d'ombre de l'autre. Les génomes sont ordonnés selon l'arbre, ou par regroupement UPGMA à
défaut d'arbre. Les valeurs « même taxon » sont en gras et une case sans valeur est marquée NA. Les paires que FastANI
ne peut pas aligner (< ~80 %) apparaissent en NA. Un astérisque marque un ANI ≥ 95 % mesuré sur moins de la moitié du
génome. L'arbre est dessiné « gravé », avec une accolade par genre et une barre d'échelle.

---

## 6. Méthodes implémentées dans le projet

### 6.1 dDDH : formule 2 du GBDP (`scripts/dddh.py`)

D'après Auch, Klenk & Göker (2010), Table 4, ligne « NCBI BLAST / Trimming (2) ».

1. **Alignement** : BLASTN de A contre B, puis de B contre A, en conservant les HSP avec une e-value ≤ 1e-2.
2. **Greedy-with-trimming** : les chevauchements entre HSP sur une même séquence sont corrigés. Les HSP les plus
   longues sont prises en premier, et des suivantes on ne compte que la partie pas encore couverte, pour qu'une région
   dupliquée ne soit comptée qu'une fois.
3. **Distance** : `d = 1 − (identités dans les HSP, deux sens) / (longueur des HSP après trimming, deux sens)`.
4. **Conversion et seuil** : `DDH ≈ 90,3998 − 438,3134 × d`, borné entre 0 et 100 %. Les espèces sont différentes si
   `d > 0,0412`, ce qui correspond à environ 70 % de DDH.

**Seuil et conversion sont deux résultats distincts** du tableau 4 d'Auch et al. (2010) :
- **le seuil 0,0412** est la distance qui classe le mieux les paires autour de 70 % de DDH (taux d'erreur 4,3 %) ;
- **la droite** vient d'une régression robuste sur toutes les paires.

Ils ne coïncident pas exactement : à d = 0,0412, la droite donne 72,3 %, et 70 % correspondrait à d = 0,0465.
L'article recommande d'utiliser **le seuil de distance** pour délimiter les espèces ; c'est lui qui décide dans le
pipeline, et le pourcentage n'est qu'une estimation. La droite n'est calibrée que près du seuil : sans borne, elle
donnerait jusqu'à −27 % entre espèces éloignées d'un même genre.

**Les autres réglages suivent aussi l'article** : BLASTN classique avec des mots de 11 nt (tableau 2), masquage des
régions peu complexes pendant la seule recherche des amorces (défaut de BLAST+), et HSP d'e-value > 1e-2 écartées.

**Choix de `-task blastn`.** Appelée sans `-task`, la commande `blastn` de BLAST+ lance megablast (mots de 28 nt), qui
n'aligne que les régions les plus conservées. Mesuré sur la paire *B. subtilis* 168 / *B. velezensis* FZB42 :

| | megablast (avant) | blastn (maintenant) |
|---|---|---|
| HSP trouvés | 934 | 2 952 |
| ADN aligné | 2,77 Mb | 4,14 Mb |
| Identité moyenne | 82,9 % | 78,3 % |
| Temps par sens | 3 s | 15 s |

Avec megablast, *S. aureus* sortait plus proche des *Bacillus* (d ≈ 0,14) que les *Bacillus* entre eux (d ≈ 0,21),
ce qui est biologiquement faux. Avec blastn, il est à 0,30 contre 0,23–0,27.

**Limites connues.**
- **Approximation du trimming** : on suppose une identité uniforme le long de chaque HSP tronquée.
- **Conversion différente de GGDC** : le service GGDC actuel (Meier-Kolthoff et al. 2013) convertit la distance avec
  un modèle logistique non publié. Les pourcentages ne sont donc pas comparables ; seules les décisions et le
  classement le sont. C'est pourquoi les figures montrent la distance.
- **Paire sans alignement** : elle est notée NA, jamais « identique ».

### 6.2 POCP (`scripts/pocp.py`)

D'après Qin et al. (2014) : `POCP = (C(A→B) + C(B→A)) / (T(A) + T(B)) × 100`, où C est le nombre de protéines
conservées et T le nombre total de protéines. Une protéine est conservée si son meilleur hit dans l'autre protéome a
une **e-value ≤ 1e-5**, une **identité ≥ 40 %** et un **alignement ≥ 50 % de sa longueur**. Deux génomes avec un
POCP > 50 % sont du même genre.

**Choix de DIAMOND `--ultra-sensitive`.** Qin et al. utilisaient BLASTP. Mesuré avec exactement les mêmes critères
(`validation/pocp/compare_engines.py`) :

| Paire | DIAMOND rapide | DIAMOND ultra-sensible | BLASTP |
|---|---|---|---|
| *B. subtilis* 168 / *B. licheniformis* | 73,11 % | 74,09 % | 74,18 % |
| *B. subtilis* 168 / *S. aureus* | 31,70 % | 32,55 % | 32,65 % |

Le mode ultra-sensible reproduit BLASTP à 0,1 point près, environ deux fois plus vite.

### 6.3 Enracinement de l'arbre (`scripts/root_tree.py`)

FastTree écrit un arbre **non enraciné**. Lire le premier nœud du fichier comme une racine faisait apparaître un
groupe externe comme frère d'un génome étudié : *S. aureus* sortait frère de *B. licheniformis*.

- **Par défaut, point médian** : la racine est placée au milieu du plus long chemin entre deux feuilles.
- **Avec `--outgroup`** : la racine est placée au milieu de la branche qui isole le groupe externe. Le module vérifie
  que ce groupe est monophylétique.
- **Supports conservés** : les valeurs de support sont stockées sur les branches, pas sur les nœuds, pour rester sur la
  bonne bipartition quand la racine change de place.
- **Arbre brut gardé** : l'arbre non enraciné reste disponible.
- **Tests** : 5 tests automatiques (`tests/test_root_tree.py`) vérifient que les distances entre feuilles et les
  supports sont conservés, couvrent le groupe externe et l'arbre déjà enraciné, et testent les messages d'erreur.

### 6.4 Fraction alignée de l'ANI

FastANI indique combien de fragments du génome requête ont été alignés. Un ANI n'est fiable que s'il porte sur une
part suffisante du génome. L'interface affiche cette fraction sous chaque ANI. Un ANI ≥ 95 % mesuré sur moins de 50 %
du génome passe en rouge (« too little to trust ») et reçoit un astérisque sur la figure.

---

## 7. Interface web

### 7.1 Backend (`webapp/backend`, FastAPI)

| Route | Rôle |
|---|---|
| `POST /api/jobs` | dépôt des génomes (`.fas`, `.fasta`, `.fna`, au moins 5) et choix des étapes ; crée un job |
| `GET /api/jobs/{id}/status` | état, étape en cours, journal |
| `POST /api/jobs/{id}/cancel` | annule le job, en arrêtant tout le groupe de processus (FastANI, CheckM2…) |
| `GET /api/jobs/{id}/results` | registre qualité (avec le nombre de protéines prédites) et liste des fichiers produits |
| `GET /api/jobs/{id}/files/{clé}` | téléchargement d'une matrice, de l'arbre ou d'une figure |

- **File d'attente** : les jobs passent un par un. Chaque job a son dossier `webapp/jobs/<id>/`, avec `genomes/`,
  `proteins/` et `results/`.
- **Après un redémarrage** : un job terminé est retrouvé sur le disque, ses résultats restent consultables.

### 7.2 Frontend (`webapp/frontend`, HTML/CSS/JS sans build)

Le site est conçu comme un **atlas d'histoire naturelle ouvert** (dans l'esprit des planches de Haeckel) :

- **Mise en page** : doubles pages avec pli central, papier, encre brune, bleu de Prusse et carmin.
- **Polices** : Bodoni Moda, Spectral et IBM Plex Mono, servies en local sous licence SIL OFL.
- **Animations** : GSAP.

Le parcours suit trois parties :

1. **The specimens** : dépôt des génomes par glisser-déposer, avec contrôle des extensions et du nombre.
2. **Contents** : les 9 étapes en cases à cocher (« Each step runs on its own »), avec les dépendances en note.
3. **In the press** : la progression. La barre est une règle graduée avec une division par étape choisie ; les étapes
   finies sont en bleu et l'étape en cours a des hachures carmin animées. Le pourcentage est en grands chiffres, la
   liste des étapes indique leur état et le journal (« Printer's log ») est consultable.
4. **Plates** : les résultats.
   - **Registre des spécimens** : taille, %GC, contigs, protéines, complétude, contamination, avec des remarques
     automatiques (complétude < 90 %, contamination > 5 %).
   - **Planches I–IV (ANI, dDDH, AAI, POCP)** : une matrice interactive avec l'arbre en marge. Survoler une case
     affiche la paire dans toutes les mesures, avec la fraction alignée pour l'ANI. Une liste montre les paires
     proches du seuil.
   - **Planche V** : l'arbre enraciné avec ses supports et les genres en accolades.
   - **Téléchargements** : chaque fichier brut et chaque figure à 300 dpi.

Le site s'adapte aux mobiles : les pages s'empilent et les tableaux défilent.

---

## 8. Identité visuelle : le logo TGA

**Idée** : les lettres **T, G, A** sont les feuilles d'un **dendrogramme**, la figure même que trace le pipeline. G et
A sont sœurs ; T les rejoint plus haut.

- **Construction** : les lettres sont tirées de Bodoni Moda (licence OFL) et converties en contours vectoriels. Il n'y
  a donc aucune police à installer.
- **Couleurs** : arbre en bleu de Prusse `#2F4A63`, lettres à l'encre `#2A211B`, sur papier `#F6EEE3`.
- **Versions** : logotype, empilée, horizontale, inversée, symbole 32–64 px, et favicon. Le favicon montre l'arbre et
  trois blocs représentant les feuilles, et reste lisible à 16 px.
- **Kit** (`branding/tga/`) : masters SVG, variantes noir, blanc et une couleur, PNG, favicons et icônes web, planche
  de présentation, guide d'usage (`GUIDELINES.md`). Tout se régénère avec `python3 final.py`.
- **Intégration** : logotype en frontispice de la page de titre, arbre dans la barre du haut, favicon et manifeste web.

---

## 9. Résultats sur le jeu de test

### 9.1 Cinq génomes complets de NCBI

| Génome | Accession | Rôle |
|---|---|---|
| *Bacillus subtilis* 168 | GCF_000009045.1 | référence |
| *Bacillus subtilis* NCIB 3610 | GCF_002055965.1 | même espèce que 168 |
| *Bacillus velezensis* FZB42 | GCF_000015785.2 | même genre, autre espèce |
| *Bacillus licheniformis* ATCC 14580 | GCF_000011645.1 | même genre, plus éloigné |
| *Staphylococcus aureus* NCTC 8325 | GCF_000013425.1 | groupe externe |

**Registre qualité**

| Génome | Taille (pb) | GC % | Contigs | Protéines | Complétude % | Contamination % |
|---|---|---|---|---|---|---|
| *B. licheniformis* ATCC 14580 | 4 222 597 | 46,2 | 1 | 4 295 | 99,99 | 0,15 |
| *B. subtilis* 168 | 4 215 606 | 43,5 | 1 | 4 226 | 100,0 | 0,05 |
| *B. subtilis* NCIB 3610 | 4 299 822 | 43,3 | 2 | 4 332 | 100,0 | 0,12 |
| *B. velezensis* FZB42 | 3 918 596 | 46,5 | 1 | 3 736 | 100,0 | 0,04 |
| *S. aureus* NCTC 8325 | 2 821 361 | 32,9 | 1 | 2 631 | 100,0 | 0,21 |

**Mesures par paire** (dDDH et POCP avec les réglages corrigés du 11 octobre ; fraction alignée la plus faible des
deux sens)

| Paire | ANI % (fraction) | dDDH d | AAI % | POCP % | Verdict |
|---|---|---|---|---|---|
| *B. subtilis* 168 / NCIB 3610 | 100,0 (98–100 %) | 0,0001 | 99,99 | 98,64 | même espèce |
| *B. subtilis* 168 / *B. velezensis* | 81,1 (47 %) | 0,2319 | 81,07 | 82,79 | autre espèce, même genre |
| *B. subtilis* 168 / *B. licheniformis* | 78,9 (23 %) | 0,2683 | 73,74 | 74,09 | autre espèce, même genre |
| *B. velezensis* / *B. licheniformis* | 78,9 (24 %) | 0,2686 | 73,05 | 72,49 | autre espèce, même genre |
| *Bacillus* / *S. aureus* | non aligné (< 80 %) | 0,30 – 0,31 | ≈ 57 | 32 – 34 | autre genre |

Les quatre critères et l'arbre enraciné concordent : *S. aureus* est le groupe externe, et les quatre *Bacillus*
forment un groupe soutenu (support 1,000).

### 9.2 Jeu de 14 génomes (étude de cas du stage)

Sur 14 génomes (*Bhargavaea*, *Candidatus Kurthia*, *Caryophanon* et quatre groupes externes), les quatre analyses
indépendantes montrent que ***Bacillus ndiopicus*** est génomiquement plus proche de *Caryophanon* (AAI ≈ 69 %, POCP ≈
50–54 %) que de *Bacillus subtilis* (AAI ≈ 60 %, POCP ≈ 38 %). Cela suggère une **erreur de classification**,
confirmée par l'ANI, l'AAI, le POCP et la topologie de l'arbre.

> Ces valeurs ont été obtenues avant les corrections du 11 octobre (dDDH, POCP, enracinement). Il faut relancer ce
> jeu avec la version actuelle avant de les publier.

---

## 10. Validation

### 10.1 POCP : DIAMOND contre BLASTP
Voir § 6.2 : l'écart est de 0,1 point avec BLASTP en mode ultra-sensible.

### 10.2 dDDH sur un jeu de référence (`validation/dddh/`)

- **Jeu de référence** : 19 génomes complets RefSeq en 4 groupes, soit **37 paires**, choisies pour encadrer le seuil
  d'espèce.
  - *Bacillus subtilis* et apparentés : *B. subtilis* ×2, *B. spizizenii*, *B. velezensis*, *B. amyloliquefaciens*,
    *B. licheniformis* ;
  - *Escherichia* / *Shigella* : *E. coli* K-12, *E. coli* O157:H7, *S. flexneri*, *E. fergusonii*, *E. albertii* ;
  - *Bacillus cereus* et apparentés : *B. cereus*, *B. anthracis*, « *B. thuringiensis* » 97-27, *B. mycoides* ;
  - *Klebsiella* : *K. pneumoniae*, *K. variicola*, *K. quasipneumoniae*, *K. oxytoca*.
- **Procédure** : `validate_dddh.py` télécharge les génomes, lance `dddh.py` tel quel et FastANI, puis compare les
  décisions. La durée est de 18 minutes.

**Accord avec l'ANI (≥ 95 %)** : **37 / 37 (100 %)**. La corrélation de rang entre la distance dDDH et l'ANI est de −0,986.

| | ANI : même espèce | ANI : espèces différentes |
|---|---|---|
| **dDDH : même espèce** | 5 | 0 |
| **dDDH : espèces différentes** | 0 | 32 |

Les cas difficiles sont bien traités :
- ***Shigella flexneri*** est classée avec *E. coli* (d = 0,022), ce qui est cohérent avec la génomique malgré le nom ;
- la souche nommée ***B. thuringiensis* 97-27** est classée avec *B. anthracis* (d = 0,026) ;
- les paires d'espèces différentes juste sous le seuil ANI sont bien séparées : *K. pneumoniae* / *K. variicola*
  (ANI 94,8 %, d = 0,058) et *B. velezensis* / *B. amyloliquefaciens* (ANI 93,9 %, d = 0,063).

**À compléter** : la comparaison avec **GGDC**, l'implémentation de référence. Le plan prévoit 8 soumissions couvrant
19 paires, dont les 10 plus proches du seuil (`ggdc_jobs.md`). Il faut reporter les valeurs dans
`ggdc_reference.tsv`, puis relancer le script.

---

## 11. Installation et utilisation

```bash
# environnement principal
mamba env create -f stage_environment.yml
conda activate stage
# CheckM2, dans un environnement séparé (conflit TensorFlow/Python)
mamba create -n checkm2env -y -c conda-forge -c bioconda python=3.12 checkm2
# EasyCGTree4 : installation manuelle dans EasyCGTree4/EasyCGTree.v4.2-Linux/ (voir README)
```

**Ligne de commande** : placer les génomes (`.fas`, `.fasta`, `.fna`, au moins 5) dans `genomes/`, puis lancer :

```bash
python scripts/pipeline.py                                   # les 9 étapes
python scripts/pipeline.py --steps qc,ani,dddh               # certaines étapes seulement
python scripts/pipeline.py --threads 8                       # threads pour FastANI et EasyCGTree
python scripts/pipeline.py --outgroup Staphylococcus-aureus-NCTC8325   # racine de l'arbre
```

**Interface web** :

```bash
cd webapp/backend && uvicorn app:app --host 127.0.0.1 --port 8000
# puis ouvrir http://127.0.0.1:8000
```

**Tests et validation** :

```bash
python tests/test_root_tree.py
python validation/dddh/validate_dddh.py
python validation/pocp/compare_engines.py A.faa B.faa
```

Pour 5 génomes d'environ 4 Mb, avec les réglages actuels : dDDH ≈ 5 min, POCP ≈ 2 min, arbre ≈ 1,5 min. CheckM2
demande environ 6 Go de RAM.

---

## 12. Organisation du dépôt

```
.
├── scripts/
│   ├── pipeline.py          # point d'entrée en ligne de commande
│   ├── pipeline_web.py      # les 9 étapes (aussi utilisé par l'interface)
│   ├── dddh.py              # dDDH, formule 2 du GBDP
│   ├── pocp.py              # POCP sur DIAMOND
│   ├── root_tree.py         # enracinement de l'arbre
│   └── plot_all_web.R       # figures
├── webapp/
│   ├── backend/             # app.py (API), job_manager.py (file d'attente)
│   └── frontend/            # index.html, css/, js/ (app.js, plates.js), fonts/, icons/, img/, vendor/ (GSAP)
├── tests/test_root_tree.py
├── validation/              # dddh/ (jeu de référence, rapport), pocp/ (comparaison des moteurs)
├── branding/tga/            # kit du logo
├── stage_environment.yml
├── tga-pipeline-flowchart.png
└── README.md
```

Les dossiers `genomes/`, `proteins/`, `results/`, `logs/`, `webapp/jobs/` et les génomes téléchargés pour la
validation sont générés localement et exclus de git.

---

## 13. Limites connues et travaux restants

### Scientifique
- **Validation GGDC** : à faire (§ 10.2).
- **Zone 95–97 % d'ANI** : non représentée dans le jeu de validation. Les distances y passent de 0,031 à 0,058 de part
  et d'autre du seuil de 0,0412. Il faut ajouter des paires dans cette zone pour vérifier le calibrage du seuil.
- **Pourcentage dDDH linéaire** : il est maintenant borné entre 0 et 100 %, mais il reste une estimation peu fiable
  loin du seuil. Pour publier, il vaut mieux donner la distance `d`.
- **POCP** : la comparaison avec BLASTP porte sur 2 paires, aucune proche du seuil de 50 %.
- **Arbre** : FastTree fournit des supports SH-like locaux. Pour une publication, IQ-TREE avec bootstrap serait plus
  solide.
- **QC** : pas de N50, et pas de seuil bloquant pour les génomes fragmentés ou contaminés.
- **Synthèse** : il manque un tableau qui réunisse ANI, dDDH, AAI, POCP et le verdict, paire par paire.
- **Étude *B. ndiopicus*** : à relancer avec la version actuelle (§ 9.2).

### Logiciel
- **Versions** : l'environnement conda ne fixe pas les versions et liste des paquets non utilisés (streamlit, seaborn,
  matplotlib, pandas, numpy, openpyxl, pheatmap, treeio).
- **Installation** : elle se fait en trois parties (deux environnements conda et EasyCGTree4 à la main). Un
  `Dockerfile` ou un script d'installation la simplifierait.
- **Versions des outils** : elles ne sont pas enregistrées dans les dossiers de résultats.
- **Threads** : fixés à 4 pour CheckM2 et BLAST.
- **Scripts obsolètes** : `etape_ani.py`, `plot_ani.R`, `plot_aai.R`, `plot_pocp.R` et `plot_tree.R` contiennent des
  chemins en dur et sont à supprimer.
- **Langue** : le français et l'anglais sont mélangés dans les noms de fonctions et de fichiers de sortie.
- **Passage à l'échelle** : dDDH et POCP traitent les paires une par une. 50 génomes font 1 225 paires.
- **Interface** : pas de limite de taille à l'envoi, et les jobs ne sont jamais nettoyés. C'est acceptable en local,
  mais pas pour un déploiement public.
- **Tests** : seul l'enracinement est testé, et il n'y a pas d'intégration continue.
- **Nom** : il n'est pas encore unifié (« Genome-Taxonomy-Pipeline », « Genome Taxonomy Explorer » dans l'API, « Taxo
  Genomics Analyses » dans l'interface).

### Publication
- **Licence** : `LICENSE` manque.
- **Citation** : il faut un `CITATION.cff`, une version `v1.0.0` et un DOI Zenodo.
- **Données d'exemple** : un jeu d'exemple avec les résultats attendus.
- **Positionnement** : un « statement of need » face à TYGS, GTDB-Tk, pyANI et POCP-nf.
- **Revue cible** : JOSS (logiciel), ou une note d'application dans *Bioinformatics* ou *BMC Bioinformatics*.

### Feuille de route du README
- module d'identification par l'ARNr 16S (barrnap et base de référence) ;
- validation croisée avec d'autres outils (pyANI-plus, IQ-TREE) ;
- déploiement public (Docker et systemd).

---

## 14. Références

- Auch AF, Klenk HP, Göker M (2010). Standard operating procedure for calculating genome-to-genome distances based on
  high-scoring segment pairs. *Standards in Genomic Sciences* 2:142–148.
- Meier-Kolthoff JP, Auch AF, Klenk HP, Göker M (2013). Genome sequence-based species delimitation with confidence
  intervals and improved distance functions. *BMC Bioinformatics* 14:60.
- Jain C et al. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries.
  *Nature Communications* (FastANI).
- Qin QL et al. (2014). A proposed genus boundary for the prokaryotes based on genomic insights. *Journal of
  Bacteriology* (POCP).
- Kim D et al. (2021). Introducing EzAAI: a pipeline for high throughput calculations of prokaryotic average amino acid
  identity. *Journal of Microbiology*.
- Chklovski A et al. (2023). CheckM2. *Nature Methods*.
- Parks DH et al. (2018). A standardized bacterial taxonomy based on genome phylogeny. *Nature Biotechnology* (GTDB,
  bac120).
- Price MN, Dehal PS, Arkin AP (2010). FastTree 2. *PLoS ONE*.
- Hyatt D et al. (2010). Prodigal. *BMC Bioinformatics*.
- Buchfink B, Reuter K, Drost HG (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature
  Methods*.
- EasyCGTree4 : https://github.com/zdf1987/EasyCGTree4
