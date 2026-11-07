# Passage du pipeline TGA à Snakemake : journal des changements

- **Date :** 22 octobre 2026
- **Branche :** `snakemake` (commit `8a374e6`, créé à partir de `main` au commit `ea6a3f7`), fusionnée dans `main` le 23 octobre 2026 (Pull Request #1, commit de fusion `76fcd0f`), CI verte sur `main`
- **Bilan :** 16 fichiers, +600 / −284 lignes

## 1. Résumé

Les neuf étapes du pipeline sont maintenant un workflow Snakemake (`workflow/Snakefile`). Les options de la ligne de commande (`--steps`, `--threads`, `--outgroup`…) et l'utilisation de l'interface web ne changent pas. Sur le jeu de test de 5 génomes, les résultats sont identiques à ceux de l'ancien pipeline (section 6).

Ce que Snakemake apporte :

- **Dépendances gérées automatiquement.** Si une étape a besoin des fichiers d'une autre et qu'ils manquent, cette autre étape tourne aussi.
- **Reprise.** Les résultats d'un run précédent sur les mêmes génomes sont réutilisés.
- **Parallélisme.** Prodigal et `EzAAI extract` tournent génome par génome, en parallèle.
- **Environnement CheckM2.** Il est déclaré dans le Snakefile (`conda: "checkm2env"`) au lieu d'être appelé par `conda run`.

## 2. Installation

Snakemake **9.27.0** (paquet `snakemake-minimal`, bioconda) a été installé dans l'environnement `stage` :

```bash
mamba install -n stage -c conda-forge -c bioconda snakemake-minimal=9.27.0
```

- **Effet sur l'environnement :** 54 paquets ajoutés et un seul rétrogradé (`packaging` 26.3 → 25.0). Les outils déjà fixés (FastANI, BLAST+, DIAMOND, EzAAI, Prodigal, R…) ne bougent pas.
- **Pourquoi la 9.27.0 :** la 9.27.1 était sortie la veille et n'a pas été retenue. Les versions bioconda antérieures à 9.10 ne s'installent pas avec Python 3.14.
- La version est fixée dans `stage_environment.yml`.

## 3. Fichiers nouveaux

| Fichier | Rôle |
|---|---|
| `workflow/Snakefile` | Les 9 étapes en règles Snakemake : `qc`, `checkm2`, `prodigal` (une tâche par génome), `proteines`, `ani`, `dddh`, `ezaai_extract` (une tâche par génome), `aai`, `pocp`, `arbre`, `figures`. |
| `scripts/qc.py` | Écrit `qc.txt` (taille, %GC, contigs, N50). Reprend la complétude et la contamination du rapport de CheckM2 quand il existe, en lisant les colonnes par leur nom. |
| `scripts/arbre.py` | Arbre phylogénomique : EasyCGTree, puis renommage des feuilles et enracinement. C'est l'ancienne fonction `etape_arbre` de `pipeline_web.py`, déplacée sans changement de logique. |
| `tests/test_workflow.py` | Vérifie à blanc (`snakemake -n`) les tâches lancées pour plusieurs choix d'étapes. |

## 4. Fichiers modifiés

### `scripts/pipeline_web.py`

**Ce qu'il garde :**
- la vérification des génomes ;
- la préparation des dossiers et l'archivage dans `previous_runs/` ;
- l'écriture de `versions.txt`, qui note maintenant aussi la version de Snakemake ;
- la synthèse `synthese_paires.tsv`.

**Ce qu'il fait en plus :**
- `etapes_prevues` calcule les étapes qui vont tourner : celles demandées, plus celles dont elles ont besoin.
- Il écrit `results/snakemake_config.json`, puis lance Snakemake en recopiant sa sortie. Les lignes inutiles de Snakemake sont filtrées.
- Si une étape échoue, il la nomme dans le message d'erreur.
- **Annulation :** il transmet le SIGTERM à Snakemake et continue de lire sa sortie jusqu'à ce que Snakemake ait fini de s'arrêter.
- **Liens vers les génomes :** les liens de `results/input_genomes/` prennent la date du génome qu'ils pointent (voir 7.2).

**Ce qui a été retiré :** les fonctions `announce`, `lire_checkm2`, `etape_qc`, `etape_checkm2`, `etape_proteines`, `etape_ani`, `etape_dddh`, `etape_aai`, `etape_pocp`, `etape_arbre`, `cleanup_easycgtree`, `nom_easycgtree` et `etape_figures`. Leur travail est maintenant fait par le Snakefile, `qc.py` et `arbre.py`.

### `webapp/backend/job_manager.py`

- Il lit le nouveau format de progression (section 5).
- L'état d'une analyse contient `running_steps` et `done_steps` à la place de `step` et `total`.
- La liste des étapes enregistrée dans `status.json` inclut les étapes ajoutées automatiquement.

### Interface web

- **`js/app.js` :** la barre et la liste des étapes affichent plusieurs étapes en cours à la fois. Le libellé devient par exemple « 1 of 7 set — Completeness & contamination, Average amino-acid identity ».
- **`css/style.css` :** les hachures rouges couvrent autant de cases que d'étapes en cours (`--r`).
- **`index.html` :** les notes de bas de page sur les dépendances entre étapes sont réécrites.

### Autres fichiers

- **`scripts/pipeline.py` :** l'aide décrit les nouvelles dépendances.
- **`tests/test_pipeline_inputs.py` :** `nom_easycgtree` est importé depuis `arbre.py`. Un test de `qc.py` est ajouté : colonnes de CheckM2 dans un autre ordre, génome absent du rapport.
- **`.github/workflows/tests.yml` :** installe Snakemake et lance `tests/test_workflow.py`.
- **`.gitignore` :** ignore `.snakemake/`.
- **`README.md` :**
  - sections Architecture, Outputs, Project structure et Tests mises à jour ;
  - comportement de `--steps` précisé ;
  - ajout de la commande de téléchargement de la base CheckM2 (~3 Go), qui manquait.
- **`RAPPORT_PROJET_TGA.md` :**
  - chronologie (ligne du 22 oct.), schéma d'architecture, progression, tableau des versions et dépendances entre étapes mis à jour ;
  - corrigé : « 4 threads » pour CheckM2 devient `--threads` ;
  - corrigé : la version de CheckM2 est fixée et la commande de téléchargement de sa base est ajoutée.

## 5. Ce qui change à l'usage

| | Avant | Maintenant |
|---|---|---|
| `--steps pocp` sans protéines | POCP sauté, avec un avertissement | Prodigal tourne d'abord, puis POCP |
| `--steps pocp` avec des protéines existantes | POCP seul | POCP seul : les protéines sont réutilisées |
| `--steps checkm2` sans `qc.txt` | CheckM2 sauté, avec un avertissement | CheckM2 tourne et `qc.txt` est refait avec ses valeurs |
| `--steps qc` après un run de CheckM2 | Complétude et contamination remises à NA | Les valeurs de CheckM2 sont conservées |
| Ordre d'exécution | Fixe, étapes 1 à 9 | Selon les dépendances, Prodigal et EzAAI en parallèle |
| Étape demandée | Relancée | Relancée (`--forcerun`) |
| Lignes de progression | `PROGRESS n/total <libellé>` | `STEPS a,b,c`, puis `PROGRESS start <étape> <libellé>` et `PROGRESS done <étape> <libellé>` |
| Fichiers de sortie en plus | — | `results/snakemake_config.json` et `.snakemake/log/` (journal de Snakemake) |
| Fichiers de sortie en moins | `results/checkm2_input_fna/` (copie des génomes) | CheckM2 lit directement `results/input_genomes/` |

Une étape choisie est toujours relancée. Une étape seulement nécessaire à une autre ne tourne que si ses fichiers manquent ou sont plus anciens que les génomes. Snakemake ne décide que sur la date des fichiers (`--rerun-triggers mtime`) : modifier le Snakefile ne refait donc pas les résultats déjà gardés.

## 6. Vérifications

**Analyse complète des 9 étapes sur les 5 génomes de `genomes/` : 9 min 40 s, aucune erreur.** Comparaison avec `results/`, produit par l'ancien pipeline au commit `ea6a3f7` :

| Fichier | Résultat |
|---|---|
| `qc.txt`, `dddh_distance.tsv`, `dddh_similarite.tsv`, `pocp_matrice.tsv`, `arbre.tree`, `arbre_non_enracine.tree`, `synthese_paires.tsv` | identiques octet par octet |
| `ani_resultats.txt` | valeurs identiques ; seul le chemin du dossier diffère |
| `aai_resultats.tsv` | valeurs identiques ; le chemin et l'ordre des lignes diffèrent, parce qu'EzAAI liste ses bases dans leur ordre de création, qui varie quand elles sont créées en parallèle |
| `heatmap_ani.png`, `heatmap_dddh.png`, `heatmap_pocp.png`, `tree.png` | identiques octet par octet |
| `heatmap_aai.png` | identique à l'œil ; les octets diffèrent à cause de l'ordre de dessin des cases |

**Reprise :** relancer `--steps pocp` sur ce même dossier a pris 1 min 03 s. Prodigal n'a pas été relancé.

**Interface web**, sur une analyse `qc,checkm2,ani,aai,pocp,figures` :
- Prodigal est ajouté automatiquement à la liste des étapes ;
- deux étapes s'affichent en cours en même temps ;
- la page de résultats est complète (registre avec complétude et nombre de protéines, planches ANI, AAI et POCP).

**Annulation pendant le dDDH :** l'analyse passe à « cancelled » et il ne reste aucun processus (Snakemake, `dddh.py`, `blastn`).

**Tests unitaires :** `test_root_tree.py`, `test_pipeline_inputs.py` (dont le nouveau `test_qc_et_checkm2`) et `test_workflow.py` (4 tests) passent tous.

## 7. Problèmes rencontrés et corrigés

1. **`sys.exit()` dans un bloc `run:` bloque Snakemake 9.** La tâche ne se termine jamais. Toutes les règles appellent donc des scripts en `shell:`, sur le modèle de `dddh.py` et `pocp.py`.
2. **Snakemake compare aussi la date des liens symboliques.** Les liens de `results/input_genomes/` sont recréés à chaque analyse, donc tous les résultats paraissaient périmés et auraient été recalculés. Chaque lien prend maintenant la date du génome qu'il pointe. C'est `test_workflow.py` qui a révélé le problème.
3. **Snakemake restait bloqué après une annulation depuis l'interface.** `pipeline_web.py` mourait au SIGTERM, et plus personne ne lisait la sortie de Snakemake. Corrigé comme décrit en section 4.
4. **EzAAI utilise `/tmp/ezaai` par défaut.** Les extractions lancées en parallèle se seraient partagé ce dossier. Chacune a maintenant son propre dossier temporaire.

## 8. Suites

Les trois points qui restaient ouverts le 22 octobre sont faits :

- **CI GitHub :** la branche a été poussée et la CI est passée sur la Pull Request #1, y compris l'installation de Snakemake par `pip` et la vérification de la syntaxe JavaScript (`node --check`), qui n'avait pas pu tourner en local.
- **Fusion :** la branche `snakemake` a été fusionnée dans `main` le 23 octobre 2026 par la Pull Request #1 (commit `76fcd0f`). La CI est verte sur `main`.
- **Correction après la fusion (30 octobre) :** `arbre.py` supposait un arbre construit sur les 120 gènes bac120. EasyCGTree en garde moins quand des génomes n'ont pas tous les gènes, et nomme alors l'arbre autrement (`.bac120.118.` avec deux *Klebsiella*). Le nombre de gènes est maintenant lu dans le nom du fichier ; un test couvre ce cas.

**Pistes possibles**, non faites :
- signaler dans l'interface qu'une étape a été ajoutée automatiquement ; aujourd'hui, seul le journal le dit ;
- utiliser `snakemake --report` pour produire un rapport HTML du run.
