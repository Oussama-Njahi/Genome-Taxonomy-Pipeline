# Genome Taxonomy Pipeline

An integrated bioinformatics pipeline for the genome-based taxonomic characterization of bacteria, combining **ANI**, **AAI**, **POCP**, and core-gene **phylogenomics** in a single automated run.

Developed during a research internship at the **Institut des Régions Arides (IRA)**, Médenine, Tunisia (Laboratoire de Recherche "Écologie Pastorale").

## What it does

Given a folder of bacterial genome assemblies (FASTA), the pipeline runs seven steps end-to-end:

| Step | Task | Tool |
|------|------|------|
| 1 | Quality control (size, GC%, contigs, N50, completeness, contamination) | custom script + [CheckM2](https://github.com/chklovski/CheckM2) |
| 2 | Protein prediction | [Prodigal](https://github.com/hyattpd/Prodigal) |
| 3 | Average Nucleotide Identity (species delineation, 95% threshold) | [FastANI](https://github.com/ParBLiSS/FastANI) |
| 4 | Average Amino acid Identity (genus-level relationships) | [EzAAI](https://github.com/endixk/ezaai) |
| 5 | Percentage of Conserved Proteins (genus delineation, 50% threshold) | custom script + [DIAMOND](https://github.com/bbuchfink/diamond) |
| 6 | Phylogenomic tree (120 bacterial core genes, GTDB bac120 set) | [EasyCGTree4](https://github.com/zdf1987/EasyCGTree4) (Perl) + FastTree |
| 7 | Figures (heatmaps + annotated tree) | R (ggplot2, ggtree) |

A single command runs the full chain:

```bash
python scripts/pipeline.py
```

## Architecture

**Python** orchestrates the whole pipeline via `subprocess`, calling out to **Perl** (EasyCGTree, for the phylogenomic tree) and **R** (ggplot2/ggtree, for the figures) — the same multi-language integration pattern used by reference pipelines such as [PGCGAP](https://github.com/liaochenlanruo/pgcgap), extended here with AAI and POCP, which PGCGAP does not provide.

## Key finding (validation dataset)

On a 14-genome test set (*Bhargavaea*, *Candidatus Kurthia*, *Caryophanon*, plus four outgroups), the four independent analyses consistently show that ***Bacillus ndiopicus*** is genomically closer to *Caryophanon* (AAI ≈ 69%, POCP ≈ 50–54%) than to *Bacillus subtilis* (AAI ≈ 60%, POCP ≈ 38%) — suggesting a taxonomic misclassification, confirmed independently by ANI, AAI, POCP, and the phylogenomic tree topology.

## Requirements

- Linux (tested on Ubuntu via WSL2)
- [Miniforge](https://github.com/conda-forge/miniforge) (conda/mamba)
- ~6 GB RAM recommended for the quality-control step (CheckM2)

## Installation

```bash
mamba env create -f stage_environment.yml
conda activate stage
```

CheckM2 requires a separate environment due to a TensorFlow/Python version conflict:

```bash
mamba create -n checkm2env -y -c conda-forge -c bioconda python=3.12 checkm2
```

EasyCGTree4 is not distributed via conda and must be installed manually — see [its repository](https://github.com/zdf1987/EasyCGTree4) for instructions (clone, unzip the Linux release and HMM profile parts, place the `.hmm` files in `EasyCGTree4/EasyCGTree.v4.2-Linux/HMM/`, `chmod +x` the `bin/` folder, and install the system package `libgomp1`).

## Usage

Place your genome assemblies (`.fas`) in `genomes/`, then run:

```bash
python scripts/pipeline.py
```

Results (quality report, ANI/AAI/POCP matrices, tree, figures) are written to `results/`.

## Project structure

```
.
├── scripts/
│   ├── pipeline.py       # main orchestrator (7 steps)
│   ├── pocp.py           # POCP calculation (DIAMOND-based)
│   └── plot_all.R        # figure generation (ggplot2 + ggtree)
├── stage_environment.yml # reproducible conda environment
└── README.md
```

*(The `genomes/`, `proteins/`, `results/` and `logs/` folders are generated locally and excluded from version control — see `.gitignore`.)*

## Roadmap

- [ ] 16S-based taxonomic identification module (barrnap + reference database)
- [ ] Cross-validation with alternative tools (pyANI-plus, IQ-TREE)
- [ ] Local web interface for interactive analysis (in progress, separate module)
- [ ] Public deployment (Docker + systemd)

## References

- Jain et al. (2018). *FastANI*. Nature Communications.
- Kim et al. (2021). *EzAAI*. Journal of Microbiology.
- Qin et al. (2014). *POCP genus boundary*. Journal of Bacteriology.
- Chklovski et al. (2023). *CheckM2*. Nature Methods.
- Parks et al. (2018). *GTDB*. Nature Biotechnology.
- Price et al. (2010). *FastTree 2*. PLoS ONE.

## Author

Oussama — Bioinformatics internship, June–August 2026, IRA Médenine (supervisor: Imed Sbissi).
