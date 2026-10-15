# Genome Taxonomy Pipeline

An integrated bioinformatics pipeline for the genome-based taxonomic characterization of bacteria, combining **ANI**, **dDDH**, **AAI**, **POCP**, and core-gene **phylogenomics** in a single automated run.

Developed during a research internship at the **Institut des Régions Arides (IRA)**, Médenine, Tunisia (Laboratoire de Recherche "Écologie Pastorale").

## What it does

Given a folder of bacterial genome assemblies (FASTA), the pipeline runs nine steps end-to-end — every step is optional, see `--steps` below:

| Step | Task | Tool |
|------|------|------|
| 1 | Quality control (size, GC%, contigs) | custom script |
| 2 | Completeness / contamination (needs step 1) | [CheckM2](https://github.com/chklovski/CheckM2) |
| 3 | Protein prediction | [Prodigal](https://github.com/hyattpd/Prodigal) |
| 4 | Average Nucleotide Identity (species delineation, 95% threshold) | [FastANI](https://github.com/ParBLiSS/FastANI) |
| 5 | digital DNA-DNA Hybridization (species delineation, GBDP formula 2, Auch et al. 2010 method) | custom script + BLASTN (`-task blastn`) |
| 6 | Average Amino acid Identity (genus-level relationships) | [EzAAI](https://github.com/endixk/ezaai) |
| 7 | Percentage of Conserved Proteins (genus delineation, 50% threshold, needs step 3) | custom script + [DIAMOND](https://github.com/bbuchfink/diamond) (`--ultra-sensitive`) |
| 8 | Phylogenomic tree (120 bacterial core genes, GTDB bac120 set), rooted at its midpoint or on an outgroup | [EasyCGTree4](https://github.com/zdf1987/EasyCGTree4) (Perl) + FastTree |
| 9 | Figures (heatmaps + annotated tree) | R (ggplot2, ggtree) |

A single command runs the full chain:

```bash
python scripts/pipeline.py
```

Run only a subset with `--steps` (comma-separated, e.g. `--steps qc,ani,dddh`) — useful to skip CheckM2 (RAM-heavy, see Requirements) while still getting the rest.

FastTree writes an unrooted tree, so the pipeline roots it before drawing it: at its midpoint by default, or on an outgroup with `--outgroup` (one or more genome names without the file extension, comma-separated, e.g. `--outgroup Staphylococcus-aureus-NCTC8325`). Prefer an outgroup when you have one. `results/arbre.tree` is the rooted tree; the tree as computed is kept in `results/arbre_non_enracine.tree`.

## Architecture

**Python** orchestrates the whole pipeline via `subprocess`, calling out to **Perl** (EasyCGTree, for the phylogenomic tree) and **R** (ggplot2/ggtree, for the figures) — the same multi-language integration pattern used by reference pipelines such as [PGCGAP](https://github.com/liaochenlanruo/pgcgap), extended here with AAI and POCP, which PGCGAP does not provide.

## Key finding (validation dataset)

On a 14-genome test set (*Bhargavaea*, *Candidatus Kurthia*, *Caryophanon*, plus four outgroups), the four independent analyses consistently show that ***Bacillus ndiopicus*** is genomically closer to *Caryophanon* (AAI ≈ 69%, POCP ≈ 50–54%) than to *Bacillus subtilis* (AAI ≈ 60%, POCP ≈ 38%) — suggesting a taxonomic misclassification, confirmed independently by ANI, AAI, POCP, and the phylogenomic tree topology.

## Notes on the methods

- **ANI and aligned fraction.** FastANI reports how many fragments of the query genome it could align. The interface shows this aligned fraction under each ANI value, and both the interface and the ANI figure mark a "same species" ANI (≥ 95%) measured on less than 50% of the genome, which is too little to conclude on its own.
- **dDDH uses classic BLASTN.** BLAST+ runs megablast when `blastn` is called without `-task`; megablast only aligns the most conserved regions and underestimates the distance between different species (on the test set, an outgroup from another genus came out *closer* than species of the same genus). `dddh.py` therefore sets `-task blastn`.
- **POCP uses DIAMOND in `--ultra-sensitive` mode.** POCP was defined with BLASTP (Qin et al. 2014); DIAMOND's default fast mode misses part of the divergent proteins and lowers POCP by about one point on the test set.

## Note on dDDH values

`scripts/dddh.py` implements the method of Auch et al. (2010), Table 4, "NCBI BLAST / Trimming (2)": GBDP formula 2 (identities within HSPs / total HSP length), greedy-with-trimming, linear conversion `DDH ≈ 90.3998 − 438.3134 × d` and species threshold `d > 0.0412` (≈ <70% DDH). Since version 2.0, the [GGDC](https://ggdc.dsmz.de) web service (Meier-Kolthoff et al. 2013) converts distances with a logistic model (GLM) and uses other settings (BLAST+ with word length 38, "coverage" algorithm); the coefficients of that model are not published in the article. The pipeline's values are therefore **not directly comparable** with those of today's GGDC, especially far from the threshold, which is why the figures show the distance rather than the estimated %DDH. To validate a conclusion, submit the genome pairs whose distance is close to 0.0412 to GGDC and compare the same-species / different-species decisions.

The distance threshold and the conversion line are two separate results of Auch et al.'s Table 4: the threshold (0.0412) is the distance that best separates pairs around 70% DDH (error ratio 4.3%), while the line comes from a robust regression over all pairs. They do not meet exactly: at d = 0.0412 the line gives 72.3%, and 70% would correspond to d = 0.0465. As the paper recommends, the species decision uses the distance threshold; the percentage is only an estimate. The line is calibrated near the threshold and would fall below 0% for distant species, so the estimate is bounded to 0–100%.

## Requirements

- Linux (tested on Ubuntu via WSL2)
- [Miniforge](https://github.com/conda-forge/miniforge) (conda/mamba)
- ~6 GB RAM recommended for the CheckM2 step (skip it with `--steps` if RAM is limited — see below)

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

Place your genome assemblies (`.fas`, `.fasta` or `.fna`, in any case; at least 2, or 5 for the phylogenomic tree) in `genomes/`, then run:

```bash
python scripts/pipeline.py
```

Results (quality report, ANI/dDDH/AAI/POCP matrices, tree, figures) are written to `results/`, predicted proteins to `proteins/`. `--steps` works as described above; `--threads` defaults to all CPUs.

Before any step, every genome is checked: it must be a non-empty DNA FASTA file, and its name (the file name without the extension) may only contain letters, digits, `.`, `_` and `-`, because it becomes the genome's label in every table and in the tree. Two files may not give the same name. The run stops with a message naming the file otherwise.

Results are kept from one run to the next, so steps can be run separately (`--steps qc`, then `--steps checkm2`), but only for the same genomes: `results/genomes.tsv` records which genomes were analysed, and when they change, the previous `results/` and `proteins/` are moved to `previous_runs/<date>/` so that old and new results are never mixed.

## Web interface

A local web UI wraps the same 9-step pipeline behind an upload form, so genomes can be dropped in a browser instead of the command line:

```bash
cd webapp/backend
uvicorn app:app --reload
```

Then open `http://127.0.0.1:8000`. Upload at least 2 genome assemblies (`.fas`, `.fasta`, `.fna`; 5 if the tree is among the steps), pick which steps to run via checkboxes, and the results (quality report, ANI/dDDH/AAI/POCP matrices, tree, heatmaps) become viewable and downloadable once the job finishes. A running job can be cancelled from the UI at any time.

## Pipeline overview

<p align="center">
  <img
    src="tga-pipeline-flowchart.png"
    alt="Genome Taxonomy Pipeline overview"
    width="100%"
  >
</p>
## Project structure

```
.
├── scripts/
│   ├── pipeline.py       # command-line entry point: runs pipeline_web.py on genomes/
│   ├── pipeline_web.py   # the 9-step pipeline (also used by the web interface)
│   ├── pocp.py           # POCP calculation (DIAMOND-based)
│   ├── dddh.py           # dDDH calculation (BLASTN-based, GBDP formula 2)
│   ├── root_tree.py      # roots the tree (midpoint or outgroup), keeping branch supports
│   └── plot_all_web.R    # figure generation (ggplot2 + ggtree)
├── webapp/
│   ├── backend/          # FastAPI server (app.py, job_manager.py)
│   └── frontend/         # upload form + results UI (HTML/CSS/JS, no build step)
├── stage_environment.yml # reproducible conda environment
└── README.md
```

*(The `genomes/`, `proteins/`, `results/` and `logs/` folders are generated locally and excluded from version control — see `.gitignore`.)*

## Roadmap

- [ ] 16S-based taxonomic identification module (barrnap + reference database)
- [ ] Cross-validation with alternative tools (pyANI-plus, IQ-TREE)
- [ ] Public deployment (Docker + systemd)

## References

- Jain et al. (2018). *FastANI*. Nature Communications.
- Auch, Klenk & Göker (2010). *Standard operating procedure for calculating genome-to-genome distances based on high-scoring segment pairs*. Standards in Genomic Sciences 2:142-148. (dDDH, GBDP formula 2)
- Meier-Kolthoff, Auch, Klenk & Göker (2013). *Genome sequence-based species delimitation with confidence intervals and improved distance functions*. BMC Bioinformatics 14:60. (current GGDC method)
- Kim et al. (2021). *EzAAI*. Journal of Microbiology.
- Qin et al. (2014). *POCP genus boundary*. Journal of Bacteriology.
- Chklovski et al. (2023). *CheckM2*. Nature Methods.
- Parks et al. (2018). *GTDB*. Nature Biotechnology.
- Price et al. (2010). *FastTree 2*. PLoS ONE.

## Author

Oussama Njahi — Bioinformatics internship, June 2026, IRA Médenine (supervisor: Imed Sbissi).
