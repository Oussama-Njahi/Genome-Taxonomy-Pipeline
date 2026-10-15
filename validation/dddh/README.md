# dDDH validation

`scripts/dddh.py` is a home-made implementation of the GBDP formula 2 of Auch et al. (2010). This folder measures
how far its species decisions can be trusted, on reference genomes chosen to frame the species threshold: strains of
the same species, closely related species and distant species (`benchmark.tsv`, 19 complete RefSeq genomes in four
groups, 37 pairs).

## Run it

In the pipeline's conda environment:

```bash
python validation/dddh/validate_dddh.py
```

The script downloads the genomes from NCBI once (`genomes/`, not tracked by git), runs `scripts/dddh.py` unchanged and
FastANI on each group (`work/`, not tracked), then writes:

- `pairs.tsv`: one line per pair, with the dDDH distance, the decision, ANI and its aligned fraction;
- `RESULTS.md`: the agreement tables and the pairs near the threshold;
- `ggdc_reference.tsv`: a table to fill with GGDC's values (created empty the first time).

It takes about 20 minutes the first time. Later runs reuse the results in `work/`; add `--force` to recompute them
after changing `dddh.py`.

## What is compared

1. **dDDH against ANI.** A pair is called "same species" by dDDH when its distance is at most 0.0412, and by ANI when
   FastANI reports at least 95 %. The two criteria are computed independently, so their agreement is a first check that
   needs no outside service.
2. **dDDH against GGDC.** This is the comparison a reviewer will ask for, because GGDC is the reference
   implementation. It needs values from the GGDC web service, which have to be requested by hand (below).

The NCBI species names are in `pairs.tsv` for information only. They are not used as a reference: the benchmark
deliberately contains cases where the name and the genome disagree (*Shigella flexneri* is genomically *Escherichia
coli*; the strain named *Bacillus thuringiensis* 97-27 sits next to *B. anthracis*).

## Getting the GGDC values

1. Open https://ggdc.dsmz.de/ggdc.php. A job takes one query genome and up to 75 reference genomes.
2. For each of the eight jobs listed in `ggdc_jobs.md`, upload the query FASTA and the reference FASTA files from
   `genomes/<group>/` (uploading is faster than giving accession numbers), keep the recommended alignment tool
   (NCBI-BLAST BLAST+), enter your e-mail address and submit.
3. From each result, copy the **formula 2** dDDH estimate (in %) into the `ggdc_ddh_formula2` column of
   `ggdc_reference.tsv`, on the line of that pair. Pairs left empty are ignored.
4. Run `validate_dddh.py` again: `RESULTS.md` then also reports the agreement with GGDC (dDDH ≥ 70 %).

GGDC converts distances to a dDDH percentage with a logistic model, while `dddh.py` uses the linear conversion of Auch
et al. (2010). The percentages are therefore not expected to match; what is compared is the same-species decision and
the ranking of the pairs.
