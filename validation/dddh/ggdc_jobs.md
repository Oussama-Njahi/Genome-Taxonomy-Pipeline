# GGDC jobs

Eight jobs on https://ggdc.dsmz.de/ggdc.php cover 19 of the 37 pairs, including the 10 pairs nearest the species
threshold. Each job takes one query genome and its reference genomes; the files are in `genomes/<group>/`.

| Job | Group folder | Query genome | Reference genomes |
|---|---|---|---|
| 1 | `escherichia` | `Escherichia-coli-K12-MG1655.fna` | `Escherichia-coli-O157H7-Sakai.fna`, `Shigella-flexneri-2a-301.fna`, `Escherichia-fergusonii-ATCC35469.fna`, `Escherichia-albertii-BIA5-2.fna` |
| 2 | `escherichia` | `Escherichia-coli-O157H7-Sakai.fna` | `Shigella-flexneri-2a-301.fna` |
| 3 | `bacillus_cereus` | `Bacillus-anthracis-Ames.fna` | `Bacillus-thuringiensis-97-27.fna`, `Bacillus-cereus-ATCC14579.fna`, `Bacillus-mycoides-KBAB4.fna` |
| 4 | `klebsiella` | `Klebsiella-pneumoniae-HS11286.fna` | `Klebsiella-variicola-At-22.fna`, `Klebsiella-quasipneumoniae-FDAARGOS1503.fna`, `Klebsiella-oxytoca-NCTC13727.fna` |
| 5 | `klebsiella` | `Klebsiella-quasipneumoniae-FDAARGOS1503.fna` | `Klebsiella-variicola-At-22.fna` |
| 6 | `bacillus_subtilis` | `Bacillus-velezensis-FZB42.fna` | `Bacillus-amyloliquefaciens-DSM7.fna`, `Bacillus-subtilis-168.fna`, `Bacillus-licheniformis-ATCC14580.fna` |
| 7 | `bacillus_subtilis` | `Bacillus-spizizenii-TU-B-10.fna` | `Bacillus-subtilis-168.fna`, `Bacillus-subtilis-NCIB3610.fna` |
| 8 | `bacillus_subtilis` | `Bacillus-subtilis-168.fna` | `Bacillus-subtilis-NCIB3610.fna`, `Bacillus-licheniformis-ATCC14580.fna` |

For every result, write the **formula 2** dDDH (in %) in `ggdc_reference.tsv` on the line of that pair, then run
`validate_dddh.py` again.
