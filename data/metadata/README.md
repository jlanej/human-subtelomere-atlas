# Sample metadata

`samples.tsv` — one row per sample in the HPRC Release 2 sample table (234 rows, including the CHM13 and GRCh38 reference entries),
built by `analyses/allele_ages/scripts/build_samples_tsv.py`.

| column | meaning |
|---|---|
| sample | HPRC sample ID |
| population / superpopulation | 1000 Genomes population and superpopulation code (AFR, AMR, EAS, EUR, SAS) |
| sex | sex as given in the HPRC Release 2 table (1kGP value used only when HPRC is blank) |
| source | which sources supplied the population labels |
| population_descriptor, biosample_id, alternative_id, paternal_id, maternal_id | copied from the HPRC Release 2 table |
| in_1kgp_3202 | sample present in the 1kGP 3,202-sample high-coverage table |
| has_telomeric_tips | sample has telomeric contig tips in this atlas (230 samples) |
| related_in_set | a parent or child of this sample is also in the table (none observed) |

## Sources (downloaded 2026-10-15; sha256 in `raw/SHA256SUMS`, URLs in `raw/DOWNLOADS.tsv`)

1. HPRC Release 2 sample metadata — https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/sample/hprc_release2_sample_metadata.csv
2. 1000 Genomes 3,202-sample population/pedigree table (Byrska-Bishop et al. 2022, doi:10.1016/j.cell.2022.08.004) — https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
3. IGSR population → superpopulation table (1000 Genomes Project Consortium 2015, doi:10.1038/nature15393) — https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv

## Manual assignments and source discrepancies

- MKK (Maasai in Kinyawa, Kenya; HapMap 3) and ASL (African Americans in St. Louis; WashU pedigree, HG06807 = PAN027) are not in the 1kGP phase-3 population table; both assigned AFR from the HPRC population descriptor.
- HG002 (GIAB Ashkenazi trio son) assigned population `ASJ_GIAB`, superpopulation EUR; HG005 (GIAB Han Chinese trio son) assigned `CHN_GIAB`, EAS. Neither has a 1kGP population code.
- The HPRC table lists HG002's parents as NA24694/NA24695 and HG005's as HG003/HG004. GIAB documents HG002's parents as HG003/HG004 and HG005's as HG006 (NA24694)/HG007 (NA24695), so the two entries look swapped. The values are copied unchanged; neither parent set is in the table, so relatedness calls are not affected.
- HG02300: sex is `female` in the HPRC table and male (code 1) in the 1kGP table; the HPRC value is reported.
- 7 of the 230 tip samples (HG005, HG01123, HG02109, HG02486, HG02559, HG03471, NA21309) are not in the 1kGP 3,202 table; their populations come from the HPRC table.
- Relatedness: no sample in the table has a parent or child also in the table (checked from the HPRC paternal/maternal IDs and 1kGP FamilyID). Cryptic second-degree relatedness was not assessed.
