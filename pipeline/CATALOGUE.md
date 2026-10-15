# Allele catalogue pipeline

`bash pipeline/run_catalogue.sh` rebuilds `data/catalogue/` (release v0.1) from public assemblies and the analysis inputs
named by `ATLAS_INPUT_DIR` (set in the git-ignored `work/config.local.sh`; see `pipeline/config.sh`). Large intermediates are
written to `work/` (git-ignored). Requirements: samtools/htslib with libcurl, minimap2, Python 3 with numpy, pandas, scipy,
matplotlib, mappy and pysam.

| step | script | what it does |
|---|---|---|
| 1a | catalogue_01_plan_windows.py | 100 kb telomere-anchored window per Release 2 telomeric tip (20,144 tips; named and unnamed contigs); reuses earlier 100 kb range reads; plans 10 kb arm probes at 250 kb and 1 Mb for unnamed-contig tips and 1,000 random named tips (calibration) |
| 1b | catalogue_02_fetch_windows.sh | range reads from the public HPRC S3 bucket (samtools faidx on the bgzipped FASTA; 4 parallel jobs, 3 retries, .ok markers); HG002 Q100 v1.1 ends |
| 1c | catalogue_03_assemble_windows.py | orients every window telomere -> centromere, adds CHM13 v2.0 (48 ends), HG002 Q100 (92) and WashU pedigree ends (368, validation only), records provenance and SHA-256, checks the first 20 kb against the earlier tip windows |
| 2a | catalogue_04_align_bins.py | 5 kb bins tiled from the telomere/subtelomere boundary plus S5/S10/S20 (and legacy t20) windows aligned to the 48 CHM13 500 kb ends (minimap2 asm20) |
| 2b | catalogue_05_arms.py | home arm: contig name, else probe mapping to CHM13 (calibrated), else mosaic rule; pooled groups acroP/PAR1/PAR2 |
| 2c | catalogue_06_mosaic.py | bin labels (home / other CHM13 end / none), smoothing, blocks, base-level switch refinement, S5/S10/S20 states |
| 2d | catalogue_08_alleles.py | core structures (blocks >= 10 kb), structural alleles by structure + switch positions (complete linkage, 10 kb), assignment of partial / low-confidence / pedigree ends |
| 2e | catalogue_07_haplogroups.py | within-arm all-vs-all k-mer ANI (FracMinHash) of the first 80 kb of subtelomere, UPGMA at ANI 0.995, medoids, alignment validation |
| 3a | catalogue_09_release.py | evidence (Verkko, read layout, HG002 reads, pedigree), flags, tiers, frequencies by superpopulation, release tables, representatives, JSON |
| 3b | catalogue_10_qc.py | QC.md, figures/qc, PROVENANCE.json, MANIFEST.tsv, SHA256SUMS |

Interpretation note: switches (junctions) mark where the sequence changes which CHM13 end it most resembles. Inside shared,
paralogous subtelomeric sequence this is a change in paralog similarity, not necessarily a structural breakpoint. 'Foreign'
is distance-dependent and is reported per 5 kb bin (`profile`) and as boundary-anchored S5/S10/S20 states.

References: CHM13 v2.0 - Nurk et al. 2022, doi:10.1126/science.abj6987. minimap2 - Li 2018, doi:10.1093/bioinformatics/bty191.
Mash distance - Ondov et al. 2016, doi:10.1186/s13059-016-0997-x.
