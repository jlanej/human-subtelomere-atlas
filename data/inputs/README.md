# Snapshot of the small derived inputs (catalogue v0.1)

Committed copies of the small tables the catalogue pipeline reads, so that the catalogue can be rebuilt from public data without
the original analysis home. Paths mirror the analysis-home layout (`pipeline/INPUTS.md` lists every input, its public source,
the upstream script that produced it, its size and SHA-256). Tables are gzip-compressed (deterministic, mtime 0) and
decompress to bytes identical to the analysis-home files; index files are uncompressed. No file contains absolute or personal
paths (checked by `pipeline/snapshot_inputs.py`). All content derives from public data (HPRC Release 2, T2T CHM13 v2.0,
HG002 Q100 v1.1, WashU pedigree v1.1 assemblies and public read sets).

| file | content |
|---|---|
| tipdup/tips_meta.tsv.gz | the 20,144 Release 2 telomeric tips (assembly, contig, side, length, contig-name arm, tip-table junk/telomere lengths, assembler) |
| tipdup/tip_placement.tsv.gz | placement of each 20 kb tip on CHM13 (junk / chimeric flags) |
| meta/r2_assemblies_index.csv.gz | public S3 paths of the Release 2 assembly FASTA files |
| meta/explorer_assemblies.csv.gz | HPRC assembly catalogue (GenBank accession, assembly name) |
| results/subtel_foreign_tips.tsv.gz | earlier 20 kb calls of named tips, Verkko concordance (vk_*), read-layout support (lay_*) |
| results/subtel_hg002_ends.tsv.gz | HG002 Q100 ends: read support of switches |
| results/subtel_pedigree_transmission.tsv.gz | WashU pedigree transmissions and intactness |
| hg002/hg002v1.1.fasta.gz.fai | HG002 Q100 v1.1 FASTA index (public) |
| subtel/pedigree/*.v1.1.fai | FASTA indexes of the 8 pedigree haplotypes (byte offsets for range reads) |
| public_files.tsv | size and SHA-256 of the large public files downloaded by `pipeline/fetch_public_inputs.sh` |
| SHA256SUMS | checksums of every file in this directory |

Regenerate (maintainer, with the analysis home): `ATLAS_INPUT_DIR=<analysis home> python pipeline/snapshot_inputs.py`.
Verify: `cd data/inputs && shasum -a 256 -c SHA256SUMS` (or `sha256sum -c SHA256SUMS`).
