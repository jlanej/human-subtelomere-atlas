# Inputs of the catalogue pipeline

Every file the catalogue pipeline (`pipeline/catalogue_*.py`, `catalogue_02_fetch_windows.sh`) and the fetch/verify scripts read from the
original analysis home (`ATLAS_INPUT_DIR`), by path relative to that directory. Upstream scripts are those of the analysis home
(`scripts/`, named by file name; not part of this repository). Sizes and SHA-256 are of the analysis-home files
(`pipeline/snapshot_inputs.py`). Everything else the pipeline uses is either public and fetched by range reads (Release 2
and HG002 Q100 windows, arm probes: steps 01-02) or produced by the pipeline itself.

Public mode (`ATLAS_INPUT_DIR` unset): `catalogue_common.inp(path)` resolves a path to `data/inputs/<path>[.gz]` (committed
snapshot; the `.gz` files decompress to bytes identical to the analysis-home files) or to `work/public/<path>` (downloaded by
`pipeline/fetch_public_inputs.sh`, checked against `data/inputs/public_files.tsv`).

| path | used by step | content | public source | upstream script(s) | bytes | sha256 | public mode |
|---|---|---|---|---|---|---|---|
| `tipdup/tips_meta.tsv` | 01, 03 | one row per Release 2 telomeric tip (20,144): assembly, contig, side, length, contig-name arm, tip-table junk/telomere lengths, assembler, HiFi input | HPRC Release 2 assemblies (459 assembly files scanned: 451 hifiasm, 8 Verkko; 455 carry telomeric tips) | stream_one.sh + stream_ends.py (outer 20 kb of every contig end -> r2scan/), telomere calls -> results/release2_all_telomeric_ends.tsv, tip_dup_queries.py | 3,028,871 | `68553412a3d4e610108064ba32a8cf7be032bafdcbccf91ed8081aa2fa3d7bf1` | snapshot |
| `tipdup/tip_placement.tsv` | 09 | placement of each 20 kb tip on CHM13 (own/foreign end, junk hits); used for the junk/chimeric flags | HPRC Release 2 assemblies; CHM13 v2.0 | tip_dup_classify.py | 4,436,613 | `b0284ca2b22f77694df07e3dad5ff319492618ab67cf607625c58f1e1c489684` | snapshot |
| `meta/r2_assemblies_index.csv` | 01 | S3 paths (submission and release2 working copy) of every Release 2 assembly FASTA | listing of the public HPRC bucket (s3://human-pangenomics/working/HPRC/*/assemblies/release2/) | bucket listing saved in the analysis home (no script) | 473,200 | `6779354f3b08c4fa274d44bdc6d20bb4b5dc2dd815730e79de4ce4efe9ada630` | snapshot |
| `meta/explorer_assemblies.csv` | 03 | assembly metadata (GenBank accession, assembly name) by assembly file | HPRC data explorer assembly catalogue (public) | downloaded catalogue (no script) | 276,878 | `bc6653f68fb069f9976b4453ba5b6e670b872acdbbfb42a77fbb12ace1de33c5` | snapshot |
| `results/subtel_foreign_tips.tsv` | 09 | earlier 20 kb calls of the named tips: t20 state/best end, prior form and junction clusters, Verkko concordance (vk_*), hifiasm read-layout support (lay_*) | HPRC Release 2 hifiasm and Verkko assemblies and assembly graphs | subtel_foreign_forms.py, which merges the outputs of subtel_foreign_verkko.py (vk_*) and subtel_foreign_layout.py (lay_*); driver subtel_foreign_reproduce.sh | 3,924,423 | `86ce6bb7ad8f6436a798dc5b952e418771776728d3285ff9be6df80225d5c294` | snapshot |
| `results/subtel_hg002_ends.tsv` | 09 | HG002 Q100 ends: read support of switches (HiFi/ONT spanning and crossing reads) | HG002 Q100 v1.1 and public HG002 read alignments | subtel_foreign_hg002_reads.py, subtel_foreign_hg002_junction.py | 22,924 | `e276e282b6e652803008c66c655243bed40ab9cd8f97ab697c05ef26f8f0181b` | snapshot |
| `results/subtel_pedigree_transmission.tsv` | 09 | WashU pedigree transmissions: child end, transmitting parental end, intactness | WashU pedigree v1.1 assemblies | subtel_ped_mendel.py | 51,023 | `9ba4bd65ce5e5c62643003b734e87deca16e01fa49ef3d3eb6119e4cfa3f86ff` | snapshot |
| `hg002/hg002v1.1.fasta.gz.fai` | 02, 03 | FASTA index of HG002 Q100 v1.1 (chromosome lengths for the 92 end windows) | https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/hg002v1.1.fasta.gz.fai | downloaded (public index) | 1,912 | `836220cfb3ca888205995c68b15c3dbcd2f467317b034c7eeaa3335d9576f73e` | snapshot |
| `subtel/pedigree/PAN027_pat.v1.1.fai` | fetch | FASTA index of a pedigree haplotype (byte offsets for range reads); 8 files, one per haplotype | https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/ | subtel_ped_remote_fai.py (PAN028_h2: subtel_ped_remote_fai_search.py) | 1,097 | `d317208c59bcd694a030e5bf7a7d323fc2ac40a313370be3130907767e019aae` | snapshot (all 8) |
| `ref/chm13v2.0.fa.gz` | 03, 04, 05 | T2T-CHM13 v2.0 (GCA_009914755.4), bgzipped analysis set: 48 end windows, 500 kb ends, probe-mapping index | https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz | downloaded | 981,185,361 | `f274b0ee8bc1bca18d231d0eb80d764650d7f02d1fe07c43819e21168224304d` | download (fetch_public_inputs.sh chm13); .fai/.gzi rebuilt by samtools |
| `subtel/pedigree/ends250k.fa` | 03 | outer 250 kb of both ends of every chromosome, 7 pedigree haplotypes (322 ends) | https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/ | subtel_ped_fetch_ends.py | 80,509,502 | `fc9bccc13d09bb8b8aa7b007ad7d78fcf9e94df6bceb1a4de31a450fa8ac0388` | range reads (fetch_public_inputs.sh pedigree; byte-identical) |
| `subtel/pedigree/ends250k_PAN028_h2.fa` | 03 | the same for PAN028 haplotype 2 (46 ends) | https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/ | subtel_ped_fetch_ends.py | 11,501,344 | `4d09ceb3a7a427804d0f2a64151a2e1a6c92e0384ecafe69ae57ea28caeee69b` | range reads (as above) |
| `tipdup/tips_masked.fa` | 03 (optional) | 20 kb tip windows with the telomere masked; only for a sanity check of the first 20 kb of every Release 2 window | HPRC Release 2 assemblies | tip_dup_queries.py | 403,851,872 | `8c171fce7c9430a93e67b06a4847da8e4e47f63138d2a2821c6dd744d520c104` | not needed (check skipped) |
| `subtel/foreign/r2_100k/*.fa (455 files)` | 01, 03 (optional) | earlier range reads of the outer 100 kb of tips; reused as a cache | HPRC Release 2 assemblies | subtel_foreign_fetch100.py, subtel_foreign_fetch100.sh | 593,831,842 | `aggregate:214b8ffae5dd00603802c7611984ad904ecc31e5e8d118825e29b8d9291de875` | not needed (re-fetched by steps 01-02) |

Pedigree index files (all in `data/inputs/subtel/pedigree/`):

| path | bytes | sha256 |
|---|---|---|
| `subtel/pedigree/PAN027_pat.v1.1.fai` | 1,097 | `d317208c59bcd694a030e5bf7a7d323fc2ac40a313370be3130907767e019aae` |
| `subtel/pedigree/PAN027_mat.v1.1.fai` | 1,097 | `b54b93cea4d3485453d20afaec40013f070e78c5027d03dcb13337f84dc3ed1f` |
| `subtel/pedigree/PAN010_h1.v1.1.fai` | 1,144 | `878c1a8d913f96a46ce90d057531edfff7664e43e381f4e63b0342559e41e76e` |
| `subtel/pedigree/PAN010_h2.v1.1.fai` | 1,144 | `503fa7d278da4c26e1c00c4911715a12acb24d0e49839e83c55552ead22e288f` |
| `subtel/pedigree/PAN011_h1.v1.1.fai` | 1,143 | `9b5972cb8fb1125dde0fb69dd17236bf9123341cfdec823e535a7ebcaf551c92` |
| `subtel/pedigree/PAN011_h2.v1.1.fai` | 1,142 | `cb50edbf08d285d7d46d8bf4c0dbae40fea2f9e6a9af6f9431bca4abbaf0912f` |
| `subtel/pedigree/PAN028_h1.v1.1.fai` | 1,143 | `919f3dc8f4d35e0a3e8756929c2859609ec662da2cb2fea019892331776a4e18` |
| `subtel/pedigree/PAN028_h2.v1.1.fai` | 1,143 | `038b2a12b19800bf04ed92056f6a6667f47085e09aec533487e2edf93480915f` |

CHM13 index files of the analysis home (`ref/chm13v2.0.fa.gz.fai`, 788 bytes, `cec84e8875f740fea679b08cb042529deae56f8276f567c2d5dca232d6172c1b`; `.gzi`, 773,592 bytes, `cb9ee1bd52f325e9ff93bb6c8d584a5954d745cb25a69882995a40e1ef528387`) are regenerated by `samtools faidx` after the download.
The 455 cached window files are summarised by an aggregate checksum (SHA-256 over the sorted list of file name and file SHA-256):
`aggregate:214b8ffae5dd00603802c7611984ad904ecc31e5e8d118825e29b8d9291de875`.

## What cannot be rebuilt without the analysis home

* **Tip discovery.** `tipdup/tips_meta.tsv` (which contig ends of the 459 Release 2 assembly files carry a telomere) comes from a
  streaming scan of every assembly (`stream_one.sh` / `stream_ends.py`: outer 20 kb of every contig end, 3.2 GB of windows),
  telomere calls on those windows and `tip_dup_queries.py`. This scan is not re-implemented here; its result is committed as the
  snapshot. To reproduce it: stream each assembly listed in `meta/r2_assemblies_index.csv` (about 411 GB of compressed FASTA by the HPRC catalogue file sizes), keep
  the outer 20 kb of each contig end, call the telomere with `subtel_atlas.telomere.find_block` and keep ends with a telomere.
  The catalogue then re-calls every telomere boundary on its own 100 kb window, so only the list of tips is inherited.
* **Prior evidence tables** (`results/subtel_foreign_tips.tsv`: Verkko concordance and read-layout support; `subtel_hg002_ends.tsv`:
  HG002 read support; `subtel_pedigree_transmission.tsv`; `tipdup/tip_placement.tsv`) come from earlier analyses of public
  assemblies, assembly graphs and read alignments (scripts listed above) and are used as evidence columns and flags only. They are
  committed as snapshots; recomputing them requires those analyses (read alignments of tens of GB per sample).
* **Optional caches** (`tipdup/tips_masked.fa`, `subtel/foreign/r2_100k/`) are not needed: windows are re-fetched, and the 20 kb
  sanity check of step 03 is skipped.
* The **arm-inference probes, all windows and the CHM13 / HG002 / pedigree sequences** are re-fetched from public sources.
