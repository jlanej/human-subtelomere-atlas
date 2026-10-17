# Read-evidence showcase (in-repository subset)

Eight exemplar claims, with up to 6 definitive reads each (chosen one per dataset first, then by number of diagnostic units).
The full evidence set (every v1-conclusive read of all 122 claims with a conclusive read) is
packaged as GitHub release assets (table below). Directory size: 16.8 MB. Checksums: [SHA256SUMS](SHA256SUMS).

Columns: *distal on D* = distal read segments (telomere boundary to the switch) whose genome-wide primary alignment is on the
donor chromosome; *proximal on R* = proximal segments aligned to the recipient chromosome at MAPQ 60; *own end* = reads whose
best own-assembly window is the claimed end (`*.segments.tsv`, `*.reads.tsv`).

| claim | individual | arm | allele | model | v1 | definitive tier | reads here | distal on D | proximal on R (MAPQ 60) | own end |
|---|---|---|---|---|---|---|---:|---:|---:|---:|
| [HG002#2#chr3_MATERNAL#R](HG002_2_chr3_MATERNAL_R/README.md) | HG002 | 3q | SA-3q-01 | CHM13 11p to 3q | confirmed | definitive multi dataset | 3 of 3 | 3/3 | 3/3 | 3/3 |
| [HG002#2#chr9_MATERNAL#L](HG002_2_chr9_MATERNAL_L/README.md) | HG002 | 9p | SA-9p-08 | CHM13 15q to 9p | confirmed | definitive multi dataset | 6 of 62 | 6/6 | 6/6 | 6/6 |
| [HG002#2#chr19_MATERNAL#R](HG002_2_chr19_MATERNAL_R/README.md) | HG002 | 19q | SA-19q-08 | CHM13 22q to 19q | confirmed | definitive multi dataset | 6 of 91 | 6/6 | 6/6 | 6/6 |
| [HG002#1#chr17_PATERNAL#L](HG002_1_chr17_PATERNAL_L/README.md) | HG002 | 17p | SA-17p-02 | CHM13 11p to 17p | confirmed | definitive multi dataset | 3 of 3 | 3/3 | 3/3 | 3/3 |
| [HG004|HG002#2#chr9_MATERNAL#L](HG004__HG002_2_chr9_MATERNAL_L/README.md) | HG004 | 9p | SA-9p-08 | CHM13 15q to 9p | confirmed | definitive multi dataset | 6 of 10 | 6/6 | 6/6 | 6/6 |
| [HG01258|20q|SA-20q-01](HG01258__20q__SA-20q-01/README.md) | HG01258 | 20q | SA-20q-01 | CHM13 3q to 20q | confirmed | definitive | 4 of 4 | 4/4 | 4/4 | 4/4 |
| [HG01258|17p|SA-17p-02](HG01258__17p__SA-17p-02/README.md) | HG01258 | 17p | SA-17p-02 | CHM13 11p to 17p | confirmed | definitive | 5 of 5 | 5/5 | 4/5 | 5/5 |
| [NA20346|2q|SA-2q-04](NA20346__2q__SA-2q-04/README.md) | NA20346 | 2q | SA-2q-04 | CHM13 7p to 2q | confirmed | definitive | 2 of 2 | 0/2 | 2/2 | 2/2 |

The trio claim (HG004, transmitting mother of HG002) uses the same unit model as HG002 9p maternal; own-assembly placement for
HG003/HG004 is against the HG002 Q100 windows (the transmitted haplotype), a proxy, because no parental assembly is used.

## Viewing in IGV

Each claim directory has two IGV sessions (IGV desktop 2.16 or later: File > Open Session):

* `<claim>.chm13.igv_session.xml`: genome = the public CHM13 v2.0 FASTA
  (`https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa`, indexed). The hosted IGV genome `hs1` (T2T
  CHM13v2.0) has the same chromosome names. Two loci open side by side: the donor block on CHM13 D and the recipient end.
  Tracks: read segments split at the switch, whole reads, and the claim model (BED; orange = donor sequence and alleles,
  blue = recipient/proximal).
* `<claim>.own.igv_session.xml`: genome = `<claim>.own_windows.fa` (the person's own assembly end windows at the arm, tip
  frame, cut from the public assembly; SHA-256 in the FASTA header and in `data/catalogue/ends.tsv`).

CRAM files carry `@SQ M5` tags and no embedded reference: decode them with the CHM13 v2.0 analysis-set FASTA
(`samtools view -T chm13v2.0.fa.gz`), or let htslib fetch the sequence by MD5. BAM copies are included for tools without CRAM
support. Read groups name the public source file (`DS` = source URL, index, region, source reference).

## Release assets (full set; not in the repository)

| file | size | SHA-256 | content |
|---|---:|---|---|
| `subtel_read_evidence_v1_1KG_ONT.tar.gz` | 23.7 MB | `b7fe1917993a945980e716e1e5837d0725555e1e7e9cbb06f0935b2b193dd324` | per-claim CRAM/BAM/window FASTA, 1KG_ONT |
| `subtel_read_evidence_v1_HG002_Q100.tar.gz` | 190.3 MB | `a27bab22dac5c1eae2d18ded829fbcdbec286e38b2d6af02c76265ae9f51f74a` | per-claim CRAM/BAM/window FASTA, HG002_Q100 |
| `subtel_read_evidence_v1_HG002_trio.tar.gz` | 16.3 MB | `0e5a08c9ceb9e95a39c5d568648d350d156815444f16442cad5a55635517f505` | per-claim CRAM/BAM/window FASTA, HG002_trio |
| `subtel_read_evidence_v1_tables.tar.gz` | 0.7 MB | `1de8d53fcec15780edf15faaaeda0292c08668f14d29ccec0a9e86406883111b` | claims.tsv, reads.tsv.gz, controls.tsv, criteria |

Each `subtel_read_evidence_v1_<cohort>.tar.gz` unpacks to `read_evidence_v1/<cohort>/<claim>/` with the same file set as here
(CHM13 CRAM + BAM, own-window BAM + FASTA, reads TSV), for every claim with at least one v1-conclusive read. Every read can
also be fetched from its public source with the `samtools_command` column of [../reads.tsv.gz](../reads.tsv.gz), and every
verdict re-derived with [../replicate_claim.sh](../replicate_claim.sh).
