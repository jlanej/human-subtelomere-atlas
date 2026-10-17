# HG002#2#chr3_MATERNAL#R

**Claim.** HG002 3q (maternal), structural allele SA-3q-01 (foreign):
the telomere-adjacent block is closer to CHM13 11p (donor D) than to CHM13 3q; the assembly places the switch back
to CHM13 3q 36,576 bp from the telomere boundary (CHM13 chr3:201064202).

**Verdicts.** v1 (pre-registered): **confirmed**, 7 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive multi dataset**, 3 definitive reads
(datasets: ont_lc24_ulk,ont_q28,ont_ul_dorado). This directory holds 3 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): 197 SNV sites in 1 homologous block(s) plus private 200 bp bins.

## Files

| file | content |
|---|---|
| `HG002_2_chr3_MATERNAL_R.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `HG002_2_chr3_MATERNAL_R.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `HG002_2_chr3_MATERNAL_R.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `HG002_2_chr3_MATERNAL_R.own.bam`, `HG002_2_chr3_MATERNAL_R.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `HG002_2_chr3_MATERNAL_R.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `HG002_2_chr3_MATERNAL_R.chm13.igv_session.xml` | IGV: two loci, donor block chr11:1-47170 and recipient end chr3:200589768-201105948 |
| `HG002_2_chr3_MATERNAL_R.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 11p for 3 of 3 reads (MAPQ 60-60).
   The proximal segment aligns to CHM13 3q for 3 of 3 reads, 3 of them at MAPQ 60. The reads therefore join
   11p-like sequence at the telomere to 3q sequence, switching inside CHM13 chr3:201062020-201069434
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 11p sequence and alleles, blue features 3q. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 3 of 3 reads place best on the claimed end window `HG002_2_chr3_MATERNAL_R` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
