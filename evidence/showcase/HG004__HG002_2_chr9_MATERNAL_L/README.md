# HG004|HG002#2#chr9_MATERNAL#L

**Claim.** HG004 9p (transmitting_parent (HG002 maternal)), structural allele SA-9p-08 (foreign):
the telomere-adjacent block is closer to CHM13 15q (donor D) than to CHM13 9p; the assembly places the switch back
to CHM13 9p 19,636 bp from the telomere boundary (CHM13 chr9:23460).

**Verdicts.** v1 (pre-registered): **confirmed**, 13 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive multi dataset**, 10 definitive reads
(datasets: HG004_ont_r10_PAO98462_c3dbcd08,HG004_ont_ul_r941). This directory holds 6 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): 100 SNV sites in 1 homologous block(s) plus private 200 bp bins.

## Files

| file | content |
|---|---|
| `HG004__HG002_2_chr9_MATERNAL_L.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `HG004__HG002_2_chr9_MATERNAL_L.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `HG004__HG002_2_chr9_MATERNAL_L.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `HG004__HG002_2_chr9_MATERNAL_L.own.bam`, `HG004__HG002_2_chr9_MATERNAL_L.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `HG004__HG002_2_chr9_MATERNAL_L.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `HG004__HG002_2_chr9_MATERNAL_L.chm13.igv_session.xml` | IGV: two loci, donor block chr15:99720713-99753195 and recipient end chr9:1-126532 |
| `HG004__HG002_2_chr9_MATERNAL_L.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 15q for 6 of 6 reads (MAPQ 46-60).
   The proximal segment aligns to CHM13 9p for 6 of 6 reads, 6 of them at MAPQ 60. The reads therefore join
   15q-like sequence at the telomere to 9p sequence, switching inside CHM13 chr9:24032-25155
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 15q sequence and alleles, blue features 9p. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 6 of 6 reads place best on the claimed end window `HG002_2_chr9_MATERNAL_L` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
