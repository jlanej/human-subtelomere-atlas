# HG002#1#chr17_PATERNAL#L

**Claim.** HG002 17p (paternal), structural allele SA-17p-02 (foreign):
the telomere-adjacent block is closer to CHM13 11p (donor D) than to CHM13 17p; the assembly places the switch back
to CHM13 17p 20,251 bp from the telomere boundary (CHM13 chr17:5582).

**Verdicts.** v1 (pre-registered): **confirmed**, 7 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive multi dataset**, 3 definitive reads
(datasets: ont_duplex,ont_ul_dorado). This directory holds 3 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): no SNV sites (no homologous block of >= 2 kb at >= 90% between the two windows), so all units are private 200 bp bins: sequence present in CHM13 11p but absent from CHM13 17p, and the reverse.

## Files

| file | content |
|---|---|
| `HG002_1_chr17_PATERNAL_L.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `HG002_1_chr17_PATERNAL_L.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `HG002_1_chr17_PATERNAL_L.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `HG002_1_chr17_PATERNAL_L.own.bam`, `HG002_1_chr17_PATERNAL_L.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `HG002_1_chr17_PATERNAL_L.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `HG002_1_chr17_PATERNAL_L.chm13.igv_session.xml` | IGV: two loci, donor block chr11:192105-227761 and recipient end chr17:1-82945 |
| `HG002_1_chr17_PATERNAL_L.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 11p for 3 of 3 reads (MAPQ 60-60).
   The proximal segment aligns to CHM13 17p for 3 of 3 reads, 3 of them at MAPQ 60. The reads therefore join
   11p-like sequence at the telomere to 17p sequence, switching inside CHM13 chr17:5943-6143
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 11p sequence and alleles, blue features 17p. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 3 of 3 reads place best on the claimed end window `HG002_1_chr17_PATERNAL_L` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
