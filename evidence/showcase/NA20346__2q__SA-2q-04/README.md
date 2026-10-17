# NA20346|2q|SA-2q-04

**Claim.** NA20346 2q (hap2), structural allele SA-2q-04 (foreign):
the telomere-adjacent block is closer to CHM13 7p (donor D) than to CHM13 2q; the assembly places the switch back
to CHM13 2q 12,773 bp from the telomere boundary (CHM13 chr2:242681777).

**Verdicts.** v1 (pre-registered): **confirmed**, 3 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive**, 2 definitive reads
(datasets: 1KG_ONT_NA20346). This directory holds 2 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): 43 SNV sites in 1 homologous block(s) plus private 200 bp bins.

## Files

| file | content |
|---|---|
| `NA20346__2q__SA-2q-04.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `NA20346__2q__SA-2q-04.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `NA20346__2q__SA-2q-04.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `NA20346__2q__SA-2q-04.own.bam`, `NA20346__2q__SA-2q-04.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `NA20346__2q__SA-2q-04.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `NA20346__2q__SA-2q-04.chm13.igv_session.xml` | IGV: two loci, donor block chr7:1-20696 and recipient end chr2:242581293-242696752 |
| `NA20346__2q__SA-2q-04.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 7p for 0 of 2 reads; elsewhere: chr5 2.
   The proximal segment aligns to CHM13 2q for 2 of 2 reads, 2 of them at MAPQ 60. The reads therefore join
   7p-like sequence at the telomere to 2q sequence, switching inside CHM13 chr2:242689560-242690543
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 7p sequence and alleles, blue features 2q. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 2 of 2 reads place best on the claimed end window `NA20346_2_JBIRDS010000056.1_R` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
