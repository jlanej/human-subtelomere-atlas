# HG01258|20q|SA-20q-01

**Claim.** HG01258 20q (hap2), structural allele SA-20q-01 (foreign):
the telomere-adjacent block is closer to CHM13 3q (donor D) than to CHM13 20q; the assembly places the switch back
to CHM13 20q 31,422 bp from the telomere boundary (CHM13 chr20:66177732).

**Verdicts.** v1 (pre-registered): **confirmed**, 4 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive**, 4 definitive reads
(datasets: 1KG_ONT_HG01258). This directory holds 4 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): no SNV sites (no homologous block of >= 2 kb at >= 90% between the two windows), so all units are private 200 bp bins: sequence present in CHM13 3q but absent from CHM13 20q, and the reverse.

## Files

| file | content |
|---|---|
| `HG01258__20q__SA-20q-01.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `HG01258__20q__SA-20q-01.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `HG01258__20q__SA-20q-01.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `HG01258__20q__SA-20q-01.own.bam`, `HG01258__20q__SA-20q-01.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `HG01258__20q__SA-20q-01.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `HG01258__20q__SA-20q-01.chm13.igv_session.xml` | IGV: two loci, donor block chr3:201060310-201105948 and recipient end chr20:66120031-66210255 |
| `HG01258__20q__SA-20q-01.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 3q for 4 of 4 reads (MAPQ 60-60).
   The proximal segment aligns to CHM13 20q for 4 of 4 reads, 4 of them at MAPQ 60. The reads therefore join
   3q-like sequence at the telomere to 20q sequence, switching inside CHM13 chr20:66176491-66186091
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 3q sequence and alleles, blue features 20q. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 4 of 4 reads place best on the claimed end window `HG01258_2_CM088150.1_R` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
