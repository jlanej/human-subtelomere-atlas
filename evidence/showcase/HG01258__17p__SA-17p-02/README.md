# HG01258|17p|SA-17p-02

**Claim.** HG01258 17p (hap1), structural allele SA-17p-02 (foreign):
the telomere-adjacent block is closer to CHM13 11p (donor D) than to CHM13 17p; the assembly places the switch back
to CHM13 17p 16,713 bp from the telomere boundary (CHM13 chr17:9825).

**Verdicts.** v1 (pre-registered): **confirmed**, 5 conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **definitive**, 5 definitive reads
(datasets: 1KG_ONT_HG01258). This directory holds 5 of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): 12 SNV sites in 1 homologous block(s) plus private 200 bp bins.

## Files

| file | content |
|---|---|
| `HG01258__17p__SA-17p-02.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `HG01258__17p__SA-17p-02.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `HG01258__17p__SA-17p-02.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `HG01258__17p__SA-17p-02.own.bam`, `HG01258__17p__SA-17p-02.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `HG01258__17p__SA-17p-02.diagnostic_sites.chm13.bed` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `HG01258__17p__SA-17p-02.chm13.igv_session.xml` | IGV: two loci, donor block chr11:192372-227028 and recipient end chr17:1-60001 |
| `HG01258__17p__SA-17p-02.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 11p for 5 of 5 reads (MAPQ 60-60).
   The proximal segment aligns to CHM13 17p for 4 of 5 reads, 4 of them at MAPQ 60. The reads therefore join
   11p-like sequence at the telomere to 17p sequence, switching inside CHM13 chr17:15670-15775
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are 11p sequence and alleles, blue features 17p. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** 5 of 5 reads place best on the claimed end window `HG01258_1_CM088120.1_L` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
