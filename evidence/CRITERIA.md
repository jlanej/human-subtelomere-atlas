# Read-level evidence for non-native terminal blocks: pre-registered criteria (v1)

Status: **pre-registered** on 2026-10-16 and **frozen after calibration** the same day (training half of HG002 and native
controls only; changes A1, C1-C3 below, no numeric threshold changed), before any read-level data (HG002 Q100, 1000 Genomes ONT, CHM13 reads) were fetched
or inspected. The thresholds are frozen in [criteria_v1.json](criteria_v1.json) (a byte-identical copy of the registered file is
kept as [criteria_v1_preregistered.json](criteria_v1_preregistered.json); the package ships a copy as
`atlas/subtel_atlas/data/criteria_v1.json`). The only information used to set them was CHM13 v2.0 itself: donor-recipient
homology, the attrition of diagnostic sites under the filters below, and the shared-sequence masks of the recipient ends.
Calibration ([calibration/CALIBRATION.md](calibration/CALIBRATION.md)) may change a threshold only on HG002 training ends and
native negative controls. Any change is recorded with its reason in `criteria_v1.json` (`calibration_changes`) before the
held-out HG002 ends or any carrier population reads are evaluated.

Implementation: `atlas/subtel_atlas/readevidence.py`, command `subtel-atlas reads-evidence`. It needs only CHM13 v2.0 (the 48
500 kb ends from `subtel-atlas fetch-refs`, plus the full genome or a minimap2 index for the anchor) and the reads. Assemblies
are used only for the secondary concordance line.

## 1. The claim

A claim concerns one chromosome end of recipient arm **R** (CHM13 end `chrR`). It is derived from the catalogue
(`readevidence.claims_from_catalogue`) or given explicitly:

* **D** (donor): the CHM13 end that best matches the claimed non-native block: the best end covering most bases of the outermost
  non-home core block (`mosaic_blocks.tsv`, raw blocks inside it).
* **P** (proximal comparator): CHM13 R when the next core block is home (or there is none), otherwise the dominant best end of
  the next core block.
* **pattern**: `DP` = terminal block (telomere, D, then P); `PDP` = internal block (telomere, P, D, P).
* **switch**: assembly switch in subtelomere coordinates (`core_switches_sub`) and, where defined, projected onto CHM13 P
  (`p_switch` = `chm13_dist_start` of the first raw block of the next core block when its best end is P).
* D window = [d_lo - 10 kb, d_hi + 10 kb] on D (d_lo/d_hi = span of the block on D); P window = [0, max(p_switch, switch,
  d_hi) + 30 kb] on P (300 kb if the switch lies beyond the assembly window).

If P differs from R, a second unit set (D vs R) is built. The read's distal segment must also favour D over R
(section 3c).

## 2. Diagnostic units (reference-only; independent of reads and assemblies)

1. **SNV sites.** Align the D window to the P window (minimap2 `asm20`, forward strand). Keep blocks of at least 2,000 bp at
   gap-compressed identity of at least 0.90 (greedy one-to-one). Inside them, every single-base difference is a candidate. A
   candidate is excluded if:
   * another difference (SNV or indel edge) lies within 3 bp;
   * either base touches a homopolymer of 3 or more bases in D or P;
   * the ±25 bp context has fewer than 20 distinct trinucleotides (low complexity);
   * the 16-mer centred on the site occurs again (or as its reverse complement) within ±2 kb in D or P (tandem or local repeat);
   * the site lies within 500 bp of the CHM13 telomere.
   Small indels are not used in v1.
2. **Private bins.** D-window (P-window) positions not covered by any such block are D-private (P-private). Tile the read
   into 200 bp bins. A bin is a D unit if at least 90% of its D projection is D-private, its identity to D is at least 0.85, and
   its identity to D exceeds its identity to P by at least 0.05. P units are defined symmetrically. Bins inside homologous blocks
   contribute only through SNV sites, so no difference is counted twice.

## 3. Per-read evidence fields

(a) **Terminal telomere.** The read is oriented telomere-first (tip frame) by the shared caller `telomere.find_block` (both
orientations tried). Terminal = a (CCCTAA)n block of at least **T = 300 bp** starting within **E = 300 bp** of the read end.
`s0_read` = offset + block length; read subtelomere coordinate = read position minus `s0_read`.

(b) **Allele calls.** The read is aligned separately to the D window and to the P window (mappy, `map-ont` or `map-hifi`
preset, forward-strand hits of at least 1 kb, greedy non-overlapping). At each SNV site the read base is taken from both
projections:
* both defined and equal, or only one defined: call **D**, **P** or **other** by comparison with the two CHM13 alleles;
* both defined but different: **other**.
Counts per read: D, P, other, uncovered.

(c) **Switch localisation.** Units (SNV and bin) are ordered by read position. Each unit matches its segment's expected label
with probability 1 - e (e = 0.05, fixed). The maximum-likelihood segmentation for the pattern (`DP`: one change point; `PDP`:
two) is compared with the no-switch alternatives all-P and all-D:
* log10 likelihood ratios `log10_lr_vs_allP` and `log10_lr_vs_allD`;
* a 95% profile-likelihood interval for each change point: all cut positions within 1.92 log units of the maximum, reported
  as read coordinates between the flanking units;
* the same interval in read subtelomere coordinates and on CHM13 P (P projection of the first proximal unit);
* per-segment unit counts and purity;
* one-sided binomial log10 p-values for the D count in the D segment, and the P count in the P segments, under the error
  rate e.
If P is not R, the D segment must also have at least 10 D-vs-R units with D purity of at least 0.80.

(d) **Anchor.** The proximal part of the read, from the upper interval bound of the last change point to the read end, is
aligned genome-wide to CHM13 v2.0 (minimap2 `map-ont` index, used for all platforms). Anchor = an alignment that is:
* at least **10 kb** long, with MAPQ of at least **60**;
* on the recipient chromosome, in the orientation of a telomere-to-centromere continuation of arm R;
* within 2 Mb of the R telomere.
`anchor_unique_bp` = anchor bases on positions of CHM13 R that no other CHM13 end covers at 95% or more over at least 1 kb
(positions beyond 500 kb count as unique). Anchor tier: `unique` (at least 1,000 such bases) or `mapq60_shared`. Because
CHM13 3q, 4q and 19p share their outer ~170-200 kb, and 9q its outer ~100 kb, with other CHM13 ends at 95% or more
(9p: outer ~45 kb), MAPQ 60 on R is required and the unique tier is
reported, not required.

(e) **Integrity.** All four conditions must hold:
* no internal telomere: no (CCCTAA)n or (TTAGGG)n block of 500 bp or more starting more than 1 kb beyond the terminal block;
* no internal adapter: no ONT ligation/rapid adapter or PacBio SMRTbell sequence at an edit fraction of 0.15 or less, outside
  200 bp from either read end;
* no foldback: no alignment of 1 kb or more at MAPQ 20 or more on the R chromosome on the opposite strand to the anchor;
* length consistent: genome-wide alignments (primary and supplementary), together with the read's D/P/R window alignments
  from (b), cover at least 90% of the read's subtelomeric part (amendment A1, below).

(f) **Verdict per read.**
* **conclusive**: terminal telomere, anchor and integrity all pass, and the switch test passes:
  * every segment of the pattern has at least 10 units at purity of at least 0.80;
  * both log10 LRs are at least 5;
  * if P is not R, the D-vs-R check in (c) passes.
* **contradicts**: terminal telomere, anchor and integrity pass, and in the claimed donor interval the read shows the native
  arm: at least 10 units with P purity of at least 0.80. The claimed donor interval is subtelomere coordinates
  [switch_sub_distal or 0, switch_sub - 2 kb].
* **consistent**: the switch test passes but terminal telomere, anchor or integrity fails (typically a read too short to
  reach the anchor).
* **non-informative**: everything else.

`switch_concordant` (secondary, assembly line) = the read's switch interval, widened by the tolerance (2 kb), overlaps the
assembly switch. The overlap is tested in subtelomere coordinates, or on CHM13 P where `p_switch` is defined.
`assembly_concordant` (secondary) = for reads aligned to the person's own assembly (HG002 Q100 alignments), the primary
alignment is on the claimed contig end. It must also have no supplementary alignment and cover at least 95% of the read.

## 4. Claim-level rule

A claim is **confirmed** when:
* at least **2 conclusive, switch-concordant reads** support it (noted when they come from at least 2 independent datasets or
  platforms), and
* there is **no conclusive contradicting read from the same haplotype-assignable set**: for HG002 Q100 this is the reads whose
  primary alignment is on the claimed haplotype contig; for unphased population reads it is undefined.

The other claim-level calls are:
* **supported**: 1 conclusive read plus at least 1 consistent read;
* **contradicted**: at least 2 contradicting reads and no conclusive read, in a haplotype-assignable set;
* **not confirmed**: otherwise.

Heterozygous carriers: unphased reads from a carrier whose other haplotype is native are expected to include contradicting
reads, so both counts are reported. For an unphased sample the conclusive fraction among conclusive plus contradicting reads
is expected near 0.5 (heterozygous) or 1 (homozygous). A claim on an unphased sample is confirmed by conclusive reads
regardless of contradicting reads, but never when the conclusive reads are not switch-concordant.

## 5. Calibration and specificity design

* **Sensitivity**: HG002 Q100 non-native ends (24 claims; ONT R10 UL dorado and HiFi Revio reads from the public Q100
  alignments, fetched by region). Ends are split into a training and a held-out half by a fixed rule (sorted end_id,
  alternating).
* **Specificity**: every claim applied to reads that should not carry it:
  * HG002 reads of the other haplotype at the same arm where that haplotype is native;
  * HG002 reads of the other haplotype where it carries a different structure (reported separately);
  * 1KG-ONT reads at the same arm from individuals native in both haplotypes per the catalogue (no carrier reads are inspected
    during calibration).
* Reported: per-read error rate at diagnostic SNV sites (mismatch to both alleles, and wrong allele on segments of known state),
  false conclusive rate with an exact binomial 95% upper bound, and sensitivity as a function of read length.

## Amendments (implementation clarifications; no threshold changed)

* **A1 (2026-10-16, before calibration).** The aligned-fraction test in (e) unions the genome-wide alignments with the
  reference-window alignments of (b). Without this, minimap2 drops the donor-block alignment of a carrier read when it overlaps
  the primary alignment on the read (`--secondary=no`), so true carriers fail the test. Found in a smoke test of the
  synthetic unit-test fixture and one HG002 claim (4q paternal, ONT UL; held-out half). The threshold (0.90) is unchanged.

* **C1 (calibration).** ONT R9.4.1 reads: the terminal-telomere caller also accepts the published R9 telomere
  basecalling artefacts in tip frame (CTTCTT, CCCTGG, TTTTAA and their rotations; Tan et al. 2022). For 1KG-ONT native-control
  reads whose CHM13 alignment reaches within 3 kb of the chromosome end, detection rose from 188/883 (21%) to 637/883 (72%).
  This affects ONT R9 only. T and E are unchanged.
* **C2 (calibration).** Switch concordance is tested against the person's assembly end window analysed with the claim's own
  unit model (an error-free "read"). It falls back to the catalogue switch when that projection shows no switch. The catalogue
  switch is a 5 kb sliding-identity boundary with a tie rule, not a diagnostic-unit change point. In the training half,
  conclusive carrier reads at 17p (both haplotypes), 19p paternal and 6q maternal placed the D/P switch unanimously 2-7 kb from
  it. Catalogue concordance is still reported (`switch_concordant_catalogue`).
* **C3 (calibration).** Negative controls must be **strictly native**: raw mosaic structure only `home`/`none` on every
  contributing haplotype. Every conclusive read among the core-native controls came from individuals whose assemblies carry
  short (<10 kb) non-home raw blocks at the arm, so the reads detect those blocks correctly:
  * HG002 19q paternal (`22q|4q|21q|home`);
  * 1KG HG00128 and HG00133 at 19q;
  * HG02976 at 9p.
  Core-native controls are reported separately.
* **T0.** No numeric threshold changed. True-carrier donor-segment purity is as low as 0.80 (HG002 15q paternal). That is
  biological: the carrier's donor-like haplotype differs from CHM13 D at some diagnostic sites. A higher purity threshold would
  discard true carriers.

## References

* Nurk et al. 2022, The complete sequence of a human genome. Science 376:44-53. doi:10.1126/science.abj6987
* Li 2018, Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100. doi:10.1093/bioinformatics/bty191
* Schloissnig, Pani et al. 2025, Structural variation in 1,019 diverse humans based on long-read sequencing. Nature.
  doi:10.1038/s41586-025-09290-7 (preprint doi:10.1101/2024.04.18.590093)
* Tan et al. 2022, Identifying and correcting repeat-calling errors in nanopore sequencing of telomeres. Genome Biology
  23:180. doi:10.1186/s13059-022-02751-6
