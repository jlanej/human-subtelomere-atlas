# Read-level evidence in independent 1000 Genomes ONT reads

Test of catalogue claims of non-native terminal blocks in single long reads that were not used to build the assemblies. The
reads come from the 1000 Genomes ONT Vienna panel: R9.4.1 reads, Guppy 6.2.1 sup basecalls, aligned to CHM13 v2.0 with minimap2
2.26 (Schloissnig, Pani et al. 2025, Nature, doi:10.1038/s41586-025-09290-7; preprint doi:10.1101/2024.04.18.590093). Per-sample
CRAMs are public at `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/<SAMPLE>.t2t.cram`. Only
small regions at chromosome ends were read, by HTTP range requests, with at most 4 concurrent requests. Use follows the panel's
data reuse statement.

The criteria are the frozen v1 criteria ([../CRITERIA.md](../CRITERIA.md), [../criteria_v1.json](../criteria_v1.json), criteria
hash `a89c5eaf023ed032`), unchanged, with platform `ONT_R9.4.1` (amendment C1 telomere-artefact motifs apply). The module is
`atlas/subtel_atlas/readevidence.py`, unchanged; its decision logic was not modified. Calibration on HG002 and on strictly native
controls is in [../calibration/CALIBRATION.md](../calibration/CALIBRATION.md).

## Claims

* **Individuals.** The 42 HPRC Release 2 individuals that are also in the 1KG-ONT panel.
* **Claims.** Every Release 2 end of these individuals with a non-native (`foreign` or `within_family`) structural allele of tier A
  or B:
  * 945 ends, giving 822 claims (individual x arm x allele) on 29 arms and 180 allele IDs;
  * 123 claims carry the same allele on both haplotypes. These are one claim each, with the model from haplotype 1; the second
    end's projection is reported as `switch_concordant_hap2` in the read table.
* **Claim model.** From `readevidence.claims_from_catalogue` (donor D, comparator P, pattern DP/PDP, assembly switch).
* **Other haplotype at the arm** (`other_haplotype` in claims.tsv). Counts:

  | other haplotype | claims |
  |---|---:|
  | same allele | 123 |
  | strictly native | 139 |
  | core-native | 140 |
  | different non-native allele | 369 |
  | not assembled | 40 |
  | several ends | 11 |

  The expected conclusive fraction among conclusive plus contradicting reads (CRITERIA.md section 4) is 1.0 for same-allele
  claims and 0.5 for heterozygotes with a native other haplotype.
* **Controls.** Individuals of the 42 that are strictly native at a claimed arm in both haplotypes (raw structure only
  home/none, criteria change C3). There are 64 sample-arms (33 individuals, 20 arms), and each is tested against the model of
  every allele claimed at that arm (283 control tests). Nine claimed arms have no strictly native individual among the 42: 12p,
  17q, 19p, 19q, 1p, 3q, 4q, 6q and Yq.

## Reads

* **Slices.** Primary alignments (`samtools view -F 0x900`, decoded with CHM13 v2.0) from two regions per claim:
  * the outer 5 kb of the recipient CHM13 end, which catches telomere-reaching reads;
  * the CHM13 donor end at the start of the claimed block (tip coordinates [d_lo - 2 kb, d_lo + 5 kb]), because a carrier read's
    primary alignment can fall on the donor end.

  Regions were merged per sample and end: 1,005 slices, 39,084 records. Reads of at least 5 kb were evaluated (as in
  calibration), deduplicated by read name per claim.
* **Read set of a claim.**
  * Reads with primary alignment in the recipient-end slice.
  * Donor-slice reads only when **conclusive**, that is, anchored at MAPQ 60 on the recipient chromosome. 53 conclusive
    rows came from donor slices.

  Donor-slice reads that are not conclusive are reported in `reads.tsv.gz` but not counted. Of the 2,905 donor-slice reads that
  passed the switch test without being conclusive, 1,093 anchor on the donor chromosome. For example, chr7p reads at the 20p
  claims of SA-20p-03 show CHM13 7p alleles distally and CHM13 20p alleles proximally, while their proximal segment aligns at
  MAPQ 60 to chr7. These reads are reads of the donor arm. Counting them as `consistent` would raise 13 claims from
  not_confirmed to supported and no claim to confirmed (column `verdict_incl_donor_slice_consistent`).
* **Verdicts.** The read verdict is `readevidence.finalize`. Switch concordance is against the claim's unit model applied to the
  person's own end window (C2). The claim verdict is `readevidence.claim_verdict` in its unphased form: confirmed requires at
  least 2 conclusive, switch-concordant reads, and contradicting reads do not block it, because 1KG reads are not phased.

## Results

Verdict classes (claims.tsv `verdict_class`):

| class | claims |
|---|---:|
| confirmed | 43 |
| supported (1 conclusive + >= 1 consistent) | 17 |
| assessed, not confirmed | 123 |
| not assessable | 639 |

Confirmed claims cover 26 individuals and 20 alleles. 37 of 595 tier A claims and 6 of 227 tier B claims were confirmed.

**Not-assessable breakdown.** First applicable reason, `assessability` column:

| reason | claims |
|---|---:|
| `diagnostic_units` | 391 |
| `read_length` | 131 |
| `anchor` | 93 |
| `no_telomere_read` | 19 |
| `no_reads` | 4 |
| `model_window_undefined` | 1 |

* **`diagnostic_units` (391).** The person's own assembly window, treated as an error-free read, does not pass the v1 switch
  test, and no read is conclusive. First failing condition (`asm_projection_fail`):

  | condition | claims |
  |---|---:|
  | fewer than 10 D units | 161 |
  | D purity < 0.80 | 118 |
  | fewer than 10 P units | 97 |
  | D-vs-R check | 13 |
  | P purity < 0.80 | 2 |

  Typical cases are short or highly D/P-similar blocks, and donor-like haplotypes that differ from CHM13 D at many sites.
* **`read_length` (131).** No telomere-reaching read spans the switch plus a 10 kb anchor.
* **`anchor` (93).** Long enough reads exist but none aligns at MAPQ 60 on the recipient chromosome. The outer 100-200 kb of 3q,
  4q and 19p and the outer ~100 kb of 9q are shared with other CHM13 ends, so these arms are affected most.
* **`model_window_undefined` (1).** HG03521 20q SA-20q-05: the raw donor blocks are not collinear on CHM13 11p (d_lo > d_hi),
  so the v1 D window is empty.

**Assessed, not confirmed (123).** At least one conclusive read exists, or telomere-reaching, anchored reads that are long enough
exist:

| reason | claims |
|---|---:|
| `switch_test` | 43 |
| `contradicting_only` | 38 |
| `one_conclusive` | 18 |
| `discordant` | 17 |
| `integrity` | 7 |

The 17 supported claims are `one_conclusive` claims and are counted under supported, not here.

**Discordant claims (17).** These have conclusive reads, but no read switch is concordant with the assembly projection; the
offsets are 2.4-31 kb. 10 of them are SA-3q-01 (donor 11p, pattern DP): 19 conclusive reads in 10 individuals place the
11p-to-3q switch at 3.5-13.7 kb from the telomere, while the assembly projections place it at 30-42 kb. All 19 anchors are in the
3q region that CHM13 shares with other ends (`mapq60_shared`). Own-assembly placement shows that these reads come from 6p, 20q or 9q
ends that carry an 11p-like tip followed by 3q-like sequence (SA-6p-01, SA-6p-118, SA-20q-01, SA-9q-14). Their switches match
those windows, and the SA-3q-01 assemblies are correct ([DIAGNOSTICS.md](DIAGNOSTICS.md), flag A). None of the 31 SA-3q-01 claims is
confirmed.

**Per allele** ([alleles.tsv](alleles.tsv)):

| allele | arm | confirmed | supported | carriers tested |
|---|---|---:|---:|---:|
| SA-20q-01 | 20q | 8 | 5 | 29 |
| SA-2q-02 | 2q | 8 | 2 | 14 |
| SA-17p-02 | 17p | 4 | 1 | 14 |
| SA-19p-04 | 19p | 3 | 2 | 15 |
| SA-19q-02 | 19q | 3 | 1 | 8 |
| SA-19q-04 | 19q | 2 | 1 | 4 |
| SA-4q-05 | 4q | 2 | 0 | 3 |
| SA-17q-02 | 17q | 1 | 0 | 18 |

12 further alleles have 1 confirmed carrier each (20 alleles with at least one confirmed carrier). alleles.tsv has 182 rows (180 allele IDs on 29 arms); for 111 of them every tested
carrier is not assessable.

**Controls.** No conclusive read among the strictly native controls:
* 0 of 307 distinct eligible control reads (telomere, anchor and integrity pass); exact one-sided 95% upper bound 0.97%;
* 0 of 1,010 eligible read x allele-model tests;
* 2,124 distinct control reads in total.

753 control rows were `contradicts`, which is the expected call for native reads. Control coverage is thin at some arms with
confirmed claims: 20q has 3 control individuals with 2 eligible reads, and 2q has 1 individual. 19p, 19q, 4q, 12p and 17q have no
strictly native control among the 42. The HG002 and 1KG specificity in CALIBRATION.md covers other arms.

**Heterozygotes.**
* Claims with a native other haplotype and at least one conclusive read (36): the median conclusive fraction is 1.0 (mean 0.76),
  above the 0.5 expected from allele dosage. Conclusive and contradicting reads have different sensitivities: `contradicts` needs
  at least 10 donor-interval units, `conclusive` needs the full switch test. The fraction is therefore not a calibrated dosage
  estimate.
* Of the 7 confirmed same-allele (homozygous) claims, none has a contradicting read.
* 8 of 123 same-allele claims have contradicting reads, 1-5 each. Examples: HG00133 1p SA-1p-37 (4 reads); NA18960 5q SA-5q-02
  (5 contradicting, 1 conclusive). These are flagged for review.

**Contradiction-heavy claims.** Some claims have contradicting reads only. SA-20p-03 (donor 7p) has 0 conclusive reads in 15
claims and 32 contradicting reads, for example:
* HG00235 20p: 18 contradicting, unique-anchor reads; the other haplotype is not assembled;
* HG00146 20p: 9 contradicting reads; the other haplotype is strictly native;
* NA18960 20p: 2 contradicting reads; homozygous claim.

Unphased reads cannot give a `contradicted` verdict under v1, so these remain `assessed_not_confirmed` (`contradicting_only`). [DIAGNOSTICS.md](DIAGNOSTICS.md) (flag B) traces these
reads to the other haplotype, or to other ends that carry 20p-derived blocks; none places on the claimed window. The allele is real
(57 of 58 ends agree with an independent assembler). The v1 model for it is not appropriate: its assembly switch is the edge of the
donor window, and native 7p sequence passes the same switch test. These claims should be read as not assessable.

**Diagnostics.** [DIAGNOSTICS.md](DIAGNOSTICS.md) places every conclusive and contradicting read on the individual's full set of
Release 2 end windows. Recommendations per claim are in [diagnostics/claim_annotations.tsv](diagnostics/claim_annotations.tsv);
verdicts in this table are unchanged.
* Of the 43 confirmed claims, 27 are robust (>= 2 allele-specific reads with a unique anchor that place on the claimed arm) and 11
  rest on anchors in CHM13-shared sequence but place on the claimed arm.
* 5 are not independently supported: HG00128, NA20503 and NA21110 19p SA-19p-04; HG00323 15q SA-15q-02; HG00235 20q SA-20q-01.
* No claim warrants a `contradicted` annotation: of 117 contradicting reads in the 38 contradicting-only claims, 2 place on the
  claimed end.
* 18 claims (14 of them SA-20p-03) have a donor-window-edge switch and are not assessable under v1.

**Allele resolution within an individual** (secondary, not part of the claim rule):
* `n_conclusive_shared_with_other_claim`. When one individual carries two different non-native alleles at one arm with the same
  donor, one read can be conclusive for both models. HG00235 20q SA-20q-01 is confirmed only by reads that are also conclusive for
  SA-20q-06 (`allele_specific_confirmed` = False). 42 of 43 confirmed claims have at least 2 allele-specific reads.
* Own-assembly placement. Conclusive reads were re-aligned (minimap2 map-ont) to the individual's own Release 2 end windows at
  the arm, both haplotypes:
  * 116 of 144 conclusive concordant reads of confirmed claims have their primary alignment on a claimed end
    (`n_conclusive_own_on_claimed_end`); 40 also meet the strict HG002 definition (no supplementary, >= 95% aligned);
  * in 7 confirmed claims fewer than 2 reads place on the claimed end (`own_other_ends`). Examples: HG03816 17p SA-17p-02, where
    5 reads place at MAPQ 60 on the haplotype 1 end, whose raw structure has a short PF03 block (cf. C3); HG00235 20q SA-20q-06;
    and NA20503 and NA21110 19p, where the reads place on a 3q window of the same individual that shares the outer sequence.

  35 claims are confirmed with at least 2 allele-specific reads and at least 2 own-assembly placements on the claimed end.

## Files

* [claims.tsv](claims.tsv): one row per claim. Includes:
  * identifiers: model end, other haplotype, D/P/R, pattern;
  * switches: assembly switch and projection, median read switch CI;
  * regions and CRAM URL;
  * read counts by stage;
  * verdict, `verdict_class`, `assessability`, `asm_projection_fail`;
  * secondary allele-resolution and own-assembly columns.
* [reads.tsv](reads.tsv): reads that reach the telomere or are not non-informative, in claim read sets and controls. Columns
  include read id, sample, CRAM URL, slice region, read group (RG), CRAM alignment, read length, orientation, terminal telomere,
  unit counts and segments, switch and 95% CI (read subtelomere and CHM13 P coordinates), concordance, anchor, integrity and
  verdict.
* `reads.tsv.gz`: every evaluated row (77,333; carrier rows include the excluded donor-slice reads), with the
  genome-wide alignments. Not in git because of its size (8.8 MB). It ships as `read_evidence_v1/population_reads_all.tsv.gz`
  in the `subtel_read_evidence_v1_tables.tar.gz` release asset (checksums in [../release_manifest.tsv](../release_manifest.tsv)),
  and `python evidence/scripts/population_report.py` rewrites it here from a population run.
* [controls.tsv](controls.tsv): one row per control sample-arm x allele model.
* [alleles.tsv](alleles.tsv): per-allele aggregates.
* Evidence BAMs (git-ignored, regenerable): `work/evidence/population/bam/<sample>.chm13.bam` and `<sample>.own.bam` hold the 210
  distinct conclusive reads (36 individuals; 72 BAMs, MANIFEST.tsv):
  * re-aligned to CHM13 v2.0 and to the individual's own Release 2 end windows (tip frame, telomere at position 1);
  * read group `ID/SM` = sample, `DS` = source CRAM URL;
  * tags `XC:Z` (claim ids) and `XS:Z` (source slice region).
* [DIAGNOSTICS.md](DIAGNOSTICS.md) and [diagnostics/](diagnostics/): root-cause analysis of the SA-3q-01 and SA-20p-03 flags,
  read placement on own-assembly windows, the donor-window-edge screen and per-claim annotations.
* Scripts:
  * [../scripts/run_population.py](../scripts/run_population.py) (plan, fetch, evaluate);
  * [../scripts/population_report.py](../scripts/population_report.py) (report, bams);
  * [../scripts/population_diagnostics.py](../scripts/population_diagnostics.py) (diagnostics);
  * [../scripts/sanitize_bam_headers.py](../scripts/sanitize_bam_headers.py) (BAM copies without local paths in headers).

## Reproduce

```
source work/config.local.sh            # ATLAS_INPUT_DIR: CHM13 v2.0 FASTA and the 48 500 kb CHM13 ends
export PYTHONPATH=atlas
python evidence/scripts/run_population.py all -j 4 -t 5     # plan, fetch, evaluate, report, bams (~3.5 h on 5 cores)
python evidence/scripts/run_population.py report            # again after bams: adds the own-assembly columns
evidence/scripts/replicate_population_claim.sh 'HG00099|20q|SA-20q-01'   # one claim, from public sources
```

The run needs the end windows `work/windows/ends.fa` (pipeline) and a minimap2 `map-ont` index of CHM13 v2.0 (default
`work/readevidence/chm13v2.map-ont.mmi`).

## Limits

* **Platform.** R9.4.1 reads have a higher error rate than R10/HiFi. The fixed error parameter e = 0.05 and the 0.80 purity are
  those calibrated on R10, HiFi and R9 controls. Telomere detection relies on the C1 artefact motifs.
* **Phasing.** Reads are unphased. The claim rule cannot use contradicting reads, and two alleles of one individual at one arm
  with a shared donor are not separable by reads alone (see allele resolution).
* **Read set.**
  * Only reads whose primary alignment lies in the recipient outer 5 kb, or in the donor interval, were fetched. Reads whose
    primary alignment is on a third end, or on P when P is not R (190 of the 945 end models), were not fetched.
  * Reads under 5 kb were not evaluated.
* **Controls.** Strictly native controls exist only within the 42 individuals and cover 20 of 29 claimed arms.
* **Sampling.** Different cell-line passages or DNA preparations of the panel and of HPRC are not distinguished. A claim that is
  not confirmed is not evidence against the assembly.
* **Diagnostic units.** 391 of 822 claims (48%) cannot be decided by any read under v1, because the diagnostic units do not
  separate D from P on the person's own haplotype.
