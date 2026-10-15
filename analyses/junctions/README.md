# Junctions track — how were foreign terminal blocks created? (hypothesis H2)

A *foreign* chromosome end carries a terminal block that matches another CHM13 end (the **donor**, D) better than its own
arm (the **recipient**, R); form `R>D`. Input calls come from the subtelomere placement of 11,998 chromosome-named
telomeric tips from 459 HPRC Release 2 assemblies (`results/subtel_foreign_tips.tsv` of the analysis home): 4,594 foreign
tips in 347 forms; 1,810 tips have a base-resolved switch and 1,756 of these fall in 383 per-form junction clusters.
This track asks three questions: is the transfer reciprocal or one-way? What sequence lies at the junction? Which junctions
best illustrate the structure at base level?

## Reproduce

```bash
# once: work/config.local.sh (git-ignored) defines ATLAS_INPUT_DIR (analysis home with tipdup/, subtel/, results/, ref/, hg002/)
bash analyses/junctions/run.sh
```
Requirements: python >= 3.10 with pandas, numpy, scipy, matplotlib, mappy (2.31), edlib; minimap2 2.31 and samtools on PATH.
Runtime about 15 min on 3 threads, plus the 1.77 GB HG002 download and two genome-wide minimap2 runs (~15 GB RAM each).
All random draws use seed 20260715 (`scripts/jcommon.py`); `scripts/checksums.py` writes `results/MANIFEST.tsv` (md5, sha256).

| script | step | outputs (`results/`) |
|---|---|---|
| `scripts/reciprocal.py` | 1: donor end in the same haplotype; breakpoint complementarity; lost recipient sequence at haplotype ends | `reciprocal_per_carrier.tsv`, `reciprocal_per_form.tsv`, `reciprocal_breakpoint_complement.tsv`, `reciprocal_summary.json` |
| `scripts/lost_sequence_hg002.py` | 1 (extension): genome-wide search of the replaced recipient sequence in the complete HG002 diploid | `hg002_lost_sequence_calls.tsv`, `hg002_lost_sequence_loci.tsv` |
| `scripts/junction_features.py` | 2: homology, microhomology and telomere-like repeats at junctions vs matched random positions; replaced length | `junction_features.tsv`, `junction_null_summary.tsv`, `junction_lost_sequence.tsv`, `junction_summary.json`, `junction_windows.fa` |
| `scripts/exemplars.py` | 3: three base-level exemplars with read support | `exemplars.tsv`, `exemplar_<id>_alignment.txt`, `exemplar_<id>_track.tsv` |
| `scripts/figure.py` | figure | `figures/junctions_mechanism.{png,pdf}`, `figures/junctions_mechanism_caption.md` |

Column definitions are in each script's docstring.

## Data and provenance

* HPRC Release 2 assemblies (s3://human-pangenomics, index `meta/r2_assemblies_index.csv` of the analysis home); tip windows
  (20 kb; 100 kb for 5,840 tips) and placement tables were produced upstream (analysis home `scripts/subtel_foreign_*.py`).
* CHM13 v2.0 (`ref/chm13v2.0.fa.gz`; outer 500 kb of each end in `subtel/foreign/ref/chm13_ends500k.fa`). Nurk et al. 2022, doi:10.1126/science.abj6987.
* HG002 Q100 v1.1 diploid assembly, s3://human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
  (1,769,801,343 bytes; Last-Modified 2024-06-23; S3 ETag `06d48c014805443c6c3bf5f6110993b3-211`; downloaded 2026-10-15 by
  `scripts/fetch_hg002_genome.sh` into `work/junctions/`). HG002 200 kb ends, the 600 kb chr16_PATERNAL p window and the ONT
  UL / PacBio Revio junction read tables (`results/subtel_hg002_junction_reads.tsv`) come from the analysis home.
* HPRC Release 2: Liao et al. 2023 (doi:10.1038/s41586-023-05896-x) describes the release-1 design; release 2 data are public at the URL above.
* DOIs above were checked against the CrossRef API on 2026-10-15 ("The complete sequence of a human genome", Nurk, 2022;
  "A draft human pangenome reference", Liao, 2023).

## Results

### 1. One-way transfer: the donor end is not affected by carrying the form (`reciprocal_*`)

For the 3,917 carrier tips of the 104 forms with >= 10 carriers (230 samples), the donor arm's end in the **same haplotype
assembly** was classified. Donor ends were assessable for 2,464 carriers (1,418 haplotypes had no named telomeric donor end;
35 unassessed). The distribution of donor-end states is the same in carriers and in haplotypes whose recipient end is
assessed but does not carry the form (pooled over forms):

| donor end in the same haplotype | carriers (n = 2,464) | non-carriers (n = 17,814) |
|---|---|---|
| native (own/tie at 20 kb, no recipient-like window) | 0.4696 | 0.4745 |
| foreign towards a third arm | 0.3973 | 0.3967 |
| contains a recipient-like 5 kb window | 0.0946 | 0.0898 |
| reciprocal (donor end is `D>R`) | 0.0386 | 0.0391 |

Chi-square on the four categories p = 0.88. Per form, no form shows reciprocal enrichment after Bonferroni correction
(one form, 10q>16q, at nominal p = 0.035 with 1/8 carriers), and no form shows depletion of native donor ends or enrichment of
absent donor ends after correction. Absence of a named donor end is similar in carriers (1,418/3,917) and non-carriers
(8,870/26,889). The terminal block of a carrier is not a recent copy of the same haplotype's donor end: median identity of
the block to the same-haplotype donor tip 0.947 (2,063 carriers; 14 at >= 99.9%), i.e. the copies have diverged since transfer.

Because the two ends are on different chromosomes, co-occurrence in one haplotype decays quickly after a reciprocal event;
complementary breakpoints were therefore also tested directly. For 233 junction clusters, 155 have at least one tip of the
reverse form `D>R` in the release; 30 clusters have reverse-form tips whose own-side and foreign-side junction coordinates are
both within 2 kb of the cluster's donor/recipient coordinates (`reciprocal_breakpoint_complement.tsv`). The largest is
20q>3q#3 (121 samples; 20q 31,464 / 3q 34,580): all 6 3q>20q tips in the release carry the complementary breakpoint.
Other complementary pairs at appreciable frequency: 5q>6q#0 (19 samples) with 6q>5q (13 samples), 6q>5q#1 (13) with
5q>6q (19), 18p>10p#1 (12) with 10p>18p (5), 10p>18p#1 (5) with 18p>10p (12). These are candidate reciprocal products (or
switches at a shared homology boundary); no null for chance coordinate agreement was computed, so they are reported as candidates.

Is the recipient's replaced sequence elsewhere in the same haplotype? At Release 2 chromosome ends (named tips, unnamed
telomeric tips, 20 kb ends of all contigs), the CHM13 recipient segment (telomere array end to the junction, 2,146 carriers
with a located switch) has a hit covering >= 50% at >= 98% identity in 761, but only 49 exceed the identity of the best
other CHM13 end (paralogue background) by > 0.05. Placement-based, 2,279/3,917 carriers (58.2%) have another named end with
a recipient-like window, versus 17,004/26,889 non-carrier haplotypes (63.2%): recipient-like windows at other ends are
paralogous background, not a moved copy. In the complete HG002 diploid (genome-wide, 18 foreign ends with a located switch;
`hg002_lost_sequence_calls.tsv`): 14 absent, 2 at paralogue level, 2 present specifically — chr9_PATERNAL q (9q>11p;
the replaced CHM13 9q segment, 90,083 bp, is at the chr16_PATERNAL p end at 99.95% identity, the end that carries the 16p>9q
block) and chr14_PATERNAL p (14p>15p; found at the chr15_PATERNAL p end, which is 15p>14p).

**Verdict (H2a, reciprocal exchange as the general mechanism): negated** — donor-end state is independent of carrier status
(chi-square p = 0.88; reciprocal 3.86% vs 3.91%), the replaced recipient sequence is generally not present elsewhere above
paralogue background (Release 2 ends 49/2,146; HG002 genome 2/18), consistent with one-way, duplicative transfer
(break-induced replication / gene conversion-like terminal copying or non-reciprocal translocation with loss of the recipient
end). A minority of clusters have complementary reverse-form alleles (30/233), including two HG002 paternal configurations;
these are candidate reciprocal events.

### 2. Junctions lie in donor-recipient homology (`junction_*`)

One representative per junction cluster (383; recomputed changepoint equals the published switch for 383/383):

* identical parental sequence at the breakpoint (CHM13 donor vs recipient, extended in both directions from the junction):
  median 85 bp; 0 bp in 33, 1-4 bp in 23, 5-25 bp in 66, 26-200 bp in 156, > 200 bp in 105 junctions. Random position
  pairs in the same two parental ends: 0 bp 56.1%, 1-4 bp 43.2%, 5-25 bp 0.62%, > 25 bp 0.001%. Mean 210.3 vs 0.69,
  permutation p = 0.0005 (2,000 draws; the minimum attainable).
* donor-recipient identity over +-500 bp (edlib global): median 0.974 vs 0.517 at random pairs (p = 0.005, 200 draws).
* 327/383 (85%) junctions are *homology-spanning* (both parental ends align collinearly >= 500 bp on both sides of the
  switch); 22 blunt, 19 microhomology (1-25 bp) and 15 short-homology junctions are non-homologous at 500 bp scale.
  Among clusters shared by >= 5 samples 60/70 are homology-spanning. After merging junctions that appear under two
  paralogous donor labels (218 physical junctions), 180 are homology-spanning and 68.4% have > 25 bp parental identity.
  The most frequent non-homologous junction is blunt: 4q>15p#1 (42 samples; same recipient coordinate as 4q>4p#0).
* interstitial telomere-like arrays (>= 3 in-phase TTAGGG/CCCTAA copies, <= 1 mismatch each) within +-1 kb: 10.97% of
  junctions in the carrier vs 1.35% at random positions (p = 0.0005), but 30 of the 42 carrier arrays are at junctions < 5 kb
  from the telomere. For the 302 junctions >= 5 kb from the telomere in carrier and both parents: 0.99% vs 1.45% (p = 0.82);
  parents 0.99% vs 1.95% (p = 0.94).
* own terminal sequence replaced (recipient distance at the switch): base-resolved tips median 14,543 bp (IQR 10,001-31,464;
  range 587-202,179; n = 1,756); cluster representatives median 12,835 bp; window-resolved tips median 31,310 bp (n = 1,261);
  1,064 tips carry foreign sequence through the whole window (lower bound median 91,449 bp). The donor block length
  (donor distance at the switch) median 15,171 bp; net length change median +216 bp (IQR -2,935 to +1,572).

**Verdicts:** homology at the junction — **confirmed** (mh_parental > 25 bp in 261/383 vs 0.001% expected; p = 0.0005);
microhomology-mediated / blunt end-joining as the main mechanism — **negated** (only 56/383 junctions are non-homologous at the 500 bp scale: 22 blunt, 19 with 1-25 bp microhomology, 15 with short > 25 bp homology);
ITS enrichment at junctions — **negated** beyond the telomere-proximal 5 kb (p = 0.82; the all-junction excess is a
telomere-proximity effect).

### 3. Exemplars (`exemplars.tsv`, `exemplar_*_alignment.txt`, `exemplar_*_track.tsv`)

| id | carrier | junction (carrier) | donor / recipient coordinate (telomere frame) | parental identity at switch | read support |
|---|---|---|---|---|---|
| E1 20q>3q | HG02647.pat CM086569.1 R end (cluster 20q>3q#3, 121 samples; 141 members at the modal coordinate) | contig 67,726,444 (tip 39,041) | 3q 34,580 / 20q 31,463 | 875 bp | Verkko junction concordant 5/5 tested in the cluster (selected tip concordant); hifiasm layout tested in 1 member (10 reads across +-1 kb) |
| E2 3q>11p | HG002 chr3_MATERNAL q | contig 201,063,742 (tip 34,220) | 11p 33,294 / 3q 37,872 | 28 bp | 42 ONT UL and 21 Revio primary alignments cross switch +-5 kb, none clipped at the switch |
| E3 16p>9q | HG002 chr16_PATERNAL p | contig 209,331 (tip 209,330) | 9q 211,735 / 16p 26,773 | 4 bp | 54 ONT UL and 22 Revio crossing alignments, none clipped |

E1 and E2 switch within long collinear donor-recipient homology (identity to the recipient 0.993-0.989 already on the donor
side). E3 is a short-homology junction (4 bp shared, recipient identity 0.51 over the 5 kb before the switch) at which a
~206 kb 9q-like block replaced 26.8 kb of 16p; in the same paternal haplotype the native 9q terminal sequence is at this 16p end
(section 1). The alignment files show donor, carrier and recipient rows over +-150 bp with diagnostic sites marked.

## Caveats

* Donor-end classification rests on 20 kb (6,158 tips) or 100 kb (5,840 tips) windows; donor ends that were never assembled
  telomere-to-telomere (36% of carriers) are not assessed. Window labels inherit the CHM13 48-end placement and its paralogy:
  related donor ends (3q/11p/19p/15q, 4p/15p, 1p/6q) are often near-equivalent, so 'replaced by a third arm' includes
  paralogous variants of the same block.
* The Release 2 lost-sequence search only covers chromosome ends; the genome-wide test is limited to HG002 (18 ends).
* Breakpoint complementarity has no explicit null model; matches may reflect shared homology boundaries.
* Base-resolved junctions are 1,810/4,594 foreign tips; window-resolved and 'through' tips are not in the junction feature
  tests. Junction coordinates are CHM13-relative.
* ITS detection uses a strict in-phase tandem definition; degenerate or short interstitial repeats beyond 1 mismatch per copy
  are not counted.
