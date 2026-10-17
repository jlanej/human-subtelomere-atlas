# Cross-validation against Tan et al. 2026

This directory compares the allele catalogue (`data/catalogue/`, release v0.1, **revision 2**: track-based mosaic with paralog-family
labels; see `data/catalogue/QC.md`, CHANGELOG) with the subtelomere sequences released by Tan et al. 2026 (bioRxiv,
doi:[10.64898/2026.09.03.748401](https://doi.org/10.64898/2026.09.03.748401)). Their sequences are distributed on Zenodo
(doi:[10.5281/zenodo.21807156](https://doi.org/10.5281/zenodo.21807156), archive `subtelomere_500kb_2kbtelo.zip`) under the
Creative Commons Attribution 4.0 licence (CC BY 4.0). We thank the authors for releasing them.

- **No Tan et al. sequence is redistributed here.** The tables contain only record names, coordinates and the calls derived from
  them.
- Tan et al. present their own catalogue (paralogy-block haplotypes S1–S7 over 500 kb), built independently with a different
  method. Their haplotype labels are not public, so the comparison is made at the sequence and structure level.
- Concurrent work: Guarracino et al. 2026, bioRxiv doi:[10.64898/2026.07.10.737660](https://doi.org/10.64898/2026.07.10.737660).

**Reproduce.**
- `bash pipeline/validation/run_tan2026.sh` is the full run. The track alignment of 35,329 records takes about 6 CPU-hours (minimap2).
- `bash pipeline/validation/run_tan2026.sh compare` re-runs only the comparison after a catalogue update. It re-labels Tan
  records from the cached alignments, which happens only if `paralog_families.tsv` changed, and then runs steps 4–7. Only
  pairwise alignments whose catalogue inputs changed are redone.
- The catalogue files used are identified by SHA-256 in `tan2026_summary.json` → `catalogue`.
- Numbers from the same cross-check against revision 1 are in `tan2026_summary.json` → `revision1_vs_revision2`.

## Files

| file | content |
|---|---|
| `tan2026_crosswalk.tsv` | one row per HPRC same-haplotype pair: catalogue end vs Tan record, alignment metrics, the calls on both sides and their agreement |
| `tan2026_summary.json` | headline numbers: catalogue SHA-256, revision 1 vs 2, same-haplotype concordance (overall, by assembler pairing, boundary agreement, tier, nativity class, arm, form, artefact candidates, junk-flagged ends, unnamed-contig arms), cohort replication, window sensitivity, novel structures |
| `tan2026_cohort_forms.tsv` | per cohort and ancestry subgroup: frequency of each tier A/B structural allele (with nativity class) and each common family-level form; 95% Wilson CI; HPRC Release 2 comparator; Fisher exact p; BH q |
| `tan2026_novel_structures.tsv` | structures in the independent cohorts that match no catalogue allele, clustered as in the catalogue; `candidate_allele` = seen in ≥ 2 unrelated samples |
| `../../figures/validation/tan2026_concordance.png` / `.pdf` | summary figure |

## Methods

Scripts are in `pipeline/validation/` (`tan2026_01` … `tan2026_07`; there is no step 03). At most 4 threads are used.

1. **Records** (`tan2026_01_extract.py`).
   - All 35,329 records were streamed from the archive.
   - The telomere block was called with the catalogue caller (`subtel_atlas.telomere`). Records are telomere-first; 7 needed
     reverse-complementing.
   - 4,333 records (12.3%) have no telomere block. Only 58% of CPC and 76% of JASA records carry one. Records without a
     telomere are classified but excluded from the primary analyses.
   - Window = telomere + junk (s0) + the first 100 kb of subtelomere.
2. **Mosaic method, catalogue revision 2** (`tan2026_02_mosaic.py`).
   - Alignment: one minimap2 pass per record of window[s0:] against the 48 CHM13 500 kb ends, using
     `catalogue_04t_align_tracks.run` and `subtel_atlas.tracks.MM2_TRACK_ARGS`. The PAF of the 10,773 records first aligned
     for the catalogue's own Tan check is reused; the other 24,556 records were aligned with the same function.
   - Labels: `subtel_atlas.tracks.mosaic_end`. This gives per-base identity tracks, a 100 bp grid anchored at the record tip,
     paralog-family labels (`paralog_families.tsv`), base-level boundaries, and S5/S10/S20 as the majority grid label.
   - Home arm = the Tan arm label, mapped to the catalogue arm group.
3. **Allele assignment.** `subtel_atlas.alleles.core` and `match_allele`, the matcher shared by the catalogue and the query
   tool: same core structure, every switch within 10 kb of the allele median, nearest median first. Applied to the full-length
   Release 2 ends, it returns each end's own allele for 99.83% of them (98.74% with the legacy most-frequent-first order).
4. **Window emulation.** A catalogue window is 100 kb from the contig tip, i.e. 100 kb minus the telomere (median 92.3 kb of
   subtelomere); Tan records give 100 kb of subtelomere. Tan blocks are clipped before the core-structure computation:
   - same-haplotype pairs: at the paired catalogue end's own subtelomere length;
   - cohort records: at a length drawn deterministically (seeded by the record name) from the arm's Release 2 subtelomere
     lengths.

   This matters because 24 tier-A/B and 46 tier-C alleles have their last switch beyond 80 kb. Such alleles can only be called
   on ends whose window reaches ≥ 10 kb past that switch, so in the catalogue they depend on telomere length. Sensitivity of
   the pooled cohort comparison to the clip (tier A/B, 293 tests):

   | window | q < 0.05 | non-native r | common tier-A alleles observed |
   |---|---|---|---|
   | sampled lengths (primary) | 4 | 0.971 | 88/88 |
   | arm-median clip | 9 | 0.959 | 87/88 |
   | unclipped 100 kb | 21 | 0.801 | 87/88 |
5. **Same-haplotype pairs** (`tan2026_04_pairs.py`).
   - Haplotype numbers do not correspond between the catalogue (mostly hifiasm) and Tan HPRCY2 (mostly Verkko) assemblies outside
     trio-phased samples. Each named catalogue end is therefore paired with the Tan record of the same sample and arm that aligns
     best (mappy asm5, either label).
   - Alignments start 500 bp before each telomere boundary, so identity is subtelomere identity. Catalogue boundaries and
     switches are projected into Tan coordinates through the alignment.
   - *Verified* = identity ≥ 0.99 over ≥ 50% of the subtelomere. *Primary* = verified and the Tan record carries a telomere.
   - S-state agreement = same state and same paralog family. Revision 1 compared state + best single CHM13 end.
6. **Cohorts** (`tan2026_06_cohorts.py`). Inputs are the HGSVC3, CPC, APC and JASA records.
   - Removed:
     - samples in the catalogue or Tan HPRCY2, including aliases: NA19240, and NA24385 (= HG002);
     - trio children: HG00514 (1000 Genomes pedigree) and APC `apr-s`;
     - 53 windows containing assembly gaps.
   - Primary records carry a telomere and ≥ 80 kb of subtelomere.
   - Ancestry comparators:
     - CPC vs EAS;
     - JASA Japanese (1000 Genomes JPT) vs EAS;
     - JASA Saudi (`ksa`) vs ALL only;
     - HGSVC3 per 1000 Genomes superpopulation (IGSR pedigree table `20130606_g1k_3202_samples_ped_population.txt`);
     - APC vs ALL only.
   - Forms are family-level (arm>family of a foreign S20 call). Allele nativity classes come from
     `alleles_structural.tsv` (`nativity_class`).

## Results (catalogue v0.1 revision 2; numbers from `tan2026_summary.json`)

**Same haplotype (HPRC).** 11,262 named Release 2 ends have a same-sample, same-arm Tan HPRCY2 record. 11,108 are verified and
10,999 are primary. Revision 1 had 10,994 primary pairs; the revision-1 telomere-boundary fix changed s0 for some ends. Haplotype
labels are the same for 7,376 verified pairs and opposite for 3,732. Median identity is 0.99999.

| primary pairs | revision 1 | revision 2 |
|---|---|---|
| S5 / S10 / S20 state | 90.1 / 92.2 / 96.4% | 85.3 / 93.8 / 97.8% |
| form | 97.2% | 98.7% |
| core structure | 92.7% | 96.3% |
| structural allele (compatible) | 93.6% (94.1%) | 96.4% (96.7%), n = 10,940 |
| Tan record unassigned | 1.65% | 0.54% |
| allele agreement, tier A / B / C | 96.3 / 81.2 / 72.8% | 98.0 / 85.6 / 74.4% |
| allele agreement by boundary difference: 0 / 1–100 bp / 101 bp–1 kb / > 1 kb | 98.8 / 96.6 / 89.8 / 73.4% | 98.1 / 97.9 / 95.9 / 87.3% |
| base-level switches identical after projection | 98.2% of 6,970 | 95.9% of 10,769 |

- **Boundary sensitivity of alleles is much reduced.** Under revision 2, allele agreement falls only from 98.1% (identical
  boundary call; 36.0% of pairs) to 87.3% (> 1 kb difference).
- **S5 remains boundary-sensitive.** S5 agreement is 99.0% when the boundary call is identical, 89.3% at 1–100 bp and 46.6%
  above 1 kb. Its first 5 kb lies next to the telomere, where the boundary call moves the window.
- By assembler, allele agreement is:
  - hifiasm vs Verkko: 96.2% (9,759 pairs);
  - Verkko vs hifiasm: 96.6% (262);
  - hifiasm vs hifiasm: 97.5% (978; effectively the same assemblies, identity 1.0 in 98.5%).
- By nativity class of the catalogue allele: native 99.6% (7,257), foreign 90.9% (3,296), within-family 83.4% (446). Labels
  inside a paralog family (e.g. 4p vs PF01, 1p/6q vs PF02) are the least reproducible. The lowest arms are 11p (85.7%), 1p, 19p
  and 6q (87–88%), and 15q.
- Of the 398 allele disagreements, 89% involve a different core structure.
- In the independent-assembler pairs, 9,764 base-level switches in 3,992 pairs were compared (figure panel b): 96% fall at the
  identical base after projection. Without projection (each record measured from its own boundary), across all primary
  pairs (10,769 base-level switches, including hifiasm-vs-hifiasm pairs; not computed for the independent-assembler subset
  alone), only 39.3% are identical, but 83.5% are within 100 bp.

**Artefact candidates and junk-flagged ends.**
- Three of the four artefact candidates reproduce exactly (same structure and allele) in Tan's Verkko assemblies: HG02647 16q
  (SA-16q-43), HG03710 4p (SA-4p-01), HG04228 6p (SA-6p-120).
- HG03486 8q (structure PF02|2q|PF02|8p; on a 90 kb contig carrying both chromosome 8 ends) has no aligned Tan 8q partner. Its
  best Tan match is a 6q record, consistent with its flag.
- Junk-flagged ends: 150 primary pairs; S20 98.0%, allele 95.3%.

**Arms of unnamed contigs** (sequence-based; unchanged from revision 1). The catalogue arm agrees with the Tan arm label for 97.6%
of the 6,707 matched unnamed-contig ends. By arm-call confidence: high 99.6% (5,087), medium 96.7% (847), low 85.3% (251),
mosaic_group 85.4% (514). The control, named ends, agrees for 99.55%.

**Independent cohorts.** 11,175 primary records from 182 samples:

| cohort | records | samples |
|---|---|---|
| HGSVC3 | 4,189 | 59 |
| CPC | 2,736 | 58 |
| APC | 3,500 | 52 |
| JASA | 750 | 13 |

- 97.7–98.7% of primary records match a catalogue allele (revision 1: 96.3–98.3%).
- All 88 tier-A alleles with Release 2 frequency ≥ 5% recur in the pooled cohorts, each in ≥ 3 samples.
- 124/129 tier-A and 128/166 tier-B alleles are observed. By nativity class: native 45/50, foreign 182/215, within-family 25/30.
- Frequencies of non-native alleles correlate with Release 2 (Pearson r 0.971, n = 253; revision 1: 0.947, n = 362).
- 14 of 5,573 allele tests and 5 of 589 form tests have BH q < 0.05 (revision 1: 6 of 7,652 and 12 of 1,113). All are against
  the all-ancestry reference; none differs from a matched superpopulation.
  - Example: SA-20q-02 (native 20q) is 0.77 in CPC and 0.08 in HGSVC3 AFR, against 0.46 overall.
  - SA-20q-01 (20q>PF04) is 0.486 pooled (95% CI 0.429–0.544) vs 0.483 in Release 2. In CPC it is 0.214, vs 0.194 in Release 2
    EAS.
  - The pooled native SA-4p-01 (0.965 vs 0.863) and within-family SA-4p-02 (PF01|home; 0.032 vs 0.128) differ. This may be
    ancestry or the within-family label sensitivity noted above.

**Novel structures.**
- 152 primary records (1.4%; revision 1: 2.9%) match no catalogue allele. They form 132 clusters.
- 15 clusters are seen in ≥ 2 unrelated samples (12 new structures, 3 known structures with new switch positions); 2 are seen in
  ≥ 3 samples.
- One of the 15 (PAR1 home|none) contains a block without CHM13 alignment.
- These are candidates only: most come from hifiasm-only cohorts without orthogonal support.

## Limits

- Switches and labels describe which CHM13 end, or paralog family, a segment most resembles inside shared paralogous sequence.
  They are not necessarily breakpoints, and the labels are relative to CHM13.
- The same-haplotype comparison is not fully independent. Tan HPRCY2 hifiasm records are effectively the Release 2 assemblies,
  and the Verkko assemblies use the same sequencing data.
- Records without a telomere (12.3%; 42% of CPC) are not assessed. Cohort frequencies count both haplotypes of each sample.
  Relatedness within CPC, APC and JASA beyond the labelled APC trio is unknown. JASA (13 samples) has wide intervals.
- Alleles whose last switch lies within ~10 kb of the window end depend on telomere length in the catalogue itself. S5 states
  remain sensitive to the telomere-boundary call.
