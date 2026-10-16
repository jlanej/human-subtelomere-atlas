# Allele catalogue pipeline

`bash pipeline/run_catalogue.sh` rebuilds `data/catalogue/` (release v0.1) from public assemblies and the analysis inputs
named by `ATLAS_INPUT_DIR` (set in the git-ignored `work/config.local.sh`; see `pipeline/config.sh`). Large intermediates are
written to `work/` (git-ignored). Requirements: samtools/htslib with libcurl, minimap2, Python 3 with numpy, pandas, scipy,
matplotlib, mappy and pysam.

The telomere caller, track alignment and labelling, paralog families, bin tiling (arm inference), arm, structure and haplogroup
functions are implemented once in the query tool package (`atlas/subtel_atlas`; tracks in `tracks.py`), which the pipeline imports
from the source tree; no installation is needed. The query tool labels a window with the same alignment settings, the published
paralog families and the same labelling code as the catalogue.

| step | script | what it does |
|---|---|---|
| 1a | catalogue_01_plan_windows.py | 100 kb telomere-anchored window per Release 2 telomeric tip (20,144 tips; named and unnamed contigs); reuses earlier 100 kb range reads; plans 10 kb arm probes at 250 kb and 1 Mb for unnamed-contig tips and 1,000 random named tips (calibration) |
| 1b | catalogue_02_fetch_windows.sh | range reads from the public HPRC S3 bucket (samtools faidx on the bgzipped FASTA; 4 parallel jobs, 3 retries, .ok markers); HG002 Q100 v1.1 ends |
| 1c | catalogue_03_assemble_windows.py | orients every window telomere -> centromere, adds CHM13 v2.0 (48 ends), HG002 Q100 (92) and WashU pedigree ends (368, validation only), calls every telomere boundary on its own window (shared caller), extends windows to boundary + 80 kb where the boundary lies beyond 20 kb (range reads), records provenance and SHA-256, checks the first 20 kb against the earlier tip windows |
| 2a | catalogue_04_align_bins.py | 5 kb bins tiled from the telomere/subtelomere boundary aligned to the 48 CHM13 500 kb ends (minimap2 asm20); used only for the arm-inference fallback of unnamed contigs (step 2b); `update` mode recomputes only ends whose boundary or window changed (bins_stamp.tsv) |
| 2a' | catalogue_04t_align_tracks.py | one alignment of each whole subtelomere window (window[s0:]) to the 48 CHM13 ends (minimap2 asm20 -c, all chains up to 200 hits with secondary-to-primary score ratio 0.01, cs tags); `shift` mode re-aligns the stability subset from s0 + d (d = -1000 ... +1000 bp) |
| 2b | catalogue_05_arms.py | home arm: contig name, else probe mapping to CHM13 (calibrated), else bin-mosaic rule; pooled groups acroP/PAR1/PAR2 |
| 2c | catalogue_06t_mosaic_tracks.py | `ties`: base-level identity tracks per CHM13 end (5 kb windows sampled every 100 bp, in the end's own coordinates and therefore independent of the tiling phase); points where several ends tie (within 0.5 percentage points of the best) counted over the defining ends. `families`: two CHM13 ends are linked when they co-tie at >= 200 points and at >= 75% of the points where the rarer of the two is tied; connected ends form a paralog family (`paralog_families.tsv`). `label`: each point labelled home (own arm or group among the tied ends), a family or a single end, or none (< 80% identity); runs < 2 kb merged into a neighbour; blocks; base-level switch refinement; S5/S10/S20 states = majority label over the first 5/10/20 kb |
| - | catalogue_06_mosaic.py | previous (bin-based) mosaic; superseded by 2c, kept for the comparison in QC section 9 |
| 2d | catalogue_08_alleles.py | core structures (blocks >= 10 kb), structural alleles by structure + switch positions (complete linkage, 10 kb), assignment of partial / low-confidence / pedigree ends (nearest median switches, then most defining ends) |
| 2e | catalogue_07_haplogroups.py | within-arm all-vs-all k-mer ANI (FracMinHash) of the first 80 kb of subtelomere, UPGMA at ANI 0.995, medoids, alignment validation |
| 2f | catalogue_08b_stable_ids.py | carries SA/HG IDs over from the previous build by reciprocal best member overlap; writes id_history.tsv |
| 2g | catalogue_08e_boundary_stable.py | per-end `boundary_stable`: allele and S20 state with the boundary moved by -/+500 bp, from the end's own tracks shifted by d (approximation; agreement with full re-alignment reported in QC section 9) |
| 2h | catalogue_12_stability_subset.py, catalogue_12_stability_eval.py | acceptance test: seeded stratified subset (>= 1,000 ends; tier x nativity x duplicon-rich arms) re-aligned at boundary shifts of -/+100, 500, 1000 bp; core-structure, S20 and allele stability with Wilson 95% intervals; family threshold scan (QC section 9) |
| 2i | catalogue_12b_offset_consensus.py | alternative design evaluated for comparison: previous bin method with consensus over five tiling offsets, 300-end subsample |
| 2j | catalogue_12c_tan_pairs.py | same-haplotype pairs with the Tan et al. 2026 records (cached cross-check in work/tan2026): allele, structure and S20 agreement before and after the redesign |
| - | catalogue_08c_boundary_shift.py, catalogue_08d_shift_control.py | previous boundary-sensitivity step (bin method); superseded by 2g |
| 3a | catalogue_09_release.py | evidence (Verkko, read layout, HG002 reads, pedigree), flags, tiers, frequencies by superpopulation, release tables, representatives, JSON |
| 3b | catalogue_10_qc.py | QC.md, figures/qc, PROVENANCE.json, MANIFEST.tsv, SHA256SUMS |
| 3c | catalogue_11_package_checksums.py, tests/regenerate_golden.py | checksums of this release in the query-tool package; query-tool golden outputs |

## Rebuild from public data

No access to the original analysis home is needed. With `ATLAS_INPUT_DIR` unset the pipeline runs in public mode: the small
derived inputs come from the committed snapshot `data/inputs/` (tip list, assembly index and catalogue, prior evidence tables,
index files; 1.8 MB) and everything large is fetched from public URLs and checked against recorded SHA-256
(`pipeline/INPUTS.md` lists every input with its source, upstream script, size and checksum).

```bash
cd human-subtelomere-atlas
unset ATLAS_INPUT_DIR; export ATLAS_THREADS=6         # optional: work/config.local.sh
bash pipeline/fetch_public_inputs.sh chm13 pedigree verify
bash pipeline/run_catalogue.sh                        # calls fetch_public_inputs.sh itself when ATLAS_INPUT_DIR is unset
python tests/test_catalogue_integrity.py && python -m pytest tests/ -q
```

Measured on a 10-core laptop (6 threads, about 5 MB/s from the public buckets); full-catalogue figures marked ~ are estimates
scaled from the subset run, except step 04t, which was measured on the full build (34 min at 8 threads).

| stage | subset run (arms 20q, 3q, 16p; 1,339 ends) | full catalogue (20,652 ends) |
|---|---|---|
| `fetch_public_inputs.sh chm13 pedigree verify` (0.98 GB CHM13, 92 MB pedigree ends, 348 verification windows) | 6.9 min | 6.9 min |
| 01-02 range reads (windows + arm probes from 455 assemblies, 4 parallel jobs) | 31 min | ~1.5-3 h |
| 03 windows, boundaries | 37 s | ~5 min |
| 04 bins (arm-inference fallback) | 108 s | ~25 min |
| 04t track alignment | 72 s | 34 min (measured, 8 threads) |
| 05 arms (incl. CHM13 asm5 index) | 61 s | ~5 min |
| 06t ties + families + label | 6 s (label only) | ~10 min |
| 07 haplogroups | 90 s | ~20 min |
| 08, 08b, 08e, 09-11 | < 1 min | ~15 min |
| 12-12c validation (optional; 12c needs the Tan et al. cache) | - | ~1 h |

Verification of the public path (`data/catalogue/qc/`):

* `public_window_verification.tsv` - 300 released windows (seeded sample stratified by set; Release 2 named / unnamed contigs,
  with and without extended windows; HG002 Q100; pedigree) plus all 48 CHM13 ends (306 distinct ends) re-fetched from their
  public URLs with `pipeline/fetch_end_window.py`: 348/348 identical to `window_sha256`. The CHM13 analysis-set download is
  identical to the file used for the release (SHA-256 in `data/inputs/public_files.tsv`), and the pedigree ends re-fetched by
  `pipeline/fetch_public_pedigree.py` are byte-identical to the original files.
* `public_rebuild_subset_check.tsv` / `.json` - a copy of the repository without `work/` and with `ATLAS_INPUT_DIR` unset ran
  steps 01, 02, 03, 04, 04t, 05, 06t label, 08, 07 and 08b for every end whose released catalogue arm is 20q, 3q or 16p
  (`ATLAS_SUBSET_ENDS`; 1,306 Release 2, 24 pedigree, 6 HG002 Q100 and 3 CHM13 ends). All 48,204 compared values (36 columns:
  window coordinates and SHA-256, telomere boundary, arm calls, profile, raw and core structures, switches, S5/S10/S20 states,
  families and best ends, structural-allele and haplogroup IDs and assignments) are identical to the release
  (`pipeline/verify_public_rebuild.py`); the 89 structural alleles and 43 haplogroups of these arms keep their released IDs.
  A subset run uses the published `paralog_families.tsv` (families are built from all defining ends by `06t ties/families` in a
  full run).

What a public rebuild inherits rather than recomputes (details in `pipeline/INPUTS.md`): the list of 20,144 telomeric tips from
the original scan of 459 Release 2 assembly files (`tipdup/tips_meta.tsv`; the catalogue re-calls every boundary on its own
window), and the prior evidence tables (Verkko concordance, read-layout support, HG002 read support, pedigree transmissions,
tip placement) used for evidence columns and flags.

Interpretation note: switches (junctions) mark where the sequence changes which CHM13 end it most resembles. Inside shared,
paralogous subtelomeric sequence this is a change in paralog similarity, not necessarily a structural breakpoint. 'Foreign'
is distance-dependent and is reported along the subtelomere (`profile`) and as boundary-anchored S5/S10/S20 states. Ends of
one paralog family are near-identical over most of their shared sequence, so a block is labelled with the family rather than an
individual end; `*_best_end` columns keep the single best CHM13 end for continuity.

References: CHM13 v2.0 - Nurk et al. 2022, doi:10.1126/science.abj6987. minimap2 - Li 2018, doi:10.1093/bioinformatics/bty191.
Mash distance - Ondov et al. 2016, doi:10.1186/s13059-016-0997-x.

Shared code: the telomere caller, track alignment/labelling, paralog families, bin tiling/summary (arm inference), structure/allele matching and sketch functions live in `atlas/subtel_atlas/` and are used by both this pipeline and the query tool.
