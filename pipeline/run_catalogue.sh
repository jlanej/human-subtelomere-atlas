#!/bin/bash
# Build the allele catalogue (release v0.1) from scratch: windows -> CHM13 tracks/mosaics -> arms -> structural alleles -> haplogroups ->
# evidence/tiers -> release tables, QC report and manifest. Idempotent: finished downloads are skipped (.ok markers).
# Requirements: samtools/htslib with libcurl, minimap2 >= 2.24, python >= 3.10 with numpy, pandas, scipy, matplotlib, mappy, pysam.
# Configuration: pipeline/config.sh (+ git-ignored work/config.local.sh: ATLAS_INPUT_DIR, ATLAS_THREADS). Without ATLAS_INPUT_DIR
# the pipeline runs in public mode (pipeline/INPUTS.md; recipe in pipeline/CATALOGUE.md).
# Usage: bash pipeline/run_catalogue.sh
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$ATLAS_REPO"; mkdir -p work/windows work/catalogue data/catalogue figures/qc
PY=python
if [ -z "${ATLAS_INPUT_DIR:-}" ]; then                      # public mode: data/inputs/ + downloads into work/public/
  bash pipeline/fetch_public_inputs.sh chm13 pedigree      # CHM13 v2.0 analysis set (0.98 GB), pedigree ends (92 MB); verified
fi
$PY pipeline/catalogue_01_plan_windows.py                 # window + probe plans (reuses earlier 100 kb fetches if available)
bash pipeline/catalogue_02_fetch_windows.sh               # S3 range reads (Release 2 windows, arm probes, HG002 Q100 ends)
$PY pipeline/catalogue_03_assemble_windows.py build       # oriented windows, provenance, sha256
if [ -s work/catalogue/bins.tsv.gz ] && [ -s work/catalogue/bins_stamp.tsv ]; then
  $PY pipeline/catalogue_04_align_bins.py update           # recompute bins only for ends whose boundary/window changed
else
  $PY pipeline/catalogue_04_align_bins.py                  # 5 kb bins vs 48 CHM13 ends (full)
fi
$PY pipeline/catalogue_04t_align_tracks.py               # one minimap2 pass per end (all chains up to 200 hits, cs) for the tracks
$PY pipeline/catalogue_05_arms.py                         # home arm (names / probes / bin-mosaic fallback)
$PY pipeline/catalogue_06t_mosaic_tracks.py ties          # per-point co-tie counts among CHM13 ends (defining ends)
$PY pipeline/catalogue_06t_mosaic_tracks.py families      # paralog families (tau 0.75) -> data/catalogue/paralog_families.tsv
$PY pipeline/catalogue_06t_mosaic_tracks.py label         # grid labels -> blocks, base-level switches, S5/S10/S20 states
$PY pipeline/catalogue_08_alleles.py                      # structural alleles
$PY pipeline/catalogue_07_haplogroups.py                  # sequence haplogroups
$PY pipeline/catalogue_08b_stable_ids.py                 # keep SA/HG IDs of the previous release (reciprocal best overlap)
$PY pipeline/catalogue_08e_boundary_stable.py            # per-end boundary_stable (allele and S20 at s0 -/+500 bp, track shift)
# Validation of the boundary-shift stability (redesign acceptance; QC section 9). Needs a previous build for the stratified subset
# (ATLAS_PREVIOUS_RELEASE); 12c needs the Tan et al. 2026 cross-check cache in work/tan2026 (pipeline/validation).
if [ -d "${ATLAS_PREVIOUS_RELEASE:-work/catalogue/v0.1_prefix/release}" ]; then
  $PY pipeline/catalogue_12_stability_subset.py           # >= 1,000 stratified ends
  $PY pipeline/catalogue_04t_align_tracks.py shift        # re-alignment at shifts -1000..+1000 bp
  $PY pipeline/catalogue_12_stability_eval.py tau 0.4,0.5,0.75,1.01   # family threshold scan
  $PY pipeline/catalogue_12_stability_eval.py final       # stability table with Wilson intervals
  $PY pipeline/catalogue_12b_offset_consensus.py          # alternative design (offset consensus of bins), 300-end subsample
fi
if [ -s work/tan2026/pairs.tsv ]; then $PY pipeline/catalogue_12c_tan_pairs.py; fi   # same-haplotype pairs, before/after
$PY pipeline/catalogue_09_release.py                      # evidence, tiers, frequencies, release tables
$PY pipeline/catalogue_10_qc.py                           # QC.md, figures/qc, MANIFEST
$PY pipeline/catalogue_11_package_checksums.py            # query-tool package: checksums of this release
$PY tests/regenerate_golden.py                            # query-tool golden outputs (embed catalogue values)
