#!/bin/bash
# Build the allele catalogue (release v0.1) from scratch: windows -> CHM13 mosaics -> arms -> structural alleles -> haplogroups ->
# evidence/tiers -> release tables, QC report and manifest. Idempotent: finished downloads are skipped (.ok markers).
# Requirements: samtools/htslib with libcurl, minimap2 >= 2.24, python >= 3.10 with numpy, pandas, scipy, matplotlib, mappy, pysam.
# Configuration: pipeline/config.sh (+ git-ignored work/config.local.sh: ATLAS_INPUT_DIR, ATLAS_THREADS).
# Usage: bash pipeline/run_catalogue.sh
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$ATLAS_REPO"; mkdir -p work/windows work/catalogue data/catalogue figures/qc
PY=python
$PY pipeline/catalogue_01_plan_windows.py                 # window + probe plans (reuses earlier 100 kb fetches)
bash pipeline/catalogue_02_fetch_windows.sh               # S3 range reads (Release 2 windows, arm probes, HG002 Q100 ends)
$PY pipeline/catalogue_03_assemble_windows.py build       # oriented windows, provenance, sha256
[ -s work/catalogue/bins.tsv.gz ] || $PY pipeline/catalogue_04_align_bins.py   # 5 kb bins vs 48 CHM13 ends
$PY pipeline/catalogue_05_arms.py                         # home arm (names / probes / mosaic fallback)
$PY pipeline/catalogue_06_mosaic.py                       # blocks + base-level switches
$PY pipeline/catalogue_08_alleles.py                      # structural alleles
$PY pipeline/catalogue_07_haplogroups.py                  # sequence haplogroups
$PY pipeline/catalogue_09_release.py                      # evidence, tiers, frequencies, release tables
$PY pipeline/catalogue_10_qc.py                           # QC.md, figures/qc, MANIFEST
