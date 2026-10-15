#!/bin/bash
# Reproduce the telomere_cis track (H1: telomere length is set in cis by the subtelomere allele).
# Requirements: ATLAS_INPUT_DIR set in work/config.local.sh (see pipeline/config.sh); python with numpy, pandas, scipy,
# statsmodels, pysam, matplotlib, pypdf, pillow; minimap2 and samtools on PATH; the 100 kb tip windows of the catalogue
# pipeline (work/windows/fetch, pipeline/catalogue_02_fetch_windows.sh). Runtime ~30-35 min with ATLAS_THREADS=3.
set -euo pipefail
cd "$(dirname "$0")"
. ../../pipeline/config.sh
S=scripts
python $S/00_karimian_digitize.py          # data/karimian2024_fig4e.tsv (skipped if the source PDF is absent)
python $S/01_baseline.py                   # arm baseline, truncation diagnostics, comparison with Karimian et al. 2024
python $S/01b_fixed_window_states.py       # telomere-length-independent allele states (S5/S10/S20) + replica of published t20
for t in S20 t20 S10 S5; do python $S/02_paired.py $t; done   # paired allele tests (S20 primary)
python $S/02b_window_robustness.py       # window artefact checks: boundary QC, W5_15, concordant pairs, TL_ext, strata
python $S/03_hg002_reads.py                # HG002 per-read telomere length (ONT UL, HiFi Revio, HiFi DC)
python $S/04_hg002_tests.py                # HG002 tests
python $S/05_figures.py                    # figures + robustness summary
python $S/06_manifest.py                   # MANIFEST.tsv (md5/sha256) + PROVENANCE.json
