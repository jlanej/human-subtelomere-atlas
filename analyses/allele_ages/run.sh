#!/bin/bash
# One-command reproduction of the allele_ages track (population distribution, single origin, relative age, apes).
# Requires: python3 with pandas, numpy, scipy, statsmodels, matplotlib, pysam, mappy; samtools (htslib with https) on PATH.
# Input location: ATLAS_INPUT_DIR = analysis home holding results/subtel_foreign_tips.tsv, results/subtel_foreign_forms.tsv,
#   tipdup/tips_meta.tsv, meta/r2_assemblies_index.csv, subtel/foreign/r2_100k/, subtel/foreign/ref/chm13_ends500k.fa.
#   Set it in the environment or in the git-ignored work/config.local.sh (export ATLAS_INPUT_DIR=...; export ATLAS_THREADS=3).
# Network: raw.githubusercontent.com, ftp.1000genomes.ebi.ac.uk, human-pangenomics.s3.amazonaws.com, api.ncbi.nlm.nih.gov,
#   eutils.ncbi.nlm.nih.gov. Runtime ~1.5 h with 3 threads (S3 range reads ~25 min, alignments ~25 min).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
[ -f "$REPO/work/config.local.sh" ] && source "$REPO/work/config.local.sh"
: "${ATLAS_INPUT_DIR:?set ATLAS_INPUT_DIR (see header)}"
cd "$REPO"
S=analyses/allele_ages/scripts
# 0. sample metadata (reuses cached raw downloads in data/metadata/raw if present; delete them to re-download)
if [ -s data/metadata/raw/SHA256SUMS ]; then python $S/build_samples_tsv.py --offline; else python $S/build_samples_tsv.py; fi
# 1. population distribution
python $S/pop_distribution.py
python $S/fig_population.py
# 2. single origin and age: fetch 100 kb windows of remaining named tips, align, summarise
python $S/fetch_native_100k.py
bash $S/fetch_native_100k.sh
python $S/origin_age.py --min-carriers 19
python $S/origin_summary.py
# 3. great apes
python $S/fetch_ape_ends.py
python $S/ape_forms.py
python $S/ape_junction_check.py
python $S/origin_summary.py   # re-run so the ape summary is included
python $S/fig_origin.py
# checksums of released result tables and figures
( cd analyses/allele_ages && shasum -a 256 results/*.tsv results/*.json figures/*.png figures/*.pdf > SHA256SUMS )
echo "allele_ages done"
