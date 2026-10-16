#!/usr/bin/env bash
# Cross-validation of the allele catalogue against Tan et al. 2026 (bioRxiv doi:10.64898/2026.09.03.748401).
# Input: subtelomere_500kb_2kbtelo.zip from Zenodo doi:10.5281/zenodo.21807156 (CC BY 4.0), placed at
#        work/tan2026/subtelomere_500kb_2kbtelo.zip (or set TAN_ZIP). The catalogue (data/catalogue/), its windows
#        (work/windows/) and CHM13 end reference (work/catalogue/chm13_ends500k.fa) must exist (pipeline/run_catalogue.sh).
# Threads: ATLAS_THREADS (capped at 4). Intermediates: work/tan2026/ (git-ignored). Outputs: data/validation/, figures/validation/.
# Usage: bash pipeline/validation/run_tan2026.sh [all|compare]. Steps are resumable/cached; 'compare' re-labels the Tan records from
# the cached track alignments and re-runs steps 4-7 (no re-alignment of Tan records against CHM13) after a catalogue update.
set -euo pipefail
cd "$(dirname "$0")/../.."
[ -f work/config.local.sh ] && source work/config.local.sh
export ATLAS_THREADS=${ATLAS_THREADS:-4}
mkdir -p work/tan2026
ZIP=${TAN_ZIP:-work/tan2026/subtelomere_500kb_2kbtelo.zip}
[ -f "$ZIP" ] || { echo "missing $ZIP: download subtelomere_500kb_2kbtelo.zip from https://doi.org/10.5281/zenodo.21807156"; exit 1; }
# 1000 Genomes pedigree/population table (IGSR), used for HGSVC3/JASA ancestry and trio children
PED=work/tan2026/igsr_3202_ped.txt
[ -s "$PED" ] || curl -fsS -o "$PED" https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
MODE=${1:-all}   # all | compare (compare = re-label from the cached Tan alignments and re-run the comparison; no re-alignment)
if [ "$MODE" = all ]; then
  [ -f work/tan2026/tan_windows.tsv ] || python pipeline/validation/tan2026_01_extract.py   # telomere-anchored windows
  python pipeline/validation/tan2026_02_mosaic.py          # track alignment (cached PAF) + labels with the published paralog families
else
  python pipeline/validation/tan2026_02_mosaic.py label    # labels only (recomputed if paralog_families.tsv changed)
fi
python pipeline/validation/tan2026_04_pairs.py     # same-haplotype alignments (cached; only pairs with changed catalogue inputs re-aligned)
python pipeline/validation/tan2026_05_compare.py   # per-record classification + crosswalk + summaries
python pipeline/validation/tan2026_06_cohorts.py   # independent cohorts, novel structures
python pipeline/validation/tan2026_07_report.py    # summary JSON + figure
