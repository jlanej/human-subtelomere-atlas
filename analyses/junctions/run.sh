#!/bin/bash
# Reproduce the junctions track (H2: mechanism that created foreign terminal blocks).
# Requires: pipeline/config.sh with ATLAS_INPUT_DIR (analysis home holding tipdup/, subtel/, results/, ref/, hg002/) set in
# work/config.local.sh; python >= 3.10 with pandas, numpy, scipy, matplotlib, mappy, edlib; minimap2 and samtools on PATH.
# Deterministic: fixed random seed (scripts/jcommon.py SEED); runtime ~15 min + HG002 download (1.77 GB) and genome alignment.
set -euo pipefail
cd "$(dirname "$0")"
. ../../pipeline/config.sh
python scripts/reciprocal.py          # step 1: donor-end status in the same haplotype, breakpoint complementarity, lost-sequence search
python scripts/junction_features.py   # step 2: microhomology / homology / ITS at base-resolved junctions vs matched random positions
python scripts/exemplars.py           # step 3: exemplar alignments, identity tracks, read support
bash scripts/fetch_hg002_genome.sh    # public HG002 Q100 v1.1 assembly -> work/junctions/
python scripts/lost_sequence_hg002.py # step 1 extension: genome-wide search of replaced recipient sequence in HG002
python scripts/figure.py
python scripts/checksums.py
