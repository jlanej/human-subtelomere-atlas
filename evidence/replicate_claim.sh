#!/usr/bin/env bash
# Re-derive one read-evidence verdict (evidence/claims.tsv or evidence/controls.tsv) from public sources only.
#
#   evidence/replicate_claim.sh 'HG002#2#chr9_MATERNAL#L'          # HG002 Q100 claim (maternal 9p)
#   evidence/replicate_claim.sh 'HG01258|20q|SA-20q-01'            # 1000 Genomes ONT claim
#   evidence/replicate_claim.sh 'control|HG00097|20q|SA-20q-01'    # strictly native control
#   evidence/replicate_claim.sh --list                             # claims with cited reads, by definitive verdict
#
# Needs: internet; samtools >= 1.16 and minimap2 on PATH; python >= 3.9 with numpy, pandas, scipy, pysam, mappy (the package
# atlas/ of this repository is used directly). First run downloads CHM13 v2.0 (~0.98 GB) into the cache
# (~/.cache/subtel-atlas/replicate, or SUBTEL_REPLICATE_CACHE); the in-memory minimap2 index needs ~8 GB RAM.
# Optional: SUBTEL_CHM13 (local CHM13 v2.0 analysis-set .fa.gz with .fai/.gzi), SUBTEL_CHM13_MMI (prebuilt map-ont index),
# SUBTEL_GRCH38 (only for HG003/HG004 ONT CRAMs). Extra arguments are passed to evidence/scripts/replicate_claim.py
# (e.g. -t 8, --no-placement). Exit status 0 = verdicts identical to the published tables.
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$here/atlas${PYTHONPATH:+:$PYTHONPATH}"
py="${PYTHON:-python3}"
for t in samtools minimap2; do command -v "$t" >/dev/null || { echo "missing: $t (conda install -c bioconda $t)" >&2; exit 2; }; done
if [ "${1:-}" = "--list" ]; then
  "$py" - "$here/evidence/claims.tsv" <<'EOF'
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
by = collections.defaultdict(list)
for r in rows:
    if int(r['n_conclusive'] or 0) or r['v1_verdict'] in ('confirmed', 'supported'): by[r['definitive_verdict']].append(r)
for k in sorted(by):
    print(f'## {k} ({len(by[k])})')
    for r in by[k]: print(f"  {r['claim_id']:45s} {r['cohort']:11s} {r['allele_id']:12s} v1={r['v1_verdict']:13s} definitive_reads={r['n_definitive']}")
EOF
  exit 0
fi
exec "$py" "$here/evidence/scripts/replicate_claim.py" "$@"
