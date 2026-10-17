#!/usr/bin/env bash
# Re-derive one population claim verdict (evidence/population/claims.tsv) from public sources: the 1000 Genomes ONT Vienna CRAM
# of the individual (IGSR FTP, HTTP range requests), CHM13 v2.0 and the atlas catalogue, with the frozen criteria v1.
# Usage: evidence/scripts/replicate_population_claim.sh '<sample>|<arm>|<allele>'   e.g. 'HG00099|20q|SA-20q-01'
# Needs: samtools >= 1.16, minimap2, python with mappy, pysam, numpy, pandas, scipy; CHM13 v2.0 (bgzip + .fai/.gzi) and the 48
# 500 kb CHM13 ends (subtel-atlas fetch-refs) under ATLAS_INPUT_DIR (or --chm13/--chm13-ends), the end windows
# work/windows/ends.fa (pipeline/catalogue_common.py) and a minimap2 map-ont index or FASTA of CHM13 v2.0 (--genome).
set -euo pipefail
here="$(cd "$(dirname "$0")/../.." && pwd)"
claim="${1:?claim id <sample>|<arm>|<allele>}"; shift || true
[ -f "$here/work/config.local.sh" ] && source "$here/work/config.local.sh"
export PYTHONPATH="$here/atlas${PYTHONPATH:+:$PYTHONPATH}"
tag="$(echo "$claim" | tr '|' '_')"
work="$here/work/evidence/replicate/$tag"
run() { python "$here/evidence/scripts/run_population.py" "$1" --work "$work" --out "$work/out" "${@:2}"; }
run plan --only "$claim" --no-controls "$@"
run fetch -j 4 "$@"
run evaluate "$@"
run report "$@"
python - "$work/out/claims.tsv" "$here/evidence/population/claims.tsv" "$claim" <<'EOF'
import csv, sys
new = {r['claim_id']: r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')}[sys.argv[3]]
ref = {r['claim_id']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}.get(sys.argv[3], {})
keys = ['verdict', 'assessability', 'n_reads', 'n_telomere_reaching', 'n_conclusive_concordant', 'n_contradicts', 'n_consistent']
for k in keys: print(f'{k:26s} replicate={new.get(k)!s:22s} published={ref.get(k)}')
sys.exit(0 if all(new.get(k) == ref.get(k) for k in keys[:1]) else 1)
EOF
