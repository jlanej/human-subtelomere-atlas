#!/bin/bash
# Step 1b: range-read planned windows from the public HPRC Release 2 bgzipped assemblies (htslib reads the remote .fai/.gzi),
# plus the 92 HG002 Q100 v1.1 ends. Polite: ATLAS_FETCH_JOBS (default 4) assemblies at a time, 3 retries, idempotent (.ok markers).
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$ATLAS_REPO/work/windows/fetch"
P="${ATLAS_FETCH_JOBS:-4}"
cat jobs.txt | xargs -P "$P" -L 1 sh -c '[ -s "$0.ok" ] && exit 0; for t in 1 2 3; do samtools faidx -r "$0.regions" "$1" > "$0.fa" 2> "$0.log" && [ "$(grep -c ">" "$0.fa")" -eq "$(wc -l < "$0.regions")" ] && { date -u +%FT%TZ > "$0.ok"; exit 0; }; sleep 10; done; echo "FAILED $0" >> failed.txt'
# arm-inference probes (same polite scheme)
cd "$ATLAS_REPO/work/windows/probes"
cat jobs.txt | xargs -P "$P" -L 1 sh -c '[ -s "$0.ok" ] && exit 0; for t in 1 2 3; do samtools faidx -r "$0.regions" "$1" > "$0.fa" 2> "$0.log" && [ "$(grep -c ">" "$0.fa")" -eq "$(wc -l < "$0.regions")" ] && { date -u +%FT%TZ > "$0.ok"; exit 0; }; sleep 10; done; echo "FAILED $0" >> failed.txt'
# HG002 Q100 v1.1: outer 100 kb of each chromosome end (regions written by catalogue_03_assemble_windows.py hg002-regions)
cd "$ATLAS_REPO/work/windows"
if [ ! -s hg002.ok ]; then
  python "$ATLAS_REPO/pipeline/catalogue_03_assemble_windows.py" hg002-regions
  for t in 1 2 3; do samtools faidx -r hg002.regions "$ATLAS_S3_HTTPS/T2T/HG002/assemblies/hg002v1.1.fasta.gz" > hg002.fa 2> hg002.log && [ "$(grep -c ">" hg002.fa)" -eq "$(wc -l < hg002.regions)" ] && { date -u +%FT%TZ > hg002.ok; break; }; sleep 10; done
fi
echo "failed window jobs: $(cat fetch/failed.txt 2>/dev/null | wc -l); failed probe jobs: $(cat probes/failed.txt 2>/dev/null | wc -l)"
