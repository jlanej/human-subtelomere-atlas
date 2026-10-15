#!/bin/bash
# Download the HG002 Q100 v1.1 diploid assembly (public T2T/HPRC bucket) into work/junctions/ for the whole-genome
# search of recipient terminal sequence (scripts/lost_sequence_hg002.py). ~1.77 GB; skipped when already present.
# Source: s3://human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz (Last-Modified 2024-06-23, ETag 06d48c014805443c6c3bf5f6110993b3-211)
set -euo pipefail
. "$(dirname "$0")/../../../pipeline/config.sh"
OUT="$ATLAS_REPO/work/junctions/hg002v1.1.fasta.gz"
URL="$ATLAS_S3_HTTPS/T2T/HG002/assemblies/hg002v1.1.fasta.gz"
if [ ! -s "$OUT" ] || [ "$(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT")" != "1769801343" ]; then
  curl -sS -L -C - -o "$OUT" "$URL"
fi
echo "hg002v1.1.fasta.gz $(stat -f%z "$OUT" 2>/dev/null || stat -c%s "$OUT") bytes"
