#!/bin/bash
# Download the large public inputs of the catalogue pipeline into work/public/ (public mode: ATLAS_INPUT_DIR unset; see
# pipeline/INPUTS.md). Every download is checked against the size/SHA-256 recorded in data/inputs/public_files.tsv.
#   chm13    : T2T-CHM13 v2.0 analysis set (bgzipped, 0.98 GB) -> work/public/ref/chm13v2.0.fa.gz (+ .fai/.gzi by samtools faidx)
#   pedigree : outer 250 kb of every chromosome end of the 8 WashU pedigree v1.1 haplotypes (range reads, 92 MB)
#              -> work/public/subtel/pedigree/ends250k.fa, ends250k_PAN028_h2.fa (pipeline/fetch_public_pedigree.py)
#   verify   : range-read the released windows of a stratified sample of ends (N=${ATLAS_VERIFY_N:-300}) plus all 48 CHM13 ends
#              from their public sources and compare with window_sha256 in data/catalogue/ends.tsv (pipeline/fetch_end_window.py)
#   windows  : the same for every end of the release (20,652 windows; long) -> work/public/windows/released_windows.fa
# Release 2 and HG002 Q100 windows for a rebuild are range-read by steps 01-02 of run_catalogue.sh (no download needed here).
# Usage: bash pipeline/fetch_public_inputs.sh [chm13] [pedigree] [verify] [windows]    (default: chm13 pedigree)
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/config.sh"
cd "$ATLAS_REPO"; P=work/public; mkdir -p "$P/ref" "$P/subtel/pedigree" "$P/windows"
want() { awk -F'\t' -v p="$1" '$1==p{print $4}' data/inputs/public_files.tsv; }
sha() { python -c 'import hashlib,sys; h=hashlib.sha256(); f=open(sys.argv[1],"rb"); [h.update(b) for b in iter(lambda: f.read(1<<22), b"")]; print(h.hexdigest())' "$1"; }
check() { local got; got=$(sha "$P/$1"); if [ "$got" = "$(want "$1")" ]; then echo "OK  $1"; else echo "SHA-256 MISMATCH $1 ($got)"; exit 1; fi; }
T=${ATLAS_THREADS:-4}
for what in ${@:-chm13 pedigree}; do
  case "$what" in
    chm13)
      f=$P/ref/chm13v2.0.fa.gz
      if [ ! -s "$f.ok" ]; then
        curl -fsSL --retry 5 -C - -o "$f.part" "$ATLAS_S3_HTTPS/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
        mv "$f.part" "$f"; check ref/chm13v2.0.fa.gz
        samtools faidx "$f"                                   # bgzipped upstream: builds .fai and .gzi
        date -u +%FT%TZ > "$f.ok"
      else echo "chm13 present ($f)"; fi ;;
    pedigree)
      if [ ! -s "$P/subtel/pedigree/ends.ok" ]; then
        python pipeline/fetch_public_pedigree.py
        check subtel/pedigree/ends250k.fa; check subtel/pedigree/ends250k_PAN028_h2.fa
        date -u +%FT%TZ > "$P/subtel/pedigree/ends.ok"
      else echo "pedigree ends present"; fi ;;
    verify)
      python pipeline/fetch_end_window.py --sample "${ATLAS_VERIFY_N:-300}" --threads "$T" --out "$P/windows/verify_sample.fa" --report "$P/windows/verify_sample.tsv"
      python - <<'EOF' > "$P/windows/chm13_ids.txt"
import pandas as pd; E = pd.read_csv('data/catalogue/ends.tsv', sep='\t', usecols=['end_id', 'set'], low_memory=False); print('\n'.join(E[E.set == 'CHM13'].end_id))
EOF
      python pipeline/fetch_end_window.py --ends-file "$P/windows/chm13_ids.txt" --threads "$T" --out "$P/windows/verify_chm13.fa" --report "$P/windows/verify_chm13.tsv" ;;
    windows)
      python pipeline/fetch_end_window.py --all --threads "$T" --out "$P/windows/released_windows.fa" --report "$P/windows/released_windows.tsv" ;;
    *) echo "unknown target $what"; exit 2 ;;
  esac
done
