#!/usr/bin/env bash
# Verify an atlas release from a fresh checkout.
#   bash verify_release.sh              offline checks (checksums, tests, report/site inputs)
#   bash verify_release.sh --network    also re-fetch a sample of windows from public sources and check SHA-256
# Requires python3 with pandas, numpy and pytest; --network also needs minimap2/samtools-free HTTP range access (S3).
set -euo pipefail
cd "$(dirname "$0")"
NET=0; [ "${1:-}" = "--network" ] && NET=1
pass=0; fail=0
ok()  { echo "  ok    $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL  $1"; fail=$((fail + 1)); }
sumcheck() {  # dir file
  if (cd "$1" && shasum -a 256 -c "$2" >/dev/null 2>&1); then ok "checksums $1/$2"; else bad "checksums $1/$2"; fi
}
echo "1. Release checksums"
sumcheck . data/catalogue/SHA256SUMS
sumcheck data/inputs SHA256SUMS
sumcheck literature SHA256SUMS
echo "2. Catalogue integrity"
if python3 tests/test_catalogue_integrity.py >/dev/null 2>&1; then ok "tests/test_catalogue_integrity.py"; else bad "tests/test_catalogue_integrity.py"; fi
echo "3. Query tool tests (offline)"
if PYTHONPATH=atlas python3 -m pytest -q -p no:cacheprovider tests/ >/dev/null 2>&1; then ok "pytest tests/"; else bad "pytest tests/"; fi
echo "4. Headline numbers"
python3 - <<'PY' && ok "headline counts match README" || bad "headline counts match README"
import pandas as pd, re
e = pd.read_csv('data/catalogue/ends.tsv', sep='\t', low_memory=False, usecols=['set', 'nativity_class'])
r2 = e[e['set'] == 'HPRC_R2']; sa = pd.read_csv('data/catalogue/alleles_structural.tsv', sep='\t'); hg = pd.read_csv('data/catalogue/haplogroups.tsv', sep='\t')
readme = open('README.md').read()
need = [f'{len(r2):,}', f'{len(sa)}', f'{len(hg)}', f'{(r2.nativity_class == "foreign").sum():,}', f'{(r2.nativity_class == "native").sum():,}']
missing = [n for n in need if n not in readme]
assert not missing, missing
PY
if [ "$NET" = 1 ]; then
  echo "5. Public re-fetch of a window sample"
  out=work/verify_release; mkdir -p "$out"
  if python3 pipeline/fetch_end_window.py --sample 30 --seed 1 --out "$out/windows.fa" --report "$out/report.tsv" >/dev/null 2>&1 \
     && python3 -c "import pandas as pd,sys; d=pd.read_csv(sys.argv[1],sep='\t'); print(f'    {(d.status==\"OK\").sum()}/{len(d)} OK'); sys.exit(0 if len(d) and (d.status=='OK').all() else 1)" "$out/report.tsv"; then
    ok "30 windows re-fetched, SHA-256 identical"; else bad "public window re-fetch (see $out)"; fi
fi
echo "passed $pass, failed $fail"
[ "$fail" = 0 ]
