#!/usr/bin/env bash
# Verify an atlas release from a fresh checkout.
#   bash verify_release.sh              offline checks (checksums, tests, report/site inputs)
#   bash verify_release.sh --network    also re-fetch a sample of windows from public sources and check SHA-256, and re-derive one
#                                       read-evidence claim from public reads (evidence/replicate_claim.sh; needs samtools, minimap2,
#                                       mappy, pysam, ~1 GB download of CHM13 v2.0 on first use and ~8 GB RAM)
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
echo "5. Read evidence"
sumcheck evidence SHA256SUMS
sumcheck evidence/showcase SHA256SUMS
python3 - <<'PY' && ok "criteria hash a89c5eaf023ed032 (repository and package copies)" || bad "criteria hash"
import json, hashlib
h = lambda p: hashlib.sha256(json.dumps(json.load(open(p))['thresholds'], sort_keys=True).encode()).hexdigest()[:16]
a, b, c = h('evidence/criteria_v1.json'), h('atlas/subtel_atlas/data/criteria_v1.json'), h('evidence/criteria_v1_preregistered.json')
assert a == b == 'a89c5eaf023ed032', (a, b)
PY
python3 - <<'PY' && ok "evidence tables consistent (claims vs reads, definitive tier only removes v1-confirmed claims)" || bad "evidence tables"
import pandas as pd
C = pd.read_csv('evidence/claims.tsv', sep='\t'); R = pd.read_csv('evidence/reads.tsv.gz', sep='\t', low_memory=False)
d = C[C.definitive_verdict.str.startswith('definitive')]
assert (d.v1_verdict == 'confirmed').all()
n = R[(R.tier == 'definitive') & (R.counted == True)].groupby('claim_id').size()
assert all(n.get(c, 0) == k for c, k in zip(C.claim_id, C.n_definitive)), 'n_definitive'
assert (C.criteria_sha == 'a89c5eaf023ed032').all()
PY
if [ "$NET" = 1 ]; then
  echo "6. Read evidence re-derived from public reads (one claim)"
  mkdir -p work
  # scratch and reference cache stay inside the repository's git-ignored work/ unless the user set SUBTEL_REPLICATE_CACHE
  if SUBTEL_REPLICATE_CACHE="${SUBTEL_REPLICATE_CACHE:-$PWD/work/replicate_cache}" \
     bash evidence/replicate_claim.sh 'HG01258|20q|SA-20q-01' -t 4 > work/verify_release_replicate.log 2>&1; then
    ok "evidence/replicate_claim.sh HG01258|20q|SA-20q-01 matches"; else bad "replicate_claim (see work/verify_release_replicate.log)"; fi
  echo "7. Public re-fetch of a window sample"
  out=work/verify_release; mkdir -p "$out"
  if python3 pipeline/fetch_end_window.py --sample 30 --seed 1 --out "$out/windows.fa" --report "$out/report.tsv" >/dev/null 2>&1 \
     && python3 -c "import pandas as pd,sys; d=pd.read_csv(sys.argv[1],sep='\t'); print(f'    {(d.status==\"OK\").sum()}/{len(d)} OK'); sys.exit(0 if len(d) and (d.status=='OK').all() else 1)" "$out/report.tsv"; then
    ok "30 windows re-fetched, SHA-256 identical"; else bad "public window re-fetch (see $out)"; fi
fi
echo "passed $pass, failed $fail"
[ "$fail" = 0 ]
