"""Step 3c: refresh the catalogue-file checksums shipped with the query-tool package (atlas/subtel_atlas/data/reference_checksums.json,
key catalogue_files) from data/catalogue/SHA256SUMS, so that `subtel-atlas fetch-refs/verify` accept exactly this release.
CHM13 end checksums and the catalogue version in that file are not touched."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
fn = os.path.join(C.REPO, 'atlas', 'subtel_atlas', 'data', 'reference_checksums.json')
X = json.load(open(fn)); sums = {l.split()[1]: l.split()[0] for l in open(f'{C.CAT}/SHA256SUMS')}
old = dict(X['catalogue_files'])
sys.path.insert(0, os.path.join(C.REPO, 'atlas')); from subtel_atlas import refs as R
X['catalogue_files'] = {f: sums[f'data/catalogue/{f}'] for f in R.CAT_FILES}
tmp = fn + '.tmp'; json.dump(X, open(tmp, 'w'), indent=1); open(tmp, 'a').write('\n'); os.replace(tmp, fn)
print('catalogue checksums updated:', [f for f in old if old[f] != X['catalogue_files'][f]])
