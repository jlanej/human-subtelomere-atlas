"""Regenerate tests/fixtures/expected_windows.<backend>.tsv after a catalogue rebuild (the golden outputs embed catalogue values such
as allele counts, frequencies and representatives). Uses the same reference set-up as tests/test_query_tool.py. Prints which output
columns changed for each fixture window. Usage: python tests/regenerate_golden.py"""
import os, sys, io, shutil, tempfile
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIX = os.path.join(REPO, 'tests', 'fixtures')
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import query as Q, bins as B, refs as R
d = tempfile.mkdtemp(dir=os.path.join(REPO, 'work')) if os.path.isdir(os.path.join(REPO, 'work')) else tempfile.mkdtemp()
shutil.copy(os.path.join(FIX, 'chm13_ends_subset.fa.gz'), os.path.join(d, 'chm13_ends500k.fa.gz'))
for f in R.CAT_FILES: shutil.copy(os.path.join(REPO, 'data', 'catalogue', f), os.path.join(d, f))
ref = Q.Reference(d, verify=False); windows = list(Q.iter_records(os.path.join(FIX, 'windows.fa.gz')))
for backend in ['mappy'] + (['minimap2'] if B.find_minimap2() else []):
    recs = Q.analyse([Q.window_end(n, c, s) for n, c, s in windows], ref, backend, threads=2)
    buf = io.StringIO(); Q.write_tsv(recs, buf, Q.meta(backend, {}, len(windows), len(recs)))
    fn = os.path.join(FIX, f'expected_windows.{backend}.tsv'); old = open(fn).read()
    def rows(t): L = [l for l in t.splitlines() if not l.startswith('#')]; h = L[0].split('\t'); return h, [dict(zip(h, l.split('\t'))) for l in L[1:]]
    h0, r0 = rows(old); h1, r1 = rows(buf.getvalue())
    for a, b in zip(r0, r1):
        diff = [c for c in h1 if a.get(c) != b.get(c)]
        print(backend, b.get('query_id'), 'changed:', ', '.join(f"{c}: {a.get(c)} -> {b.get(c)}" for c in diff) if diff else 'none')
    open(fn + '.tmp', 'w').write(buf.getvalue()); os.replace(fn + '.tmp', fn)
