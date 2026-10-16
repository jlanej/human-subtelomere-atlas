"""Snapshot the small derived inputs of the catalogue pipeline into data/inputs/ and inventory every input (pipeline/INPUTS.md).
Run once by the maintainer with the original analysis home (ATLAS_INPUT_DIR set); outside users never need to run it.
  * small tables are copied under the same relative path as in the analysis home, gzip-compressed deterministically (mtime 0);
    every file is checked for absolute/personal paths (any occurrence of the analysis-home path is rewritten to <analysis_home>,
    any other absolute user path aborts the snapshot);
  * index files needed for public range reads (HG002 Q100 .fai, pedigree .fai) are copied uncompressed;
  * large public files are not copied: their size and SHA-256 are recorded in data/inputs/public_files.tsv so that
    pipeline/fetch_public_inputs.sh can verify the downloads;
  * data/inputs/SHA256SUMS covers every snapshot file; work/inputs_inventory.tsv lists every analysis-home input (size, SHA-256).
Usage: ATLAS_INPUT_DIR=<analysis home> python pipeline/snapshot_inputs.py
"""
import os, sys, re, io, gzip, glob, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C

A = C.need_input(); OUT = os.path.join(C.REPO, 'data', 'inputs')
TABLES = ['tipdup/tips_meta.tsv', 'tipdup/tip_placement.tsv', 'meta/r2_assemblies_index.csv', 'meta/explorer_assemblies.csv',
          'results/subtel_foreign_tips.tsv', 'results/subtel_hg002_ends.tsv', 'results/subtel_pedigree_transmission.tsv']
PED_HAPS = ['PAN027_pat', 'PAN027_mat', 'PAN010_h1', 'PAN010_h2', 'PAN011_h1', 'PAN011_h2', 'PAN028_h1', 'PAN028_h2']
INDEXES = ['hg002/hg002v1.1.fasta.gz.fai'] + [f'subtel/pedigree/{h}.v1.1.fai' for h in PED_HAPS]
LARGE = [('ref/chm13v2.0.fa.gz', f'{C.S3}/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'),
         ('subtel/pedigree/ends250k.fa', 'range reads of https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/ (7 haplotypes)'),
         ('subtel/pedigree/ends250k_PAN028_h2.fa', 'range reads of https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/ (PAN028 hap2)')]
OPTIONAL = ['tipdup/tips_masked.fa']
BAD = re.compile(r'(/Users/|/home/[a-z]|C:\\\\Users)')

def sha(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''): h.update(b)
    return h.hexdigest()

def clean(text, name):
    text = text.replace(A.rstrip('/'), '<analysis_home>')
    m = BAD.search(text)
    if m: raise SystemExit(f'{name}: absolute user path in content near {text[max(0, m.start() - 40):m.end() + 40]!r}')
    return text

def write_gz(rel, text):
    fn = os.path.join(OUT, rel + '.gz'); os.makedirs(os.path.dirname(fn), exist_ok=True)
    buf = io.BytesIO()
    with gzip.GzipFile(filename='', mode='wb', fileobj=buf, mtime=0, compresslevel=9) as g: g.write(text.encode())
    open(fn + '.tmp', 'wb').write(buf.getvalue()); os.replace(fn + '.tmp', fn); return fn

inv, snap = [], []
for rel in TABLES:
    src = os.path.join(A, rel); text = clean(open(src).read(), rel)
    fn = write_gz(rel, text); snap.append(fn)
    inv.append(dict(path=rel, bytes=os.path.getsize(src), sha256=sha(src), snapshot=os.path.relpath(fn, C.REPO)))
for rel in INDEXES:
    src = os.path.join(A, rel); dst = os.path.join(OUT, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
    text = clean(open(src).read(), rel); open(dst + '.tmp', 'w').write(text); os.replace(dst + '.tmp', dst); snap.append(dst)
    inv.append(dict(path=rel, bytes=os.path.getsize(src), sha256=sha(src), snapshot=os.path.relpath(dst, C.REPO)))
pub = ['path\turl\tbytes\tsha256']
for rel, url in LARGE:
    src = os.path.join(A, rel); s_ = sha(src)
    inv.append(dict(path=rel, bytes=os.path.getsize(src), sha256=s_, snapshot='public download (work/public/' + rel + ')'))
    pub.append(f'{rel}\t{url}\t{os.path.getsize(src)}\t{s_}')
open(os.path.join(OUT, 'public_files.tsv'), 'w').write('\n'.join(pub) + '\n'); snap.append(os.path.join(OUT, 'public_files.tsv'))
for rel in OPTIONAL + ['ref/chm13v2.0.fa.gz.fai', 'ref/chm13v2.0.fa.gz.gzi']:
    src = os.path.join(A, rel)
    if os.path.exists(src): inv.append(dict(path=rel, bytes=os.path.getsize(src), sha256=sha(src), snapshot='not needed (optional or regenerated)'))
# optional cache of earlier 100 kb windows (per-assembly FASTA): aggregate checksum over sorted (name, sha256)
fas = sorted(glob.glob(os.path.join(A, 'subtel/foreign/r2_100k/*.fa')))
agg = hashlib.sha256(''.join(f'{os.path.basename(f)}\t{sha(f)}\n' for f in fas).encode()).hexdigest()
inv.append(dict(path=f'subtel/foreign/r2_100k/*.fa ({len(fas)} files)', bytes=sum(os.path.getsize(f) for f in fas), sha256='aggregate:' + agg,
                snapshot='not needed (re-fetched by steps 01-02)'))
with open(os.path.join(OUT, 'SHA256SUMS.tmp'), 'w') as f:
    for fn in sorted(snap): f.write(f'{sha(fn)}  {os.path.relpath(fn, OUT)}\n')
os.replace(os.path.join(OUT, 'SHA256SUMS.tmp'), os.path.join(OUT, 'SHA256SUMS'))
import pandas as pd
C.write_tsv_atomic(pd.DataFrame(inv), os.path.join(C.WORK, 'inputs_inventory.tsv'))
print('snapshot files', len(snap), 'total MB', round(sum(os.path.getsize(f) for f in snap) / 1e6, 2))
