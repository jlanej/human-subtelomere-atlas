"""Shared paths, constants and helpers for the allele-catalogue pipeline (pipeline/catalogue_*.py).
Input locations come from the environment (see pipeline/config.sh): ATLAS_INPUT_DIR = analysis home holding the
Release 2 tip tables (tipdup/), prior subtelomere results (results/, subtel/), CHM13 (ref/) and HG002 files (hg002/).
All outputs are repo-relative: large intermediates under work/, release tables under data/catalogue/.
"""
import os, hashlib, gzip

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT = os.environ.get('ATLAS_INPUT_DIR')
THREADS = int(os.environ.get('ATLAS_THREADS', '3'))
S3 = os.environ.get('ATLAS_S3_HTTPS', 'https://human-pangenomics.s3.amazonaws.com')
WORK = os.path.join(REPO, 'work')
WIN = os.path.join(WORK, 'windows')          # oriented 100 kb end windows (git-ignored)
CW = os.path.join(WORK, 'catalogue')         # catalogue intermediates (git-ignored)
CAT = os.path.join(REPO, 'data', 'catalogue')
FIGQC = os.path.join(REPO, 'figures', 'qc')
WINDOW = 100_000                             # telomere-anchored window length
HG002_URL = f'{S3}/T2T/HG002/assemblies/hg002v1.1.fasta.gz'
CHM13_URL = f'{S3}/T2T/CHM13/assemblies/chm13v2.0.fa'   # range-readable public copy (windows are read from a local copy)

def need_input():
    if not INPUT:
        raise SystemExit('ATLAS_INPUT_DIR is not set (source pipeline/config.sh)')
    return INPUT

def inp(*p): return os.path.join(need_input(), *p)

_comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(_comp)[::-1]

def readfa(fn):
    op = gzip.open if fn.endswith('.gz') else open
    name = None; buf = []
    with op(fn, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name is not None: yield name, ''.join(buf)
                name = line[1:].strip().split()[0]; buf = []
            else: buf.append(line.strip())
    if name is not None: yield name, ''.join(buf)

def sha256_str(s): return hashlib.sha256(s.encode()).hexdigest()
def sha256_file(fn, bs=1 << 20):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(bs), b''): h.update(b)
    return h.hexdigest()

def write_tsv_atomic(df, fn, **kw):
    tmp = fn + '.tmp'; df.to_csv(tmp, sep='\t', index=False, **kw); os.replace(tmp, fn)

def arm_of(chm13_end):
    """'chr20q' -> '20q'"""
    return chm13_end[3:] if isinstance(chm13_end, str) and chm13_end.startswith('chr') else chm13_end

# Catalogue arm groups: ends whose distal sequence is shared between chromosomes are catalogued as pooled groups.
#   acroP = 13p, 14p, 15p, 21p, 22p (acrocentric short arms); PAR1 = Xp, Yp; PAR2 = Xq, Yq. All other arms are their own group.
ARM_GROUPS = {'acroP': ['13p', '14p', '15p', '21p', '22p'], 'PAR1': ['Xp', 'Yp'], 'PAR2': ['Xq', 'Yq']}
_G = {a: g for g, arms in ARM_GROUPS.items() for a in arms}
def group_of(arm): return _G.get(arm, arm) if isinstance(arm, str) and arm else 'unassigned'
def group_ends(group):
    """CHM13 end names ('chr13p', ...) that count as 'home' for a catalogue arm group"""
    if group == 'unassigned': return []
    return ['chr' + a for a in ARM_GROUPS.get(group, [group])]
DEFINING_CONFIDENCE = {'named', 'high', 'medium', 'mosaic_group'}
