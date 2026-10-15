"""Shared paths and helpers for the telomere_cis track (hypothesis H1: telomere length is set in cis by the subtelomere allele).

Inputs are read from ATLAS_INPUT_DIR (analysis home holding tipdup/, results/, subtel/, hg002/, meta/; set it in
work/config.local.sh, see pipeline/config.sh). Outputs are repo-relative: analyses/telomere_cis/{results,figures,data};
large intermediates under work/telomere_cis/ (git-ignored).
"""
import os, sys, hashlib
import numpy as np, pandas as pd

TRACK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(os.path.dirname(TRACK))
RES = os.path.join(TRACK, 'results'); FIG = os.path.join(TRACK, 'figures'); DAT = os.path.join(TRACK, 'data')
WORK = os.path.join(REPO, 'work', 'telomere_cis')
for d in (RES, FIG, DAT, WORK): os.makedirs(d, exist_ok=True)
sys.path.insert(0, os.path.join(REPO, 'pipeline'))
from catalogue_telomere import find_block, CAN   # same telomere caller/parameters as the Release 2 tip table
SEED = 20261015

def inp(*p):
    d = os.environ.get('ATLAS_INPUT_DIR')
    if not d: raise SystemExit('ATLAS_INPUT_DIR is not set (source pipeline/config.sh)')
    return os.path.join(d, *p)

_comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(_comp)[::-1]

def canon_frac(s):
    """fraction of bases covered by a canonical CCCTAA-rotation 6-mer (tip frame)"""
    if len(s) < 6: return 0.0
    m = np.zeros(len(s), bool)
    for i in range(len(s) - 5):
        if s[i:i + 6] in CAN: m[i:i + 6] = True
    return float(m.mean())

def arm_order(a):
    """sort key for arm names like '1p', '12q', 'Xp'"""
    c, s = a[:-1], a[-1]
    c = {'X': 23, 'Y': 24}.get(c, None) if not c.isdigit() else int(c)
    return (c if c is not None else 99, s)

ACRO_P = {'13p', '14p', '15p', '21p', '22p'}

def sha256(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def load_tips():
    """Named Release 2 tips with telomere length and allele state. One row per telomeric contig tip assigned to an arm."""
    F = pd.read_csv(inp('results', 'subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    keep = ['q', 'asm', 'sample', 'hap', 'assembler', 'contig', 'end', 'own', 'hifi_input', 'contig_len', 'junk_len',
            'telo_len', 't20_state', 'foreign_arm', 'form', 'junction_cluster_samples', 'vk_junction_concordant',
            'lay_reads_span_1kb']
    F = F[keep].rename(columns={'own': 'arm'})
    F['hapn'] = F.hap.map({'hap1': 1, 'mat': 1, 'hap2': 2, 'pat': 2})
    return F
