"""Shared helpers for the junctions track (H2: mechanism that created foreign terminal blocks).

Input locations come from the environment (pipeline/config.sh; machine-specific values in work/config.local.sh):
  ATLAS_INPUT_DIR  analysis home holding tipdup/, subtel/, results/, hg002/ (prior subtelomere placement results)
  ATLAS_THREADS    threads for minimap2 (default 3)
Outputs are repo-relative: analyses/junctions/{results,figures}/ (small, released) and work/junctions/ (large, git-ignored).

Coordinate conventions (inherited from the placement step that produced results/subtel_foreign_tips.tsv):
  tip frame      : tip sequence oriented telomere first (R contig ends reverse-complemented); s0 = masked telomere+outer junk
  subtel coord   : tip-frame coordinate minus s0 (x = switch_final - s0)
  CHM13 ends     : subtel/foreign/ref/chm13_ends500k.fa, outer 500 kb of each CHM13 v2.0 end in GENOME orientation
                   (p ends: telomere at 0; q ends: telomere at 500,000). Distance from telomere: p = pos, q = 500,000 - pos.
  junction       : sw_own_dist = distance from the own (recipient) CHM13 telomere at the switch; sw_for_dist = same on the
                   donor (foreign) CHM13 end.
Terminology: carrier = haplotype tip carrying a foreign terminal block; recipient (R) = the carrier's own arm; donor (D) = the
CHM13 end that the terminal block matches best (form 'R>D').
"""
import os, re, numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TRACK = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(TRACK))
RES = os.path.join(TRACK, 'results'); FIG = os.path.join(TRACK, 'figures')
WORK = os.path.join(REPO, 'work', 'junctions')
for d in (RES, FIG, WORK): os.makedirs(d, exist_ok=True)
THREADS = int(os.environ.get('ATLAS_THREADS', '3'))
EL = 500_000
SEED = 20260715

def inp(*p):
    base = os.environ.get('ATLAS_INPUT_DIR')
    if not base: raise SystemExit('ATLAS_INPUT_DIR is not set (source pipeline/config.sh)')
    return os.path.join(base, *p)

_comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(_comp)[::-1]

def readfa(fn):
    name = None; buf = []
    for line in open(fn):
        if line[0] == '>':
            if name is not None: yield name, ''.join(buf)
            name = line[1:].strip().split()[0]; buf = []
        else: buf.append(line.strip())
    if name is not None: yield name, ''.join(buf)

_CHM = None
def chm13_ends():
    """dict end name ('chr3q') -> genome-orientation 500 kb sequence (upper case)"""
    global _CHM
    if _CHM is None: _CHM = {n: s.upper() for n, s in readfa(inp('subtel/foreign/ref/chm13_ends500k.fa'))}
    return _CHM
def chm13_teloframe(end):
    """CHM13 end oriented telomere first (index = distance from the telomere)"""
    s = chm13_ends()[end]; return s if end.endswith('p') else rc(s)
def dist_from_telo(end, pos):
    return int(pos) if end.endswith('p') else EL - int(pos)

def load_tip_seqs(qs):
    """tip-frame sequences (telomere masked with N) for named Release 2 tips from subtel/foreign/r2/tips_input.tsv"""
    qs = set(qs); out = {}
    with open(inp('subtel/foreign/r2/tips_input.tsv')) as f:
        h = f.readline().rstrip('\n').split('\t'); iq, isq, is0, isrc = h.index('q'), h.index('seq'), h.index('s0'), h.index('src')
        for line in f:
            p = line.rstrip('\n').split('\t')
            if p[iq] in qs: out[p[iq]] = dict(seq=p[isq], s0=int(p[is0]), src=p[isrc])
    return out

# ---- per-base alignment tracks (same definitions as the placement step, mappy asm20, hits >= 2 kb) ----
_AL = {}
def aligner(end):
    import mappy
    if end not in _AL: _AL[end] = mappy.Aligner(seq=chm13_ends()[end], preset='asm20')
    return _AL[end]
def hits(end, seq, minlen=2000):
    return [(h.q_st, h.q_en, h.strand, h.r_st, h.r_en, h.cs) for h in aligner(end).map(seq, cs=True) if h.q_en - h.q_st >= minlen]
def err_track_cs(h, n):
    """per-query-base error indicator (1 = mismatch/indel/unaligned) and query->target position map for one mappy hit"""
    qs, qe, strand, ts, te, cs = h
    e = np.ones(n, np.int8); tpos = np.full(n, -1, np.int64); t = ts
    q = qs if strand == 1 else qe - 1; step = 1 if strand == 1 else -1
    for op, val in re.findall(r'([:*+\-~])([0-9]+|[a-z]+)', cs):
        if op == ':':
            l = int(val); idx = np.arange(q, q + step * l, step); e[idx] = 0; tpos[idx] = np.arange(t, t + l); q += step * l; t += l
        elif op == '*':
            if 0 <= q < n: e[q] = 1; tpos[q] = t
            q += step; t += 1
        elif op == '+':
            l = len(val); idx = np.arange(q, q + step * l, step); e[idx] = 1; tpos[idx] = t; q += step * l
        elif op == '-':
            if 0 <= q < n: e[q] = 1
            t += len(val)
    return e, tpos
def tracks(end, seq):
    """min-over-hits error track, 'aligned' indicator and per-hit (e, tpos) list of seq against one CHM13 end"""
    H = hits(end, seq); n = len(seq)
    if not H: return np.ones(n, np.int8), np.zeros(n, bool), []
    E = [err_track_cs(h, n) for h in H]
    e = np.min(np.vstack([x[0] for x in E]), axis=0); al = np.max(np.vstack([x[1] >= 0 for x in E]), axis=0)
    return e, al, list(zip(H, E))
def changepoint(ef, eo):
    cf = np.concatenate([[0], np.cumsum(ef)]); co = np.concatenate([[0], np.cumsum(eo)])
    cost = cf + (co[-1] - co); x = int(np.argmin(cost)); return x, int(cost[x])

# ---- telomere-like repeats ----
HEX = ('TTAGGG', 'CCCTAA')
def _deg_set():
    s = set()
    for m in HEX:
        for i in range(6):
            for b in 'ACGT': s.add(m[:i] + b + m[i + 1:])
    return s
DEG = _deg_set()
def its_tracks(seq):
    """per-position indicators: exact TTAGGG/CCCTAA hexamer start; degenerate (<= 1 mismatch) hexamer start; start of an
    in-phase tandem run of >= 3 degenerate hexamers (>= 18 bp ITS-like array) -- same-strand motif family within a run"""
    s = seq.upper(); n = len(s)
    ex = np.zeros(n, np.int8); dg = np.zeros(n, np.int8); fam = np.zeros(n, np.int8)
    dg_f = set(); dg_r = set()
    for i in range(6):
        for b in 'ACGT': dg_f.add('TTAGGG'[:i] + b + 'TTAGGG'[i + 1:]); dg_r.add('CCCTAA'[:i] + b + 'CCCTAA'[i + 1:])
    for i in range(n - 5):
        k = s[i:i + 6]
        if k in ('TTAGGG', 'CCCTAA'): ex[i] = 1
        if k in dg_f: dg[i] = 1; fam[i] = 1
        elif k in dg_r: dg[i] = 1; fam[i] = 2
    run = np.zeros(n, np.int8)
    # run of >= 3 consecutive in-phase degenerate hexamers of one family, any rotation not required (phase fixed by start)
    for i in range(n - 17):
        f = fam[i]
        if f and fam[i + 6] == f and fam[i + 12] == f: run[i] = 1
    return ex, dg, run
def wsum(cs, a, b):
    a = max(a, 0); b = min(b, len(cs) - 1); return int(cs[b] - cs[a]) if b > a else 0

def md5(fn):
    import hashlib
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
def sha256(fn):
    import hashlib
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
