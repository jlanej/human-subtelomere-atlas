"""Boundary-anchored 5 kb bins of an end window and their best CHM13 end (catalogue step 2a).

Used unchanged by pipeline/catalogue_04_align_bins.py (catalogue build) and by the query tool.
Queries per end (tip frame, s0 = junk_len + telo_len = telomere/subtelomere boundary):
  bins k = tip[s0 + 5000 k, s0 + 5000 (k + 1)) (the last bin kept if >= 2 kb);
  S5 / S10 / S20 = tip[s0, s0 + 5/10/20 kb) (kept if >= 1 kb); t20 = tip[s0, 20 kb) (earlier definition, comparison only).
Alignment: minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 against the 48 CHM13 v2.0 500 kb ends. Per query and CHM13 end:
  best alignment score (AS) among alignments covering >= 50% of the query (summarise()).
Backends: 'minimap2' (the minimap2 executable; exactly the catalogue computation) or 'mappy' (Python binding; same preset and
  -N 50, but mappy cannot set -p, so secondary alignments scoring < 0.8 x the primary are not reported, and AS / de are recomputed
  from the CIGAR / cs string with the asm20 scoring). The mappy backend is an approximation; see atlas/README.md.
"""
import os, re, io, shutil, subprocess, tempfile
import pandas as pd

BIN, MINW, EL = 5000, 2000, 500_000
FIXED = (('S5', 5000), ('S10', 10000), ('S20', 20000))
MM2_ARGS = ['-x', 'asm20', '-c', '--secondary=yes', '-N', '50', '-p', '0.2']
COLS = ['end_id', 'bin', 'qstart', 'qend', 'best_end', 'best_AS', 'best_ident', 'best_dist0', 'best_dist1', 'best_strand', 'cand', 'n_ends']
BIN_DTYPES = {'bin': str, 'cand': str, 'best_end': str}


def tile(s0, L):
    """Query spans of one end: list of (tag, start, end) in tip frame; tags '0', '1', ... (5 kb bins), 'S5', 'S10', 'S20', 't20'."""
    out = []
    for k in range((L - s0 + BIN - 1) // BIN):
        a, b = s0 + k * BIN, min(s0 + (k + 1) * BIN, L)
        if b - a >= MINW: out.append((str(k), a, b))
    for tag, w in FIXED:
        b = min(s0 + w, L)
        if b - s0 >= 1000: out.append((tag, s0, b))
    b = min(20000, L)
    if b - s0 >= 1000: out.append(('t20', s0, b))
    return out


def summarise(qn, al, span):
    """Best CHM13 end of one query. al = alignments (qlen, qs, qe, strand, target, ts, te, AS, ident)."""
    best = {}
    for a in al:
        if a[2] - a[1] < 0.5 * a[0]: continue
        if a[4] not in best or a[7] > best[a[4]][7]: best[a[4]] = a
    end_id, k = qn.rsplit('@', 1)
    rec = dict(end_id=end_id, bin=k, qstart=span[0], qend=span[1])
    if not best:
        rec.update(best_end='', best_AS=0, best_ident='', best_dist0='', best_dist1='', best_strand='', cand='', n_ends=0); return rec
    b = max(best.values(), key=lambda a: (a[7], a[4]))
    d0, d1 = (b[5], b[6]) if b[4].endswith('p') else (EL - b[6], EL - b[5])
    cand = sorted((a for a in best.values() if a[7] >= 0.5 * b[7]), key=lambda a: (-a[7], a[4]))
    rec.update(best_end=b[4], best_AS=b[7], best_ident=round(b[8], 4), best_dist0=d0, best_dist1=d1, best_strand=b[3],
               cand=';'.join(f'{a[4]}:{a[7]}:{a[8]:.4f}' for a in cand), n_ends=len(best))
    return rec


def format_row(r): return '\t'.join(str(r[c]) for c in COLS)


def paf_alignment(line):
    """(query name, alignment tuple) from one minimap2 PAF line with AS:i and de:f tags."""
    f = line.split('\t')
    AS = next(int(t[5:]) for t in f[12:] if t.startswith('AS:i:'))
    de = next(float(t[5:]) for t in f[12:] if t.startswith('de:f:'))
    return f[0], (int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), AS, 1 - de)


# ---------------------------------------------------------------- alignment backends
def find_minimap2(exe=None):
    exe = exe or os.environ.get('SUBTEL_ATLAS_MINIMAP2') or shutil.which('minimap2')
    return exe if exe and shutil.which(exe) else None


def align_minimap2(ref_fa, queries, threads=1, exe='minimap2'):
    """queries: list of (name, seq). Returns {name: [alignment tuples]} (PAF order)."""
    out = {}
    with tempfile.TemporaryDirectory(prefix='subtel_atlas_') as td:
        qfa = os.path.join(td, 'q.fa')
        with open(qfa, 'w') as f:
            for n, s in queries: f.write(f'>{n}\n{s}\n')
        p = subprocess.run([exe] + MM2_ARGS + ['-t', str(threads), ref_fa, qfa], capture_output=True, text=True)
        if p.returncode != 0: raise RuntimeError('minimap2 failed: ' + p.stderr[-2000:])
        for line in p.stdout.splitlines():
            q, a = paf_alignment(line); out.setdefault(q, []).append(a)
    return out


# asm20 scoring (minimap2 2.24-2.31): match 1, mismatch 4, gap open/extend 6/2 and 26/1 (dual affine), ambiguous base 1
_SC = dict(A=1, B=4, q=6, e=2, q2=26, e2=1, ambi=1)
_CS = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})')


def score_cs(cs):
    """(AS, de) recomputed from a minimap2 cs string: DP score with the asm20 dual-affine scoring and gap-compressed divergence."""
    s = 0; m = mm = ngo = 0
    for tok in _CS.findall(cs):
        op = tok[0]
        if op == ':': l = int(tok[1:]); s += _SC['A'] * l; m += l
        elif op == '*':
            if 'n' in tok[1:]: s -= _SC['ambi']
            else: s -= _SC['B']; mm += 1
        elif op in '+-':
            l = len(tok) - 1; s -= min(_SC['q'] + _SC['e'] * l, _SC['q2'] + _SC['e2'] * l); ngo += 1
    den = m + mm + ngo
    return s, (float(f'{(mm + ngo) / den:.4g}') if den else 0.0)


def align_mappy(ref_fa, queries, threads=1, aligner=None):
    """mappy approximation of align_minimap2 (see module docstring)."""
    import mappy
    from concurrent.futures import ThreadPoolExecutor
    a = aligner or mappy.Aligner(ref_fa, preset='asm20', best_n=50)
    def one(nq):
        n, s = nq; buf = mappy.ThreadBuffer(); res = []
        for h in a.map(s, buf=buf, cs=True):
            AS, de = score_cs(h.cs)
            res.append((len(s), h.q_st, h.q_en, '+' if h.strand == 1 else '-', h.ctg, h.r_st, h.r_en, AS, 1 - de))
        return n, res
    out = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for n, res in ex.map(one, queries):
            if res: out[n] = res
    return out


def bin_table(ends, ref_fa, backend='minimap2', threads=1, exe='minimap2'):
    """ends: list of (end_id, tip-frame window, s0). Returns the bins table exactly as catalogue_04 writes and catalogue_05/06 read
    it (a DataFrame parsed from the same TSV text, with the same dtypes)."""
    queries, spans = [], {}
    for e, s, s0 in ends:
        for tag, a, b in tile(int(s0), len(s)):
            qn = f'{e}@{tag}'; queries.append((qn, s[a:b])); spans[qn] = (a, b)
    al = align_minimap2(ref_fa, queries, threads, exe) if backend == 'minimap2' else align_mappy(ref_fa, queries, threads)
    lines = ['\t'.join(COLS)] + [format_row(summarise(qn, al.get(qn, []), spans[qn])) for qn, _ in queries]
    return pd.read_csv(io.StringIO('\n'.join(lines) + '\n'), sep='\t', dtype=BIN_DTYPES)
