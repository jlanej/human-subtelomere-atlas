"""Step 2c (track-based): paralog families and mosaics from one alignment pass per end (atlas/subtel_atlas/tracks.py).
Modes
  ties     : co-tie statistics (tracks.tie_counts) over the defining ends (full-length HPRC_R2 / CHM13 / HG002_Q100 ends with a reliable
             arm call) -> work/catalogue/tie_counts.tsv (single and pair counts).
  families : paralog families at co-tie threshold TAU (default tracks.FAM_TAU; env ATLAS_FAM_TAU) -> data/catalogue/paralog_families.tsv
             and work/catalogue/paralog_family_edges.tsv.
  label    : mosaics of all ends with the published families -> work/catalogue/mosaic_blocks.tsv, mosaic_ends.tsv (same layout as the
             earlier step 06, with label/family/best-end columns).
Input  : work/catalogue/tracks.paf.gz (catalogue_04t_align_tracks.py), ends_windows.tsv, arms.tsv.
Usage  : python catalogue_06t_mosaic_tracks.py ties|families|label [paf] [prefix]
"""
import os, sys, gzip, collections, multiprocessing as mp, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas import tracks as TK
NAMES = [l.split('\t')[0] for l in open(f'{C.CW}/chm13_ends500k.fa.fai')] if os.path.exists(f'{C.CW}/chm13_ends500k.fa.fai') else []
_G = {}

def groups(paf):
    """yield (query name, [alignment tuples]) from a PAF in query order"""
    cur, al = None, []
    with (gzip.open(paf, 'rt') if paf.endswith('.gz') else open(paf)) as fh:
        for line in fh:
            q, n, qs, qe, st, tn, ts, te, cs = TK.parse_paf(line)
            if q != cur:
                if cur is not None: yield cur, al
                cur, al = q, []
            al.append((qs, qe, st, tn, ts, te, cs))
    if cur is not None: yield cur, al

def _ties(item):
    q, al = item; e, d = q.rsplit('@', 1); w = _G['W'][e]
    if not w['defining'] or d != '0': return None
    tr = TK.tracks(w['L'] - w['s0'], al); g, I = TK.identity_grid(tr, w['s0'], w['L'], NAMES)
    return TK.tie_counts(I, NAMES)

def _label(item):
    q, al = item; e, d = q.rsplit('@', 1); w = _G['W'][e]; s0 = w['s0'] + int(d)
    bl, rec = TK.mosaic_end(e, al, s0, w['L'], w['grp'], _G['fam'], NAMES)
    rec['shift'] = int(d); rec['query'] = q
    for b in bl: b['shift'] = int(d)
    return bl, rec

def winfo():
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t'); A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
    d = {}
    for r in W.itertuples():
        conf = A.at[r.end_id, 'arm_confidence']
        d[r.end_id] = dict(s0=int(r.junk_len + r.telo_len), L=int(r.window_len), grp=A.at[r.end_id, 'catalogue_arm'],
                           defining=(r.set in ('HPRC_R2', 'CHM13', 'HG002_Q100') and conf in C.DEFINING_CONFIDENCE
                                     and r.window_len >= C.WINDOW and r.window_len - r.junk_len - r.telo_len >= 80000))
    return d

def ties(paf, out):
    _G['W'] = winfo(); single, pair = collections.Counter(), collections.Counter()
    with mp.get_context('fork').Pool(C.THREADS) as pool:
        for r in pool.imap(_ties, groups(paf), chunksize=8):
            if r: single.update(r[0]); pair.update(r[1])
    rows = [dict(a=a, b='', count=c) for a, c in sorted(single.items())] + [dict(a=a, b=b, count=c) for (a, b), c in sorted(pair.items())]
    C.write_tsv_atomic(pd.DataFrame(rows), out); print('tie statistics: ends', len(single), 'pairs', len(pair))

def read_ties(fn):
    T = pd.read_csv(fn, sep='\t', keep_default_na=False)
    single = collections.Counter({r.a: r.count for r in T[T.b == ''].itertuples()})
    pair = collections.Counter({(r.a, r.b): r.count for r in T[T.b != ''].itertuples()})
    return single, pair

def families(tau):
    single, pair = read_ties(f'{C.CW}/tie_counts.tsv')
    fam, rows, edges = TK.build_families(single, pair, tau=tau)
    F = pd.DataFrame(rows, columns=['family', 'n_members', 'members', 'n_edges', 'min_cotie_rate', 'median_cotie_rate', 'cotie_points', 'tied_points']); F['tau'] = tau; F['min_points'] = TK.FAM_MIN_POINTS; F['margin'] = TK.MARGIN; F['window_bp'] = TK.WIN
    C.write_tsv_atomic(F, f'{C.CAT}/paralog_families.tsv')
    C.write_tsv_atomic(pd.DataFrame(edges, columns=['a', 'b', 'cotie_points', 'cotie_rate']), f'{C.CW}/paralog_family_edges.tsv')
    print(F[['family', 'n_members', 'members', 'min_cotie_rate']].to_string(index=False))

def label(paf, prefix, fam_path=None, only=None):
    _G['W'] = winfo(); _G['fam'] = TK.load_families(fam_path or f'{C.CAT}/paralog_families.tsv')
    BR, ER, seen = [], [], set()
    def items():
        for q, al in groups(paf):
            seen.add(q); yield q, al
    with mp.get_context('fork').Pool(C.THREADS) as pool:
        for bl, rec in pool.imap(_label, items(), chunksize=8): BR += bl; ER.append(rec)
    if only is None:                                     # ends without any alignment
        for e in _G['W']:
            if f'{e}@0' not in seen: bl, rec = _label((f'{e}@0', [])); BR += bl; ER.append(rec)
    B = pd.DataFrame(BR); E = pd.DataFrame(ER)
    B = B.rename(columns={'source': 'best_end'}); B['source'] = B.label                 # source = label token (family/arm/home/none)
    B.loc[B.is_home, 'source'] = B.loc[B.is_home, 'best_end']                            # home blocks: best home end (as before)
    B = B.sort_values(['end_id', 'shift', 'block']); E = E.sort_values(['end_id', 'shift'])
    C.write_tsv_atomic(B, f'{prefix}_blocks.tsv'); C.write_tsv_atomic(E, f'{prefix}_ends.tsv'); print('ends', E.end_id.nunique(), 'records', len(E))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'ties': ties(f'{C.CW}/tracks.paf.gz', f'{C.CW}/tie_counts.tsv')
    elif mode == 'families': families(float(os.environ.get('ATLAS_FAM_TAU', TK.FAM_TAU)))
    elif mode == 'label':
        paf = sys.argv[2] if len(sys.argv) > 2 else f'{C.CW}/tracks.paf.gz'; prefix = sys.argv[3] if len(sys.argv) > 3 else f'{C.CW}/mosaic'
        label(paf, prefix)
