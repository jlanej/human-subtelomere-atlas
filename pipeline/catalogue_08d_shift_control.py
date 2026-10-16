"""Control for step 08c: the subset code path (re-tiling, re-alignment, shared mosaic and core-structure functions) is run with a
zero boundary shift on the first N_CTRL subset ends; it must reproduce the catalogue profile and core structure exactly.
Output: work/catalogue/boundary_shift_zero_control.tsv"""
import os, sys, types, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_08c_boundary_shift as X, catalogue_common as C
from subtel_atlas.bins import tile, summarise
from subtel_atlas.mosaic import RefEnds, mosaic_end, FIX
from subtel_atlas import alleles as AL
N_CTRL = 120
T = pd.read_csv(f'{C.CW}/boundary_shift_subset.tsv', sep='\t'); sub = sorted(T.end_id)[:N_CTRL]
W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t').set_index('end_id'); A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
SM = pd.read_csv(f'{C.CW}/sa_members.tsv', sep='\t').set_index('end_id'); ME = pd.read_csv(f'{C.CW}/mosaic_ends.tsv', sep='\t').set_index('end_id')
seqs = {n: s for n, s in C.readfa(f'{C.WIN}/ends.fa') if n in set(sub)}
s0 = (W.junk_len + W.telo_len).astype(int)
q, spans = [], {}
for e in sub:
    for tg, a, b in tile(s0[e], len(seqs[e])): qn = f'{e}#0@{tg}'; q.append((qn, seqs[e][a:b])); spans[qn] = (a, b)
al = X.align(q, 'zero'); refs = RefEnds(dict(C.readfa(X.REF))); out = []
for e in sub:
    rows = []
    for tg, a, b in tile(s0[e], len(seqs[e])):
        qn = f'{e}#0@{tg}'; rr = summarise(qn, al.get(qn, []), spans[qn]); rr['end_id'] = e; rr['bin'] = tg; rows.append(rr)
    B = pd.DataFrame(rows); g = B[~B.bin.isin(FIX)].assign(k=lambda x: x.bin.astype(int)).sort_values('k')
    fixed = {t: types.SimpleNamespace(**B[B.bin == t].iloc[0].to_dict()) if (B.bin == t).any() and B[B.bin == t].iloc[0].best_end else None for t in FIX}
    g = g.assign(best_end=g.best_end.replace('', np.nan))
    BR, er = mosaic_end(refs, e, g, fixed, s0[e], A.at[e, 'catalogue_arm'], seqs[e])
    _, st, sw = AL.core([types.SimpleNamespace(**b) for b in BR])
    out.append(dict(end_id=e, profile_identical=er['profile'] == ME.at[e, 'profile'], structure_identical=st == SM.at[e, 'structure']))
O = pd.DataFrame(out); O.to_csv(f'{C.CW}/boundary_shift_zero_control.tsv', sep='\t', index=False)
print('zero-shift control: profile identical', int(O.profile_identical.sum()), '/', len(O), '; core structure identical', int(O.structure_identical.sum()), '/', len(O))
