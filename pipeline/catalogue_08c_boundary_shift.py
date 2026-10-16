"""Step 2h: sensitivity of calls to the telomere-boundary position (bins and fixed windows are tiled from s0).
(a) All ends: the S20 window is re-aligned with the boundary moved by -SHIFT and +SHIFT bp (500 bp; clamped to the window) and the
    S20 state (and best CHM13 end if foreign) is recomputed with the catalogue rule. boundary_stable = both shifted S20 calls equal
    the call at s0. (Cheap: two 20 kb queries per end.)
(b) Seeded subset (N_SUB defining Release 2 ends, seed 20261016): every bin is re-tiled and re-aligned at s0 -/+ SHIFT, the mosaic,
    core structure and structural allele (shared matching rule against the catalogue alleles) are recomputed and compared with
    the same rule applied at s0. Allele-level stability for all ends would need a full re-alignment of every bin at each shift and
    is not computed.
Outputs: work/catalogue/boundary_shift_S20.tsv (all ends), boundary_shift_subset.tsv, boundary_shift_summary.tsv
"""
import os, sys, subprocess, types, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas.bins import MM2_ARGS, COLS, tile, summarise, paf_alignment
from subtel_atlas.mosaic import RefEnds, mosaic_end, t20_call, FIX
from subtel_atlas import alleles as AL
SHIFT, N_SUB, SEED = 500, 600, 20261016
REF = f'{C.CW}/chm13_ends500k.fa'

def align(queries, tag):
    qfa = f'{C.CW}/shift_{tag}.queries.fa'
    with open(qfa, 'w') as f:
        for qn, sq in queries: f.write(f'>{qn}\n{sq}\n')
    p = subprocess.run(['minimap2', *MM2_ARGS, '-t', str(C.THREADS), REF, qfa], capture_output=True, text=True, check=True)
    al = {}
    for line in p.stdout.splitlines(True):
        qn, a = paf_alignment(line); al.setdefault(qn, []).append(a)
    return al

def main():
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t').set_index('end_id')
    A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
    ME = pd.read_csv(f'{C.CW}/mosaic_ends.tsv', sep='\t').set_index('end_id')
    SM = pd.read_csv(f'{C.CW}/sa_members.tsv', sep='\t').set_index('end_id'); SA = pd.read_csv(f'{C.CW}/sa_alleles.tsv', sep='\t')
    seqs = {n: s for n, s in C.readfa(f'{C.WIN}/ends.fa')}
    s0 = (W.junk_len + W.telo_len).astype(int)
    # (a) S20 at shifted boundaries
    q, spans = [], {}
    for e, sq in seqs.items():
        for d in (-SHIFT, SHIFT):
            a = min(max(0, s0[e] + d), len(sq)); b = min(a + 20000, len(sq))
            if b - a >= 1000: qn = f'{e}@S20{d:+d}'; q.append((qn, sq[a:b])); spans[qn] = (a, b)
    al = align(q, 'S20')
    rows = []
    for e in seqs:
        homes = C.group_ends(A.at[e, 'catalogue_arm']); rec = dict(end_id=e)
        for d in (-SHIFT, SHIFT):
            qn = f'{e}@S20{d:+d}'
            if qn not in spans: rec[f'S20{d:+d}_state'] = 'none'; rec[f'S20{d:+d}_best_end'] = ''; continue
            r = summarise(qn, al.get(qn, []), spans[qn]); r = types.SimpleNamespace(**r) if r['best_end'] else None
            c = t20_call(r, homes); rec[f'S20{d:+d}_state'] = c['t20_state']; rec[f'S20{d:+d}_best_end'] = c['t20_best_end']
        rows.append(rec)
    S = pd.DataFrame(rows).set_index('end_id')
    ref_state = ME.S20_state.reindex(S.index); ref_end = ME.S20_best_end.reindex(S.index).fillna('')
    def same(d):
        st = S[f'S20{d:+d}_state']; be = S[f'S20{d:+d}_best_end'].fillna('')
        return (st == ref_state) & ((st != 'foreign') | (be == ref_end))
    S['boundary_stable'] = same(-SHIFT) & same(SHIFT)
    S.reset_index().to_csv(f'{C.CW}/boundary_shift_S20.tsv', sep='\t', index=False)
    # (b) subset: full re-tiling at shifted boundaries
    cand = sorted(SM.index[(SM.set == 'HPRC_R2') & SM.defining & SM.full & (SM.assignment == 'full')])
    sub = sorted(np.random.default_rng(SEED).choice(cand, min(N_SUB, len(cand)), replace=False).tolist())
    q, spans = [], {}
    for e in sub:
        sq = seqs[e]
        for d in (-SHIFT, SHIFT):
            for tg, a, b in tile(max(0, s0[e] + d), len(sq)):
                qn = f'{e}#{d:+d}@{tg}'; q.append((qn, sq[a:b])); spans[qn] = (a, b)
    al = align(q, 'subset')
    refs = RefEnds(dict(C.readfa(REF)))
    alle = {}
    for r in SA.itertuples():
        sw = [int(x) for x in str(r.switches_median).split(',')] if isinstance(r.switches_median, str) and r.switches_median else []
        alle.setdefault(r.arm, []).append(dict(allele_id=r.allele_id, structure=r.structure, sw=sw, n=int(r.n_full)))
    out = []
    for e in sub:
        grp = A.at[e, 'catalogue_arm']; sw0 = [int(x) for x in str(SM.at[e, 'switches']).split(',')] if isinstance(SM.at[e, 'switches'], str) else []
        base, _ = AL.match_allele(SM.at[e, 'structure'], sw0, int(SM.at[e, 'sub_len']), True, alle.get(grp, []))
        rec = dict(end_id=e, arm=grp, structure=SM.at[e, 'structure'], allele=SM.at[e, 'allele_id'], allele_rematch=base['allele_id'] if base else '',
                   S20=ME.at[e, 'S20_state'])
        for d in (-SHIFT, SHIFT):
            sh = s0[e] + d; rows_ = []
            for tg, a, b in tile(max(0, sh), len(seqs[e])):
                qn = f'{e}#{d:+d}@{tg}'; rr = summarise(qn, al.get(qn, []), spans[qn]); rr['end_id'] = e; rr['bin'] = tg; rows_.append(rr)
            B = pd.DataFrame(rows_)
            g = B[~B.bin.isin(FIX)].assign(k=lambda x: x.bin.astype(int)).sort_values('k')
            fixed = {t: types.SimpleNamespace(**B[B.bin == t].iloc[0].to_dict()) if (B.bin == t).any() and B[B.bin == t].iloc[0].best_end else None for t in FIX}
            g = g.assign(best_end=g.best_end.replace('', np.nan))
            BR, er = mosaic_end(refs, e, g, fixed, max(0, sh), grp, seqs[e])
            blk = [types.SimpleNamespace(**b) for b in BR]
            _, st, sw = AL.core(blk)
            best, _ = AL.match_allele(st, sw, len(seqs[e]) - max(0, sh), True, alle.get(grp, []))
            rec.update({f'structure{d:+d}': st, f'allele{d:+d}': best['allele_id'] if best else '', f'S20{d:+d}': er['S20_state']})
        out.append(rec)
    T = pd.DataFrame(out); T.to_csv(f'{C.CW}/boundary_shift_subset.tsv', sep='\t', index=False)
    summ = [dict(test='S20 stable at -/+500 bp (all ends)', n=len(S), stable=int(S.boundary_stable.sum())),
            dict(test='subset: S20 stable', n=len(T), stable=int(((T['S20-500'] == T.S20) & (T['S20+500'] == T.S20)).sum())),
            dict(test='subset: core structure stable', n=len(T), stable=int(((T['structure-500'] == T.structure) & (T['structure+500'] == T.structure)).sum())),
            dict(test='subset: structural allele stable (vs re-match at s0)', n=len(T), stable=int(((T['allele-500'] == T.allele_rematch) & (T['allele+500'] == T.allele_rematch)).sum())),
            dict(test='subset: re-match at s0 equals catalogue allele', n=len(T), stable=int((T.allele_rematch == T.allele).sum()))]
    Su = pd.DataFrame(summ); Su['fraction'] = (Su.stable / Su.n).round(4); Su.to_csv(f'{C.CW}/boundary_shift_summary.tsv', sep='\t', index=False)
    print(Su.to_string(index=False))

if __name__ == '__main__':
    main()
