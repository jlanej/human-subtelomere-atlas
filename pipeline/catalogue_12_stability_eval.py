"""Boundary-shift stability of the track-based mosaic (redesign acceptance test).
Input : work/catalogue/stability_subset.tsv (catalogue_12_stability_subset.py), work/catalogue/tracks_shift.paf.gz (subset windows
        re-aligned from shifted boundaries s0 + d, d in SHIFTS; catalogue_04t_align_tracks.run), families, and (allele level) the
        structural alleles of the current build (work/catalogue/sa_alleles.tsv).
For every end and shift: core structure (alleles.core), switches, S20 state (+ family) and structural allele (alleles.match_allele
against the catalogue alleles of the end's arm); compared with the same end at d = 0. Wilson 95% intervals.
Modes : tau  -> structure-level stability for a grid of family thresholds (work/catalogue/stability_tau_scan.tsv)
        final -> full table for the published families (work/catalogue/stability_results.tsv, stability_summary.tsv)
"""
import os, sys, types, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C, catalogue_06t_mosaic_tracks as M
from subtel_atlas import tracks as TK, alleles as AL
SHIFTS = (-1000, -500, -100, 0, 100, 500, 1000)

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(c - h, 4), round(c + h, 4)

def per_end(prefix, alleles=None):
    E = pd.read_csv(f'{prefix}_ends.tsv', sep='\t'); B = pd.read_csv(f'{prefix}_blocks.tsv', sep='\t')
    W = M.winfo(); rows = []
    for (e, d), g in B.groupby(['end_id', 'shift']):
        bl = [types.SimpleNamespace(is_home=r.is_home, source=r.source, sub_start=r.sub_start, sub_end=r.sub_end) for r in g.sort_values('block').itertuples()]
        _, st, sw = AL.core(bl)
        w = W[e]; sub_len = w['L'] - (w['s0'] + d); full = w['L'] >= C.WINDOW and sub_len >= AL.FULL_SUB - max(0, d)
        a = ''
        if alleles is not None:
            best, _ = AL.match_allele(st, sw, sub_len, full, alleles.get(w['grp'], [])); a = best['allele_id'] if best else ''
        rows.append(dict(end_id=e, shift=d, structure=st, switches=','.join(map(str, sw)), allele=a))
    R = pd.DataFrame(rows).merge(E[['end_id', 'shift', 'S20_state', 'S20_family', 'S20_best_end']], on=['end_id', 'shift'])
    return R

def compare(R, S):
    z = R[R['shift'] == 0].set_index('end_id'); out = []
    for d in SHIFTS:
        if d == 0: continue
        x = R[R['shift'] == d].set_index('end_id').reindex(S.end_id)
        zz = z.reindex(S.end_id)
        same_struct = (x.structure == zz.structure).values
        same_all = (x.allele == zz.allele).values
        s20 = ((x.S20_state == zz.S20_state) & ((x.S20_state != 'foreign') | (x.S20_family == zz.S20_family))).values
        out.append(pd.DataFrame(dict(end_id=S.end_id.values, shift=d, same_structure=same_struct, same_allele=same_all, same_S20=s20)))
    return pd.concat(out)

def summarise(Cmp, S, by):
    rows = []
    x = Cmp.merge(S, on='end_id')
    for (d, *k), g in x.groupby(['shift'] + by):
        r = dict(shift=d); r.update(dict(zip(by, k))); r['n'] = len(g)
        for m in ('same_structure', 'same_allele', 'same_S20'):
            kk = int(g[m].sum()); r[m] = round(kk / len(g), 4); lo, hi = wilson(kk, len(g)); r[m + '_ci'] = f'{lo}-{hi}'
        rows.append(r)
    return pd.DataFrame(rows)

def load_alleles():
    SA = pd.read_csv(f'{C.CW}/sa_alleles.tsv', sep='\t'); al = {}
    for r in SA.itertuples():
        sw = [int(v) for v in str(r.switches_median).split(',')] if isinstance(r.switches_median, str) and r.switches_median else []
        al.setdefault(r.arm, []).append(dict(allele_id=r.allele_id, structure=r.structure, sw=sw, n=int(r.n_full), native=r.is_native))
    return al

if __name__ == '__main__':
    mode = sys.argv[1]; S = pd.read_csv(f'{C.CW}/stability_subset.tsv', sep='\t'); paf = f'{C.CW}/tracks_shift.paf.gz'
    if mode == 'tau':
        single, pair = M.read_ties(f'{C.CW}/tie_counts.tsv'); res = []
        for tau in [float(t) for t in sys.argv[2].split(',')]:
            fam, rows, _ = TK.build_families(single, pair, tau=tau)
            fp = f'{C.CW}/fam_tau{tau}.tsv'; pd.DataFrame(rows if rows else [dict(family='', n_members=0, members='')]).to_csv(fp, sep='\t', index=False)
            pre = f'{C.CW}/stab_tau{tau}'; M.label(paf, pre, fam_path=fp, only=True)
            Cmp = compare(per_end(pre), S); sm = summarise(Cmp, S, ['native'])
            for r in sm.itertuples(): res.append(dict(tau=tau, n_families=len(rows), largest=max([r_['n_members'] for r_ in rows], default=1), shift=r.shift, native=r.native, n=r.n, same_structure=r.same_structure, same_S20=r.same_S20))
            print(tau, len(rows), sm[sm['shift'].abs() == 500][['shift', 'native', 'same_structure', 'same_S20']].to_dict('records'), flush=True)
        C.write_tsv_atomic(pd.DataFrame(res), f'{C.CW}/stability_tau_scan.tsv')
    else:
        pre = f'{C.CW}/stab_final'; M.label(paf, pre, only=True)
        R = per_end(pre, load_alleles()); Cmp = compare(R, S); R.to_csv(f'{C.CW}/stability_results.tsv', sep='\t', index=False)
        SA = pd.read_csv(f'{C.CW}/sa_alleles.tsv', sep='\t').set_index('allele_id')
        z = R[R['shift'] == 0].set_index('end_id').allele.reindex(S.end_id)
        S = S.assign(new_native=z.map(SA.is_native).fillna(False).values, new_tier=z.map(SA.tier if 'tier' in SA else pd.Series(dtype=str)).values)
        parts = [summarise(Cmp, S, []).assign(stratum='all'), summarise(Cmp, S, ['native']).assign(stratum='previous nativity'),
                 summarise(Cmp, S, ['new_native']).assign(stratum='new nativity'), summarise(Cmp, S, ['tier']).assign(stratum='previous tier'),
                 summarise(Cmp, S, ['arm_class']).assign(stratum='arm class')]
        out = pd.concat(parts); C.write_tsv_atomic(out, f'{C.CW}/stability_summary.tsv'); Cmp.to_csv(f'{C.CW}/stability_compare.tsv', sep='\t', index=False)
        print(out[out['shift'].abs() == 500].to_string(index=False))
