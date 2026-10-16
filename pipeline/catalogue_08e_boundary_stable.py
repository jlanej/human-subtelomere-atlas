"""Step 2h (track-based): per-end boundary stability. For every end the identity tracks of window[s0:] (work/catalogue/tracks.paf.gz)
are shifted by d = -SHIFT and +SHIFT bp (tracks.shift_tracks: no re-alignment; d < 0 adds unmatched bases) and the mosaic, core
structure, S20 state (+ family) and structural allele (alleles.match_allele against the catalogue alleles) are recomputed.
boundary_stable_S20 / boundary_stable_allele: the call at both shifts equals the call at s0; boundary_stable = both.
The approximation is checked against real re-alignment of shifted windows on the stability subset (QC.md).
Output: work/catalogue/boundary_stable.tsv"""
import os, sys, types, multiprocessing as mp, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C, catalogue_06t_mosaic_tracks as M, catalogue_12_stability_eval as SE
from subtel_atlas import tracks as TK, alleles as AL
SHIFT = 500
_G = {}

def calls(e, tr, s0, L, grp, d):
    bl, rec = TK.mosaic_from_tracks(e, tr, s0 + d, L, grp, _G['fam'], M.NAMES)
    blk = [types.SimpleNamespace(is_home=b['label'] == 'home', source=b['label'] if b['label'] != 'home' else b['source'],
                                 sub_start=b['sub_start'], sub_end=b['sub_end']) for b in bl]
    _, st, sw = AL.core(blk) if blk else (None, '', [])
    sub_len = L - s0 - d; full = L >= C.WINDOW and sub_len >= AL.FULL_SUB - max(0, d)
    best, _ = AL.match_allele(st, sw, sub_len, full, _G['al'].get(grp, []))
    return st, best['allele_id'] if best else '', rec['S20_state'], rec['S20_family']

def one(item):
    q, al = item; e, d = q.rsplit('@', 1)
    if d != '0': return None
    w = _G['W'][e]; tr = TK.tracks(w['L'] - w['s0'], al); r = {}
    for dd in (-SHIFT, 0, SHIFT):
        r[dd] = calls(e, tr if dd == 0 else TK.shift_tracks(tr, dd), w['s0'], w['L'], w['grp'], dd)
    sa = all(r[dd][1] == r[0][1] for dd in (-SHIFT, SHIFT)); s20 = all(r[dd][2] == r[0][2] and r[dd][3] == r[0][3] for dd in (-SHIFT, SHIFT))
    return dict(end_id=e, boundary_stable_allele=sa, boundary_stable_S20=s20, boundary_stable=sa and s20,
                allele_minus500=r[-SHIFT][1], allele_plus500=r[SHIFT][1], S20_state_minus500=r[-SHIFT][2], S20_state_plus500=r[SHIFT][2])

def main(paf=f'{C.CW}/tracks.paf.gz', out=f'{C.CW}/boundary_stable.tsv'):
    _G['W'] = M.winfo(); _G['fam'] = TK.load_families(f'{C.CAT}/paralog_families.tsv'); _G['al'] = SE.load_alleles()
    rows = []
    with mp.get_context('fork').Pool(C.THREADS) as pool:
        for r in pool.imap(one, M.groups(paf), chunksize=8):
            if r: rows.append(r)
    seen = {r['end_id'] for r in rows}
    for e in _G['W']:                                   # ends without alignments: all-'none' mosaics, stable by construction
        if e not in seen: rows.append(dict(end_id=e, boundary_stable_allele=True, boundary_stable_S20=True, boundary_stable=True))
    R = pd.DataFrame(rows).sort_values('end_id'); C.write_tsv_atomic(R, out)
    print('ends', len(R), 'stable S20', round(R.boundary_stable_S20.mean(), 4), 'allele', round(R.boundary_stable_allele.mean(), 4), 'both', round(R.boundary_stable.mean(), 4))

if __name__ == '__main__':
    main()
