"""Same-haplotype check of the track-based mosaic against Tan et al. 2026 records (cached by the cross-check in work/tan2026:
tan_tips.fa, tan_windows.tsv, pairs.tsv; no Tan sequence is redistributed). Primary pairs as in data/validation (verified: identity
>= 0.99 and >= 50% of the subtelomere aligned; Tan record has a telomere). The Tan windows are aligned and labelled with the catalogue
method (tracks.py, published paralog families), Tan core structures are clipped at the shorter subtelomere of the pair, and assigned
with alleles.match_allele against the catalogue alleles of the arm. Agreement = Tan allele equals the catalogue end's allele, over
pairs assigned on both sides (as in the cross-check), and with unassigned Tan records counted as disagreement.
Outputs: work/catalogue/tan_pairs_tracks.tsv, tan_pairs_tracks.json"""
import os, sys, json, types, multiprocessing as mp, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C, catalogue_04t_align_tracks as A4, catalogue_06t_mosaic_tracks as M, catalogue_12_stability_eval as SE
from subtel_atlas import tracks as TK, alleles as AL
TW = f'{C.WORK}/tan2026'; _G = {}

def lab(item):
    q, al = item; e = q.rsplit('@', 1)[0]; w = _G['W'][e]
    bl, rec = TK.mosaic_end(e, al, w['s0'], w['L'], w['grp'], _G['fam'], M.NAMES)
    return e, bl, rec

def core_clip(bl, L):
    b = [x for x in bl if x['sub_start'] < L]
    if not b: return '', []
    b[-1] = dict(b[-1], sub_end=min(b[-1]['sub_end'], L))
    blk = [types.SimpleNamespace(is_home=x['label'] == 'home', source=x['label'] if x['label'] != 'home' else x['source'], sub_start=x['sub_start'], sub_end=x['sub_end']) for x in b]
    _, st, sw = AL.core(blk); return st, sw

def truncate(st, sw, L):
    keep = [s for s in sw if s < L - AL.MIN_BLOCK]; return '|'.join(st.split('|')[:len(keep) + 1]) if st else '', keep

def main():
    P = pd.read_csv(f'{TW}/pairs.tsv', sep='\t', low_memory=False); TWn = pd.read_csv(f'{TW}/tan_windows.tsv', sep='\t')
    P = P.merge(TWn[['end_id', 'has_telomere', 'tan_arm']].rename(columns={'end_id': 'tan_end_id'}), on='tan_end_id', how='left')
    prim = P[(P.ident >= 0.99) & (P.sub_aligned_frac >= 0.5) & P.has_telomere.astype(bool)].copy()
    paf = f'{TW}/tan_tracks.paf.gz'
    if not os.path.exists(paf): A4.run(f'{TW}/tan_tips.fa', f'{TW}/tan_windows.tsv', paf, C.THREADS, only=set(prim.tan_end_id))
    _G['W'] = {r.end_id: dict(s0=int(r.junk_len + r.telo_len), L=int(r.window_len), grp=C.group_of(r.tan_arm)) for r in TWn.itertuples()}
    _G['fam'] = TK.load_families(f'{C.CAT}/paralog_families.tsv')
    res = {}
    with mp.get_context('fork').Pool(C.THREADS) as pool:
        for e, bl, rec in pool.imap(lab, M.groups(paf), chunksize=8): res[e] = (bl, rec)
    al = SE.load_alleles()
    SM = pd.read_csv(f'{C.CW}/sa_members.tsv', sep='\t').set_index('end_id'); ME = pd.read_csv(f'{C.CW}/mosaic_ends.tsv', sep='\t').set_index('end_id')
    rows = []
    for r in prim.itertuples():
        c = r.end_id; t = r.tan_end_id; arm = SM.at[c, 'arm'] if c in SM.index else ''
        if t not in res or c not in SM.index: rows.append(dict(end_id=c, tan_end_id=t, status='missing')); continue
        bl, rec = res[t]; Lt = _G['W'][t]['L'] - _G['W'][t]['s0']; Lc = int(SM.at[c, 'sub_len']); L = min(Lt, Lc)
        st, sw = core_clip(bl, L); full = L >= AL.FULL_SUB
        best, cand = AL.match_allele(st, sw, L, full, al.get(arm, []))
        cs = SM.at[c, 'structure'] if isinstance(SM.at[c, 'structure'], str) else ''
        csw = [int(v) for v in str(SM.at[c, 'switches']).split(',')] if isinstance(SM.at[c, 'switches'], str) and SM.at[c, 'switches'] else []
        a1, _ = truncate(cs, csw, L); b1, _ = truncate(st, sw, L)
        rows.append(dict(end_id=c, tan_end_id=t, arm=arm, boundary_delta=r.boundary_delta, sub_len_common=L, full=full,
                         allele_c=SM.at[c, 'allele_id'], assignment_c=SM.at[c, 'assignment'], allele_t=best['allele_id'] if best else '',
                         structure_c=cs, structure_t=st, structure_agree=a1 == b1,
                         S20_c=f"{ME.at[c, 'S20_state']}:{ME.at[c, 'S20_family'] if isinstance(ME.at[c, 'S20_family'], str) else ''}",
                         S20_t=f"{rec['S20_state']}:{rec['S20_family'] or ''}"))
    X = pd.DataFrame(rows); X.to_csv(f'{C.CW}/tan_pairs_tracks.tsv', sep='\t', index=False)
    both = X[(X.assignment_c == 'full') & X.full & (X.allele_t != '')]
    strict = X[(X.assignment_c == 'full') & X.full]
    agree = lambda d: float((d.allele_c == d.allele_t).mean()) if len(d) else float('nan')
    bd = both.boundary_delta.abs()
    cats = {'exact': bd == 0, '1-100 bp': (bd > 0) & (bd <= 100), '101 bp-1 kb': (bd > 100) & (bd <= 1000), '>1 kb': bd > 1000}
    S0 = json.load(open(f'{C.REPO}/data/validation/tan2026_summary.json'))['same_haplotype']
    before_bd = {k: v.get('sa_agree') for k, v in S0.get('by_boundary_delta', {}).items()}
    J = dict(n_primary=int(len(X)), n_compared=int(len(both)), sa_agree=round(agree(both), 4), sa_agree_strict=round(agree(strict), 4), n_strict=int(len(strict)),
             sa_agree_before=S0['primary']['sa_agree'], structure_agree=round(float(X.structure_agree.mean()), 4),
             structure_agree_before=S0['primary']['core_structure_agree'], S20_agree=round(float((X.S20_c == X.S20_t).mean()), 4),
             S20_agree_before=S0['primary']['S20_state_and_end_agree'],
             sa_agree_by_boundary={k: round(agree(both[m]), 4) for k, m in cats.items()}, n_by_boundary={k: int(m.sum()) for k, m in cats.items()},
             sa_agree_by_boundary_before=before_bd)
    json.dump(J, open(f'{C.CW}/tan_pairs_tracks.json', 'w'), indent=1); print(json.dumps(J, indent=1))

if __name__ == '__main__':
    main()
