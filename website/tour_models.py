"""Formal claim models (criteria v1) for the guided read tour, in the tip frame of the tour's CHM13 end contigs.

For every definitive claim: the catalogue model end (claims.tsv model_end) -> subtel_atlas.readevidence.ClaimContext, the same
code that produced the verdicts (evidence/scripts/replicate_claim.py). Writes docs/evidence/tour/data/claims/<slug>/model.json
and model.bed (features on the tour reference's D and R contigs only): paralog-diagnostic SNV sites (tip-frame position on the donor end D and on the proximal comparator end P, with the
D and P letters), homologous blocks and D-/P-private sequence (>= 200 bp runs). Positions are 0-based tip frame (0 = telomeric
tip), the coordinates of the CHM13_<end>_end contigs of the tour.

Needs the CHM13 v2.0 analysis set (bgzip, ~1 GB; downloaded into the cache) and the CHM13 end references
(subtel_atlas.refs.fetch_refs). Usage: python website/tour_models.py   (cache: $TOUR_CACHE or work/tour/)
"""
import os, sys, json, subprocess, urllib.request, shutil
import numpy as np, pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.join(REPO, 'evidence', 'scripts'))
from subtel_atlas import readevidence as RE, refs  # noqa: E402
import evidence_lib as L  # noqa: E402
CACHE = os.environ.get('TOUR_CACHE', os.path.join(REPO, 'work', 'tour'))
OUT = os.path.join(REPO, 'docs', 'evidence', 'tour', 'data', 'claims')


def download(url, dst):
    if os.path.exists(dst): return dst
    os.makedirs(os.path.dirname(dst), exist_ok=True); tmp = dst + '.part'
    with urllib.request.urlopen(url, timeout=900) as r, open(tmp, 'wb') as fh: shutil.copyfileobj(r, fh, 1 << 22)
    os.replace(tmp, dst); return dst


def runs(mask, offset):
    m = np.asarray(mask, bool)
    if not m.any(): return []
    d = np.diff(np.r_[0, m.astype(int), 0]); st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    return [(int(a) + offset, int(b) + offset) for a, b in zip(st, en) if b - a >= 200]


def main():
    chm13 = os.path.join(CACHE, 'chm13', 'chm13v2.0.fa.gz')
    for ext in ('.gzi', ''): download(L.CHM13_FA_URL + ext, chm13 + ext)
    if not os.path.exists(chm13 + '.fai'): subprocess.run(['samtools', 'faidx', chm13], check=True)
    rdir = os.path.join(CACHE, 'refs')
    if not (os.path.exists(os.path.join(rdir, refs.CHM13_FA)) and all(refs.verify(rdir, quiet=True).values())):
        refs.fetch_refs(rdir, chm13=chm13, catalogue=os.path.join(REPO, 'data', 'catalogue'), log=print)
    ends = RE.Chm13Ends(os.path.join(rdir, refs.CHM13_FA))
    crit = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    cat = os.path.join(REPO, 'data', 'catalogue')
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False)
    C = C[C.definitive_verdict.astype(str).str.startswith('definitive')]
    models = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'),
                                      end_ids=sorted(set(C.model_end)))
    M = {m['end_id']: m for m in models}
    for c in C.itertuples():
        model = M[c.model_end]
        ctx = RE.ClaimContext(model, ends, crit, preset='map-ont', mask_cache=os.path.join(CACHE, 'masks'))
        D, P = ctx.c['D'], ctx.c['P']
        sites = [dict(d_pos=int(dp), p_pos=int(pp), d_letter=da, p_letter=pa)
                 for (dp, pp), (da, pa) in zip(ctx.pm.sites.tolist(), ctx.pm.site_alleles)]
        blocks = [dict(d_start=int(b['d_start']), d_end=int(b['d_end']), p_start=int(b['p_start']), p_end=int(b['p_end']),
                       ident=float(b['ident'])) for b in ctx.pm.blocks]
        priv = dict(donor=runs(ctx.pm.d_private, ctx.d_win[0]), proximal=runs(ctx.pm.p_private, ctx.p_win[0]))
        out = dict(claim_id=c.claim_id, model_end=c.model_end, D=D, P=P, R=model['R'], pattern=model.get('pattern'),
                   d_win=[int(x) for x in ctx.d_win], p_win=[int(x) for x in ctx.p_win], sites=sites, blocks=blocks,
                   private=priv, criteria=RE.criteria_hash(crit))
        od = os.path.join(OUT, c.slug); os.makedirs(od, exist_ok=True)
        json.dump(out, open(os.path.join(od, 'model.json'), 'w'))
        cD, cP = f'CHM13_{D}_end', f'CHM13_{P}_end'
        rows = []
        show = {cD} | ({cP} if P == model['R'] else set())   # contigs of the tour reference (D and R ends)
        for s in sites:
            lab = f'{D[:-1]}:{s["d_letter"]} | {P[:-1]}:{s["p_letter"]}'
            rows.append((cD, s['d_pos'], lab, '213,94,0', 1)); rows.append((cP, s['p_pos'], lab, '0,114,178', 1))
        with open(os.path.join(od, 'model.bed'), 'w') as fh:
            fh.write('track name="diagnostic sites (claim model)" itemRgb="On"\n')
            for b in blocks:
                for cn, a, e in ((cD, b['d_start'], b['d_end']), (cP, b['p_start'], b['p_end'])):
                    if cn not in show: continue
                    fh.write(f'{cn}\t{a}\t{e}\tshared block ({b["ident"]:.1%} identical)\t0\t.\t{a}\t{a}\t170,170,170\n')
            for side, cn, col in (('donor', cD, '213,94,0'), ('proximal', cP, '0,114,178')):
                for a, e in priv[side]:
                    if cn not in show: continue
                    fh.write(f'{cn}\t{a}\t{e}\t{"donor" if side == "donor" else "proximal"}-only sequence\t0\t.\t{a}\t{a}\t{col}\n')
            for cn, pos, lab, col, _ in sorted(rows, key=lambda r: (r[0], r[1])):
                if cn not in show: continue
                fh.write(f'{cn}\t{pos}\t{pos + 1}\t{lab}\t0\t.\t{pos}\t{pos + 1}\t{col}\n')
        print(c.slug, D, P, model['R'], len(sites), 'sites', len(blocks), 'blocks', flush=True)


if __name__ == '__main__':
    main()
