"""Step 2f: structural alleles per arm = recurring mosaic structures.
Core structure of an end (from work/catalogue/mosaic_blocks.tsv): blocks shorter than MIN_BLOCK (subtelomere coordinates) are
  dropped - their span is split at its midpoint between the neighbouring blocks - and adjacent blocks with the same source are
  merged. The structure string lists block sources telomere -> centromere ('home' = the end's own arm, '3q' = a block best
  matching CHM13 3q, 'none' = no CHM13 alignment). Switch vector = subtelomere coordinates of the core block boundaries.
Definition set: full-length ends (window_len = 100 kb and subtelomere >= FULL_SUB) from HPRC_R2, CHM13 and HG002_Q100 whose
  catalogue arm is reliable (arm_confidence in named/high/medium/mosaic_group; catalogue_common.DEFINING_CONFIDENCE).
  Within each (arm, structure) the switch vectors are clustered by complete linkage on the Chebyshev distance, cut at
  SW_TOL (every pair of members has every switch within SW_TOL). Allele ID SA-<arm>-NN (NN by number of ends, ties by structure).
  The native allele is the structure 'home' (no foreign block).
Assignment of other ends: partial ends (shorter windows/subtelomere) - the alleles whose structure, truncated at the end's observed
  subtelomere length, matches the end's structure with all observed switches within SW_TOL; unique -> 'partial_unique', several
  -> the nearest by median switch positions, then the most frequent ('partial_ambiguous', candidates listed); none -> 'unassigned_partial'. Full-length pedigree_aux ends and
  full-length ends with a low-confidence/conflicting arm call are matched against the defined alleles of their catalogue arm
  (same structure, every switch within SW_TOL): 'aux' / 'arm_lowconf', else 'unassigned_aux' / 'unassigned_arm_lowconf'.
Representative (medoid): the full member with the smallest summed Chebyshev distance to the other members (ties: end_id).
Shared code: atlas/subtel_atlas/alleles.py (core structure, clustering, matching; also used by the query tool).
Outputs: work/catalogue/sa_members.tsv, sa_alleles.tsv, core_blocks.tsv
"""
import os, sys, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas import alleles as AL     # structure rules shared with the query tool
MIN_BLOCK = int(os.environ.get('ATLAS_MIN_BLOCK', '10000')); SW_TOL = int(os.environ.get('ATLAS_SW_TOL', '10000'))
FULL_SUB = AL.FULL_SUB
DEF_SETS = AL.DEF_SETS

def core(bl): return AL.core(bl, MIN_BLOCK)

def prefix_compatible(struct, sw, a_struct, a_sw, Lobs): return AL.prefix_compatible(struct, sw, a_struct, a_sw, Lobs, MIN_BLOCK, SW_TOL)

def main():
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t').set_index('end_id')
    A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
    Bk = pd.read_csv(f'{C.CW}/mosaic_blocks.tsv', sep='\t')
    E = pd.read_csv(f'{C.CW}/mosaic_ends.tsv', sep='\t').set_index('end_id')
    rows, cb = [], []
    for e, g in Bk.groupby('end_id', sort=True):
        bl, st, sw = core(list(g.sort_values('block').itertuples()))
        sub_len = int(W.at[e, 'window_len'] - E.at[e, 's0'])
        full = W.at[e, 'window_len'] >= C.WINDOW and sub_len >= FULL_SUB
        rows.append(dict(end_id=e, arm=A.at[e, 'catalogue_arm'], set=W.at[e, 'set'], sample=W.at[e, 'sample'], structure=st,
                         switches=sw, n_core_blocks=len(bl), sub_len=sub_len, full=full,
                         defining=A.at[e, 'arm_confidence'] in C.DEFINING_CONFIDENCE))
        for i, b in enumerate(bl):
            cb.append(dict(end_id=e, core_block=i + 1, source=b['hsrc'] if b['source'] == 'home' else b['source'],
                           is_home=b['source'] == 'home', sub_start=int(round(b['s'])), sub_end=int(round(b['e']))))
    R = pd.DataFrame(rows); alle, mem, diag = [], [], []
    for arm, g in R.groupby('arm', sort=True):
        d = g[g.full & g.set.isin(DEF_SETS) & g.defining]; cl = []
        for st, h in d.groupby('structure', sort=True):
            for c in AL.cluster_structure(list(h.end_id), list(h.switches), SW_TOL):
                c['structure'] = st; cl.append(c)
        cl.sort(key=lambda c: (-len(c['members']), c['structure'], c['sw']))
        for r, c in enumerate(cl, 1):
            c['allele_id'] = f'SA-{arm}-{r:02d}'
            for m in c['members']: mem.append(dict(end_id=m, allele_id=c['allele_id'], assignment='full'))
        # partial and auxiliary ends
        for c in cl: c['n'] = len(c['members'])
        own = {m: c['allele_id'] for c in cl for m in c['members']}
        for x in d.itertuples():          # tie-break diagnostic: re-assign every defining end as if it were not defining
            bn, cn = AL.match_allele(x.structure, x.switches, x.sub_len, True, cl, MIN_BLOCK, SW_TOL)
            bl_, _ = AL.match_allele(x.structure, x.switches, x.sub_len, True, cl, MIN_BLOCK, SW_TOL, legacy=True)
            diag.append(dict(end_id=x.end_id, arm=arm, kind='defining_self', own=own[x.end_id], n_candidates=len(cn),
                             nearest=bn['allele_id'] if bn else '', legacy=bl_['allele_id'] if bl_ else ''))
        for x in g[~g.end_id.isin(d.end_id)].itertuples():
            kind = ('aux' if x.set == 'pedigree_aux' else 'arm_lowconf') if x.full else ('partial' if x.set != 'pedigree_aux' else 'aux_partial')
            best, cand = AL.match_allele(x.structure, x.switches, x.sub_len, x.full, cl, MIN_BLOCK, SW_TOL)
            if cand:
                bl_, _ = AL.match_allele(x.structure, x.switches, x.sub_len, x.full, cl, MIN_BLOCK, SW_TOL, legacy=True)
                diag.append(dict(end_id=x.end_id, arm=arm, kind=kind, own='', n_candidates=len(cand), nearest=best['allele_id'], legacy=bl_['allele_id']))
                asg = kind if kind in ('aux', 'arm_lowconf') else (f'{kind}_unique' if len(cand) == 1 else f'{kind}_ambiguous')
                mem.append(dict(end_id=x.end_id, allele_id=best['allele_id'], assignment=asg, candidates=','.join(c['allele_id'] for c in cand)))
            else:
                mem.append(dict(end_id=x.end_id, allele_id='', assignment='unassigned_aux' if kind.startswith('aux') else f'unassigned_{kind}'))
        for c in cl:
            alle.append(dict(allele_id=c['allele_id'], arm=arm, structure=c['structure'], is_native=c['structure'] == 'home',
                             n_core_blocks=c['structure'].count('|') + 1, switches_median=','.join(map(str, c['sw'])),
                             switches_min=','.join(map(str, c['swmin'])), switches_max=','.join(map(str, c['swmax'])),
                             n_full=len(c['members']), representative_end_id=c['medoid']))
    Mm = R[['end_id', 'arm', 'set', 'sample', 'structure', 'n_core_blocks', 'sub_len', 'full', 'defining']].assign(
        switches=[','.join(map(str, s)) for s in R.switches]).merge(pd.DataFrame(mem), on='end_id', how='left')
    C.write_tsv_atomic(Mm, f'{C.CW}/sa_members.tsv'); C.write_tsv_atomic(pd.DataFrame(alle), f'{C.CW}/sa_alleles.tsv')
    C.write_tsv_atomic(pd.DataFrame(cb), f'{C.CW}/core_blocks.tsv')
    D = pd.DataFrame(diag); C.write_tsv_atomic(D, f'{C.CW}/sa_tiebreak_diagnostics.tsv')
    sd = D[D.kind == 'defining_self']; nd = D[D.kind != 'defining_self']
    print(f"tie-break: defining ends re-assigned to their own allele: nearest-median {(sd.nearest == sd.own).mean():.4f}, "
          f"legacy most-frequent {(sd.legacy == sd.own).mean():.4f} (n = {len(sd)}); non-defining ends whose allele differs between rules: "
          f"{int((nd.nearest != nd.legacy).sum())} of {len(nd)}")
    print('alleles', len(alle), Mm.assignment.value_counts(dropna=False).to_dict())

if __name__ == '__main__':
    main()
