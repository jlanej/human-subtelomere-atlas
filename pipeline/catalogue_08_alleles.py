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
  -> the most frequent ('partial_ambiguous', candidates listed); none -> 'unassigned_partial'. Full-length pedigree_aux ends and
  full-length ends with a low-confidence/conflicting arm call are matched against the defined alleles of their catalogue arm
  (same structure, every switch within SW_TOL): 'aux' / 'arm_lowconf', else 'unassigned_aux' / 'unassigned_arm_lowconf'.
Representative (medoid): the full member with the smallest summed Chebyshev distance to the other members (ties: end_id).
Outputs: work/catalogue/sa_members.tsv, sa_alleles.tsv, core_blocks.tsv
"""
import os, sys, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
MIN_BLOCK = int(os.environ.get('ATLAS_MIN_BLOCK', '10000')); SW_TOL = int(os.environ.get('ATLAS_SW_TOL', '10000'))
FULL_SUB = 80_000
DEF_SETS = {'HPRC_R2', 'CHM13', 'HG002_Q100'}

def core(bl):
    bl = [dict(source='home' if b.is_home else b.source, hsrc=b.source, s=float(b.sub_start), e=float(b.sub_end)) for b in bl]
    while len(bl) > 1:
        L = [b['e'] - b['s'] for b in bl]; i = int(np.argmin(L))
        if L[i] >= MIN_BLOCK: break
        if i == 0: bl[1]['s'] = bl[0]['s']
        elif i == len(bl) - 1: bl[-2]['e'] = bl[-1]['e']
        else:
            m = (bl[i]['s'] + bl[i]['e']) / 2; bl[i - 1]['e'] = m; bl[i + 1]['s'] = m
        del bl[i]
        j = 1
        while j < len(bl):
            if bl[j]['source'] == bl[j - 1]['source']: bl[j - 1]['e'] = bl[j]['e']; del bl[j]
            else: j += 1
    toks = [b['source'] if b['source'] in ('home', 'none') else C.arm_of(b['source']) for b in bl]
    return bl, '|'.join(toks), [int(round(b['s'])) for b in bl[1:]]

def prefix_compatible(struct, sw, a_struct, a_sw, Lobs):
    """end (struct, sw) observed up to Lobs; allele (a_struct, a_sw median). Allele truncated at Lobs must equal the end."""
    at = a_struct.split('|'); keep = [s for s in a_sw if s < Lobs - MIN_BLOCK]
    if struct.split('|') != at[:len(keep) + 1]: return False
    return len(sw) == len(keep) and all(abs(a - b) <= SW_TOL for a, b in zip(sw, keep))

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
    R = pd.DataFrame(rows); alle, mem = [], []
    for arm, g in R.groupby('arm', sort=True):
        d = g[g.full & g.set.isin(DEF_SETS) & g.defining]; cl = []
        for st, h in d.groupby('structure', sort=True):
            ids = list(h.end_id); X = np.array([s for s in h.switches], float).reshape(len(h), -1)
            if len(h) == 1 or X.shape[1] == 0: lab = np.ones(len(h), int)
            else: lab = fcluster(linkage(X, 'complete', metric='chebyshev'), t=SW_TOL, criterion='distance')
            for l in sorted(set(lab)):
                ix = np.where(lab == l)[0]; Xi = X[ix]
                if len(ix) > 1 and Xi.shape[1]:
                    D = np.abs(Xi[:, None, :] - Xi[None, :, :]).max(2).sum(1); med = sorted(range(len(ix)), key=lambda j: (D[j], ids[ix[j]]))[0]
                else: med = sorted(range(len(ix)), key=lambda j: ids[ix[j]])[0]
                cl.append(dict(structure=st, members=[ids[j] for j in ix], sw=np.median(Xi, 0).round().astype(int).tolist() if Xi.shape[1] else [],
                               swmin=Xi.min(0).astype(int).tolist() if Xi.shape[1] else [], swmax=Xi.max(0).astype(int).tolist() if Xi.shape[1] else [],
                               medoid=ids[ix[med]]))
        cl.sort(key=lambda c: (-len(c['members']), c['structure'], c['sw']))
        for r, c in enumerate(cl, 1):
            c['allele_id'] = f'SA-{arm}-{r:02d}'
            for m in c['members']: mem.append(dict(end_id=m, allele_id=c['allele_id'], assignment='full'))
        # partial and auxiliary ends
        for x in g[~g.end_id.isin(d.end_id)].itertuples():
            if x.full:
                cand = [c for c in cl if c['structure'] == x.structure and len(c['sw']) == len(x.switches)
                        and all(abs(a - b) <= SW_TOL for a, b in zip(x.switches, c['sw']))]
                kind = 'aux' if x.set == 'pedigree_aux' else 'arm_lowconf'
            else:
                cand = [c for c in cl if prefix_compatible(x.structure, x.switches, c['structure'], c['sw'], x.sub_len)]
                kind = 'partial' if x.set != 'pedigree_aux' else 'aux_partial'
            if cand:
                swd = lambda c: max((abs(a - b) for a, b in zip(c['sw'], x.switches)), default=0)
                best = sorted(cand, key=lambda c: (-len(c['members']), swd(c), c['allele_id']))[0]
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
    print('alleles', len(alle), Mm.assignment.value_counts(dropna=False).to_dict())

if __name__ == '__main__':
    main()
