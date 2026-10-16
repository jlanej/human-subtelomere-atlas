"""Seeded, stratified subset for the boundary-shift stability test (redesign acceptance). Strata: previous-build structural-allele tier
(A/B/C) x nativity (native 'home' structure vs non-native) x arm class (duplicon-rich arms 6p, 8p, 9q, 11p, 15q, 16p, 16q, 19p vs
others). Eligible: full-length Release 2 ends that defined an allele in the previous build. Up to PER_STRATUM ends per stratum
(all if fewer). Previous build = work/catalogue/v0.1_prefix/release (or ATLAS_PREVIOUS_RELEASE).
Output: work/catalogue/stability_subset.tsv"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
PER_STRATUM, SEED = 110, 20261017
DUP = {'6p', '8p', '9q', '11p', '15q', '16p', '16q', '19p'}
PREV = os.environ.get('ATLAS_PREVIOUS_RELEASE', f'{C.CW}/v0.1_prefix/release')
E = pd.read_csv(f'{PREV}/ends.tsv', sep='\t', low_memory=False)
SA = pd.read_csv(f'{PREV}/alleles_structural.tsv', sep='\t').set_index('allele_id')
x = E[(E.set == 'HPRC_R2') & (E.sa_assignment == 'full')].copy()
x['tier'] = x.structural_allele_id.map(SA.tier); x['native'] = x.structural_allele_id.map(SA.is_native)
x['arm_class'] = np.where(x.catalogue_arm.isin(DUP), 'duplicon_rich', 'other')
rng = np.random.default_rng(SEED); out = []
for k, g in x.groupby(['tier', 'native', 'arm_class'], sort=True):
    ids = sorted(g.end_id); take = ids if len(ids) <= PER_STRATUM else sorted(rng.choice(ids, PER_STRATUM, replace=False).tolist())
    out.append(g[g.end_id.isin(take)])
S = pd.concat(out)[['end_id', 'catalogue_arm', 'tier', 'native', 'arm_class', 'structural_allele_id', 'core_structure', 'form']]
S.to_csv(f'{C.CW}/stability_subset.tsv', sep='\t', index=False)
print(len(S)); print(S.groupby(['tier', 'native', 'arm_class']).size().to_string())
