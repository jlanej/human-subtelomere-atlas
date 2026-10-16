"""Step 5: per-record classification table for all Tan et al. 2026 records, and the same-haplotype crosswalk with the catalogue.
Per Tan record (work/tan2026/tan_classified.tsv): S5/S10/S20 states and families (track-based mosaic, tan2026_02_mosaic.py),
  form = arm>S20 family when S20 is foreign (as catalogue_09_release), core structure and switches (subtel_atlas.alleles.core:
  blocks < 10 kb dropped), nearest structural allele with the catalogue matcher (subtel_atlas.alleles.match_allele via
  tan2026_common.assign), tier (and nativity_class when the release provides it) of that allele.
Window emulation: a catalogue window is 100 kb from the contig tip, i.e. 100 kb minus the telomere (median 92.3 kb of
  subtelomere); Tan records give 100 kb of subtelomere. Tan mosaic blocks are therefore clipped before the core computation
  (tan2026_common.core_of): for crosswalk pairs at the paired catalogue end's own subtelomere length; for all records at a length
  drawn (deterministically, seeded by the record name) from the Release 2 subtelomere lengths of the arm (sub_len_eff), so that
  the cohort records see the catalogue's distribution of window ends (alleles whose last switch lies within 10 kb of the window
  end can only be called on ends with a short telomere). Sensitivity variants: clip at the arm median (*_medclip) and unclipped
  100 kb (*_100k).
Crosswalk (data/validation/tan2026_crosswalk.tsv): one row per key pair from tan2026_04_pairs.py (catalogue named end vs the Tan
  HPRCY2 record of the same sample and arm that matches it best, either haplotype label; hap_label = same/opposite label), with alignment metrics and the agreement of every catalogue call.
  S5/S10/S20 agreement: state and paralog family (primary; family = the label of a foreign window) and state and best single
  CHM13 end (secondary; the best end flips among near-identical paralogs).
  verified = identity >= 0.99 and >= 50% of the catalogue subtelomere aligned (same haplotype sequence);
  structure_agree_common = core structures equal after both are truncated to the shorter observed subtelomere;
  switch deltas: Tan minus catalogue switch position, raw (each in its own boundary-anchored coordinates) and after projecting
  the catalogue switch through the pairwise alignment (removes boundary-call differences);
  sa_agree = same structural allele ID; sa_compatible = the catalogue allele is among the Tan record's candidate alleles.
Also: arm labels of Tan records versus the catalogue arm for unnamed-contig Release 2 ends (best match from tan2026_04_pairs.py).
Outputs: work/tan2026/tan_classified.tsv, data/validation/tan2026_crosswalk.tsv, work/tan2026/summary_crosswalk.json,
         work/tan2026/switch_deltas_primary.tsv
Usage  : python tan2026_05_compare.py
"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
import catalogue_common as C

STATES = ('S5', 'S10', 'S20')

def nativity(SA):
    """allele -> nativity class: the release column nativity_class when present, else native/non-native from is_native"""
    if 'nativity_class' in SA.columns: return dict(zip(SA.allele_id, SA.nativity_class.fillna('')))
    return dict(zip(SA.allele_id, np.where(SA.is_native.astype(bool), 'native', 'non_native')))

def classify():
    out = f'{T.TW}/tan_classified.tsv'; SA0 = pd.read_csv(f'{T.REPO}/data/catalogue/alleles_structural.tsv', sep='\t', low_memory=False)
    deps = [f'{T.TW}/tan_tracks_ends.tsv', f'{T.TW}/tan_tracks_ends.tsv.families_sha256', os.path.join(T.HERE, 'tan2026_common.py'), os.path.abspath(__file__), f'{T.REPO}/data/catalogue/alleles_structural.tsv',
            f'{T.REPO}/data/catalogue/ends.tsv']
    if os.path.exists(out) and os.path.getmtime(out) > max(os.path.getmtime(f) for f in deps):
        return pd.read_csv(out, sep='\t', low_memory=False), SA0
    W = pd.read_csv(f'{T.TW}/tan_windows.tsv', sep='\t')
    ME = pd.read_csv(f'{T.TW}/tan_tracks_ends.tsv', sep='\t', low_memory=False)
    MB = pd.read_csv(f'{T.TW}/tan_tracks_blocks.tsv', sep='\t', low_memory=False)
    SA, BY = T.load_alleles(); tier = dict(zip(SA.allele_id, SA.tier))
    E = pd.read_csv(f'{T.REPO}/data/catalogue/ends.tsv', sep='\t', low_memory=False, usecols=['set', 'catalogue_arm', 'sub_len'])
    R2 = E[E.set == 'HPRC_R2']; Lm = R2.groupby('catalogue_arm').sub_len.median().to_dict()
    Ld = {a: np.sort(g.sub_len.values) for a, g in R2.groupby('catalogue_arm')}
    def draw(e, a, full):   # deterministic draw from the arm's Release 2 subtelomere lengths (seed = record name)
        if a not in Ld: return full
        rng = np.random.default_rng(int(C.sha256_str(e)[:12], 16)); return min(full, int(rng.choice(Ld[a])))
    wl = ME.end_id.map(dict(zip(W.end_id, W.window_len)))
    ME = ME.assign(sub_len_eff=[draw(e, a, w - s0) for e, w, s0, a in zip(ME.end_id, wl, ME.s0, ME.catalogue_arm)],
                   sub_len_median=[min(w - s0, Lm.get(a, w - s0)) for w, s0, a in zip(wl, ME.s0, ME.catalogue_arm)])
    clipL = dict(zip(ME.end_id, ME.sub_len_eff)); clipM = dict(zip(ME.end_id, ME.sub_len_median))
    G = list(MB.groupby('end_id', sort=False))
    core = {e: T.core_of(g, clip=clipL[e]) for e, g in G}
    coreM = {e: T.core_of(g, clip=clipM[e]) for e, g in G}
    core100 = {e: T.core_of(g) for e, g in G}
    D = W.merge(ME[['end_id', 'catalogue_arm', 's0', 'n_blocks', 'profile', 'structure', 'switches_sub', 'switch_res']
                   + [f'{s}_{k}' for s in STATES for k in ('state', 'best_end', 'family')] + ['S20_best_ident', 'sub_len_eff', 'sub_len_median']], on='end_id', how='left')
    D['sub_len'] = D.window_len - D.s0
    import pysam, re
    fa = pysam.FastaFile(f'{T.TW}/tan_tips.fa')   # assembly gaps (N runs) in the classified window
    D['n_gap_bases'] = [sum(m.end() - m.start() for m in re.finditer('N+', fa.fetch(e).upper())) for e in D.end_id]
    D['core_structure'] = [core[e][0] if e in core else '' for e in D.end_id]
    D['core_switches_sub'] = [','.join(map(str, core[e][1])) if e in core else '' for e in D.end_id]
    D['n_core_blocks'] = [core[e][2] if e in core else 0 for e in D.end_id]
    D['core_structure_100k'] = [core100[e][0] if e in core100 else '' for e in D.end_id]
    D['core_switches_sub_100k'] = [','.join(map(str, core100[e][1])) if e in core100 else '' for e in D.end_id]
    D['form'] = np.where(D.S20_state == 'foreign', D.tan_arm + '>' + D.S20_family.fillna('').astype(str), '')
    asg = [T.assign(a, s, T.ints(sw), L, BY) if s else ('', 'no_mosaic', '') for a, s, sw, L in
           zip(D.catalogue_arm, D.core_structure, D.core_switches_sub, D.sub_len_eff)]
    D['structural_allele_id'] = [a[0] for a in asg]; D['sa_assignment'] = [a[1] for a in asg]; D['sa_candidates'] = [a[2] for a in asg]
    D['structural_allele_id_100k'] = [T.assign(a, s, T.ints(sw), L, BY)[0] if s else '' for a, s, sw, L in
                                      zip(D.catalogue_arm, D.core_structure_100k, D.core_switches_sub_100k, D.sub_len)]
    D['structural_allele_id_medclip'] = [T.assign(a, coreM[e][0], coreM[e][1], L, BY)[0] if e in coreM and coreM[e][0] else ''
                                         for e, a, L in zip(D.end_id, D.catalogue_arm, D.sub_len_median)]
    D['sa_tier'] = D.structural_allele_id.map(tier).fillna('')
    D['sa_nativity_class'] = D.structural_allele_id.map(nativity(SA)).fillna('')
    C.write_tsv_atomic(D, f'{T.TW}/tan_classified.tsv'); return D, SA

def switch_cmp(r):
    """per-switch deltas for pairs with identical raw mosaic structures"""
    a, b = T.ints(r.switches_sub_c), T.ints(r.switches_sub_t); pr = str(r.switches_projected) if isinstance(r.switches_projected, str) else ''
    pj = [int(x) if x not in ('', 'nan') else None for x in pr.split(',')] if pr else []
    ra = str(r.switch_res_c).split(',') if isinstance(r.switch_res_c, str) else []; rb = str(r.switch_res_t).split(',') if isinstance(r.switch_res_t, str) else []
    keep = [i for i, y in enumerate(b) if y < r.sub_len_c]           # Tan switches beyond the catalogue window are not compared
    st_t = '|'.join(str(r.structure_t).split('|')[:len(keep) + 1]); b = [b[i] for i in keep]; rb = [rb[i] for i in keep if i < len(rb)]
    if r.structure_c != st_t or not a or len(a) != len(b): return []
    out = []
    for i, (x, y) in enumerate(zip(a, b)):
        p = pj[i] if i < len(pj) else None
        out.append(dict(d_raw=y - x, d_proj=(y - p) if p is not None else np.nan, base=(i < len(ra) and i < len(rb) and ra[i] == 'base' and rb[i] == 'base')))
    return out

def summ(X, swl):
    n = len(X); f = lambda m: round(float(m.mean()), 4) if len(m) else None
    o = dict(n_pairs=int(n))
    if n == 0: return o
    o.update(identity_median=round(float(X.ident.median()), 6), frac_identity_ge_0999=f(X.ident >= 0.999), frac_identity_ge_09999=f(X.ident >= 0.9999),
             sub_aligned_frac_median=round(float(X.sub_aligned_frac.median()), 4), frac_sub_aligned_ge_095=f(X.sub_aligned_frac >= 0.95))
    bd = X.boundary_delta.dropna().abs()
    o.update(n_boundary_projected=int(len(bd)), frac_boundary_exact=f(bd == 0), frac_boundary_le_10bp=f(bd <= 10), frac_boundary_le_100bp=f(bd <= 100))
    for s in STATES:
        o[f'{s}_state_agree'] = f(X[f'{s}_state_c'] == X[f'{s}_state_t'])
        o[f'{s}_state_and_family_agree'] = f((X[f'{s}_state_c'] == X[f'{s}_state_t']) & (X[f'{s}_family_c'].fillna('').astype(str) == X[f'{s}_family_t'].fillna('').astype(str)))
        o[f'{s}_state_and_end_agree'] = f((X[f'{s}_state_c'] == X[f'{s}_state_t']) & (X[f'{s}_best_end_c'].fillna('') == X[f'{s}_best_end_t'].fillna('')))
    o['S20_foreign_c'] = int((X.S20_state_c == 'foreign').sum()); o['S20_foreign_t'] = int((X.S20_state_t == 'foreign').sum())
    o['form_agree'] = f(X.form_c.fillna('') == X.form_t.fillna(''))
    o['core_structure_agree'] = f(X.core_structure_c == X.core_structure_t); o['core_structure_agree_common_length'] = f(X.structure_agree_common)
    sa = X[(X.sa_assignment_c == 'full') & (X.sa_assignment_t == 'full')]
    o.update(n_sa_compared=int(len(sa)), sa_agree=f(sa.sa_agree), sa_compatible=f(sa.sa_compatible),
             frac_tan_unassigned=f(X.sa_assignment_t.isin(['unassigned', 'unassigned_partial'])))
    S = pd.DataFrame([s for l in swl for s in l]) if any(swl) else pd.DataFrame(columns=['d_raw', 'd_proj', 'base'])
    o['n_switches_compared'] = int(len(S))
    if len(S):
        Sb = S[S.base]; o['n_switches_base_both'] = int(len(Sb))
        for nm, Y in (('all', S), ('base', Sb)):
            if len(Y):
                o[f'switch_{nm}_frac_raw_exact'] = f(Y.d_raw == 0); o[f'switch_{nm}_frac_raw_le_100bp'] = f(Y.d_raw.abs() <= 100)
                o[f'switch_{nm}_frac_proj_exact'] = f(Y.d_proj.dropna() == 0); o[f'switch_{nm}_frac_proj_le_10bp'] = f(Y.d_proj.dropna().abs() <= 10)
                o[f'switch_{nm}_frac_proj_le_1kb'] = f(Y.d_proj.dropna().abs() <= 1000); o[f'switch_{nm}_median_abs_proj'] = float(Y.d_proj.dropna().abs().median())
    return o

def main():
    D, SA = classify()
    E = pd.read_csv(f'{T.REPO}/data/catalogue/ends.tsv', sep='\t', low_memory=False); E = E[E.set == 'HPRC_R2']
    E['s0'] = E.junk_len + E.telo_len; E['S20_best_ident'] = E['S20_best_ident']
    P = pd.read_csv(f'{T.TW}/pairs.tsv', sep='\t', dtype={'switches_projected': str, 'core_switches_projected': str}).drop(columns=['pansn_hap', 'tan_hap'])
    cc = ['end_id', 'sample', 'pansn_hap', 'catalogue_arm', 'assembler', 's0', 'sub_len', 'structure', 'switches_sub', 'switch_res', 'core_structure',
          'core_switches_sub', 'form', 'structural_allele_id', 'sa_assignment', 'sa_candidates', 'flags'] + [f'{s}_{k}' for s in STATES for k in ('state', 'best_end', 'family')]
    tc = ['end_id', 'assembler', 'hap', 'has_telomere', 'telo_len', 's0', 'sub_len', 'structure', 'switches_sub', 'switch_res', 'core_structure',
          'core_switches_sub', 'form', 'structural_allele_id', 'sa_assignment', 'sa_candidates', 'sa_tier'] + [f'{s}_{k}' for s in STATES for k in ('state', 'best_end', 'family')]
    X = P.merge(E[cc].add_suffix('_c').rename(columns={'end_id_c': 'end_id'}), on='end_id').merge(
        D[tc].add_suffix('_t').rename(columns={'end_id_t': 'tan_end_id'}), on='tan_end_id')
    # Tan core structure and allele re-derived with the blocks clipped at the paired catalogue end's subtelomere length
    MB = pd.read_csv(f'{T.TW}/tan_tracks_blocks.tsv', sep='\t', low_memory=False); MB = MB[MB.end_id.isin(set(X.tan_end_id))]
    g = dict(list(MB.groupby('end_id', sort=False))); _, BY = T.load_alleles(); re = []
    for r in X.itertuples():
        L = min(r.sub_len_t, r.sub_len_c)
        if r.tan_end_id not in g: re.append(('', '', '', 'no_mosaic', '')); continue
        st, sw, _n = T.core_of(g[r.tan_end_id], clip=L); a = T.assign(r.catalogue_arm_c, st, sw, L, BY)
        re.append((st, ','.join(map(str, sw)), a[0], a[1], a[2]))
    for i, c in enumerate(['core_structure_t', 'core_switches_sub_t', 'structural_allele_id_t', 'sa_assignment_t', 'sa_candidates_t']):
        X[c] = [v[i] for v in re]
    X['tier_c'] = X.structural_allele_id_c.map(dict(zip(SA.allele_id, SA.tier))).fillna('')
    X['nativity_class_c'] = X.structural_allele_id_c.map(nativity(SA)).fillna('')
    X['comparison'] = X.assembler_c + '_vs_' + X.assembler_t
    X['independent_assembler'] = X.assembler_c != X.assembler_t
    X['verified'] = (X.ident >= 0.99) & (X.sub_aligned_frac >= 0.5)
    for c in ('core_structure_c', 'core_structure_t', 'structure_c', 'structure_t'): X[c] = X[c].fillna('')
    sac = []
    for r in X.itertuples():
        L = min(r.sub_len_c, r.sub_len_t)
        a = T.truncate(r.core_structure_c, T.ints(r.core_switches_sub_c), L); b = T.truncate(r.core_structure_t, T.ints(r.core_switches_sub_t), L)
        sac.append(a[0] == b[0])
    X['structure_agree_common'] = sac
    X['sa_agree'] = (X.structural_allele_id_c.fillna('') == X.structural_allele_id_t.fillna('')) & (X.structural_allele_id_c.fillna('') != '')
    X['sa_compatible'] = [isinstance(a, str) and a != '' and a in str(c).split(',') for a, c in zip(X.structural_allele_id_c, X.sa_candidates_t)]
    swl = [switch_cmp(r) for r in X.itertuples()]
    X['n_switch_compared'] = [len(s) for s in swl]
    X['max_abs_switch_delta_raw'] = [max((abs(d['d_raw']) for d in s), default=np.nan) for s in swl]
    X['max_abs_switch_delta_projected'] = [max((abs(d['d_proj']) for d in s if d['d_proj'] == d['d_proj']), default=np.nan) for s in swl]
    fl = X.flags_c.fillna('')
    X['artefact_candidate'] = fl.str.contains('artefact_candidate'); X['junk_flagged'] = fl.str.contains('junk_')
    X = X.rename(columns={'sample_c': 'sample', 'pansn_hap_c': 'hap', 'catalogue_arm_c': 'arm'})
    out_cols = ['sample', 'hap', 'arm', 'end_id', 'tan_end_id', 'assembler_c', 'assembler_t', 'comparison', 'independent_assembler', 'hap_label', 'n_candidates',
                'ident', 'nm_per_kb', 'sub_aligned_frac', 'tan_sub_aligned_frac', 'n_aln', 'verified', 's0_c', 's0_t', 'has_telomere_t', 'boundary_delta',
                'diag_offset'] + [f'{s}_{k}_{w}' for s in STATES for k in ('state', 'family', 'best_end') for w in 'ct'] + [
                'form_c', 'form_t', 'core_structure_c', 'core_structure_t', 'core_switches_sub_c', 'core_switches_sub_t', 'core_switches_projected',
                'structure_agree_common', 'structure_c', 'structure_t', 'switches_sub_c', 'switches_sub_t', 'switches_projected', 'switch_res_c',
                'switch_res_t', 'n_switch_compared', 'max_abs_switch_delta_raw', 'max_abs_switch_delta_projected', 'structural_allele_id_c',
                'sa_assignment_c', 'tier_c', 'nativity_class_c', 'structural_allele_id_t', 'sa_assignment_t', 'sa_candidates_t', 'sa_agree', 'sa_compatible',
                'artefact_candidate', 'junk_flagged', 'flags_c']
    X = X.sort_values(['sample', 'hap', 'arm', 'end_id', 'tan_end_id'])
    C.write_tsv_atomic(X.rename(columns={'end_id': 'catalogue_end_id'})[[c if c != 'end_id' else 'catalogue_end_id' for c in out_cols]],
                       f'{T.VAL}/tan2026_crosswalk.tsv')
    # summaries: primary set = verified pairs whose Tan record carries a telomere
    X['swl'] = swl; prim = X[X.verified & X.has_telomere_t.astype(bool)]
    SWD = pd.DataFrame([dict(end_id=r.end_id, comparison=r.comparison, **d) for r in prim.itertuples() for d in r.swl])
    SWD.to_csv(f'{T.TW}/switch_deltas_primary.tsv', sep='\t', index=False)   # per-switch deltas of primary pairs (figure panel b)
    S = dict(n_key_pairs=int(len(X)), n_keys=int(X[['sample', 'hap', 'arm']].drop_duplicates().shape[0]), n_verified=int(X.verified.sum()),
             n_tan_no_telomere=int((~X.has_telomere_t.astype(bool)).sum()), n_primary=int(len(prim)), n_samples=int(X['sample'].nunique()),
             hap_label_relation=X[X.verified].hap_label.value_counts().to_dict(),
             self_assignment_baseline=self_baseline(E))
    S['all_pairs'] = summ(X, list(X.swl)); S['primary'] = summ(prim, list(prim.swl))
    S['by_comparison'] = {k: summ(g, list(g.swl)) for k, g in prim.groupby('comparison')}
    bd = prim.boundary_delta.abs()
    S['by_boundary_delta'] = {k: summ(prim[m], list(prim[m].swl)) for k, m in
                              (('exact', bd == 0), ('1-100bp', (bd > 0) & (bd <= 100)), ('101bp-1kb', (bd > 100) & (bd <= 1000)),
                               ('>1kb', bd > 1000), ('not_projected', bd.isna()))}
    S['by_comparison_and_boundary'] = {f'{c}|{"exact" if e else "differs"}': summ(g, list(g.swl)) for (c, e), g in
                                       prim.assign(_e=bd == 0).groupby(['comparison', '_e'])}
    S['by_tier'] = {(k or 'none'): summ(g, list(g.swl)) for k, g in prim.groupby('tier_c')}
    S['by_nativity_class'] = {(k or 'none'): summ(g, list(g.swl)) for k, g in prim.groupby('nativity_class_c')}
    S['nativity_class_source'] = 'alleles_structural.nativity_class' if 'nativity_class' in SA.columns else 'alleles_structural.is_native (nativity_class absent)'
    S['by_arm'] = {k: summ(g, list(g.swl)) for k, g in prim.groupby('arm')}
    top = prim.form_c.fillna('').value_counts(); top = [f for f in top.index if f][:25]
    S['by_form_top25'] = {k: summ(prim[prim.form_c == k], list(prim[prim.form_c == k].swl)) for k in top}
    S['by_form_class'] = {k: summ(g, list(g.swl)) for k, g in prim.groupby(np.where(prim.form_c.fillna('') == '', 'no_S20_foreign', 'S20_foreign'))}
    S['artefact_candidates'] = [dict(catalogue_end_id=r.end_id, tan_end_id=r.tan_end_id, ident=r.ident, sub_aligned_frac=r.sub_aligned_frac,
                                     boundary_delta=None if pd.isna(r.boundary_delta) else int(r.boundary_delta),
                                     core_structure_c=r.core_structure_c, core_structure_t=r.core_structure_t, sa_c=r.structural_allele_id_c,
                                     sa_t=r.structural_allele_id_t, form_c=r.form_c if isinstance(r.form_c, str) else '',
                                     form_t=r.form_t if isinstance(r.form_t, str) else '') for r in X[X.artefact_candidate].itertuples()]
    art = set(E.end_id[E['flags'].fillna('').str.contains('artefact_candidate')])
    S['artefact_candidates_without_pair'] = sorted(art - set(X.end_id))
    S['junk_flagged'] = summ(X[X.junk_flagged], list(X[X.junk_flagged].swl))
    S['junk_flagged_primary'] = summ(prim[prim.junk_flagged], list(prim[prim.junk_flagged].swl))
    S['unnamed_arm_check'] = unnamed_arms(E, D)
    json.dump(S, open(f'{T.TW}/summary_crosswalk.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps({k: S['primary'][k] for k in ('n_pairs', 'identity_median', 'S20_state_agree', 'core_structure_agree', 'sa_agree')}))

def self_baseline(E):
    SA, BY = T.load_alleles(); r = E[E.sa_assignment == 'full']
    re = [T.assign(a, s if isinstance(s, str) else 'home', T.ints(sw), L, BY)[0] for a, s, sw, L in zip(r.catalogue_arm, r.core_structure, r.core_switches_sub, r.sub_len)]
    re2 = [T.assign(a, s if isinstance(s, str) else 'home', T.ints(sw), L, BY, rule='most_frequent')[0] for a, s, sw, L in zip(r.catalogue_arm, r.core_structure, r.core_switches_sub, r.sub_len)]
    return dict(n=int(len(r)), frac_rule_returns_own_allele=round(float(np.mean(np.array(re) == r.structural_allele_id.values)), 4),
                frac_catalogue_most_frequent_rule_returns_own_allele=round(float(np.mean(np.array(re2) == r.structural_allele_id.values)), 4))

def unnamed_arms(E, D):
    BM = pd.read_csv(f'{T.TW}/best_match.tsv', sep='\t').merge(E[['end_id', 'catalogue_arm', 'arm_confidence', 'arm_source']], on='end_id')
    BM = BM.merge(D[['end_id', 'tan_arm']].rename(columns={'end_id': 'best_tan'}), on='best_tan', how='left')
    BM['matched'] = (BM.best_ident >= 0.99) & (BM.best_sub_frac >= 0.5)
    BM['tan_group'] = [C.group_of(a) if isinstance(a, str) else '' for a in BM.tan_arm]
    BM['agree'] = BM.tan_group == BM.catalogue_arm
    BM.to_csv(f'{T.TW}/best_match_arms.tsv', sep='\t', index=False)
    out = {}
    for k, g in BM.groupby(np.where(BM.arm_confidence == 'named', 'named(control)', 'unnamed:' + BM.arm_confidence.astype(str))):
        m = g[g.matched]
        out[k] = dict(n_ends=int(len(g)), n_matched=int(len(m)), arm_agree=round(float(m.agree.mean()), 4) if len(m) else None,
                      n_disagree=int((~m.agree).sum()))
    m = BM[BM.matched & (BM.arm_confidence != 'named')]
    out['unnamed_all'] = dict(n_ends=int((BM.arm_confidence != 'named').sum()), n_matched=int(len(m)),
                              arm_agree=round(float(m.agree.mean()), 4) if len(m) else None, n_disagree=int((~m.agree).sum()))
    return out

if __name__ == '__main__':
    main()
