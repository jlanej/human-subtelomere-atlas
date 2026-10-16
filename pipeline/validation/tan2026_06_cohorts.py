"""Step 6: replication of catalogue forms and structural alleles in the non-HPRC cohorts of Tan et al. 2026 (HGSVC3, CPC, APC, JASA).
Records: work/tan2026/tan_classified.tsv (tan2026_05_compare.py). Primary set = records with a telomere block at the tip (the catalogue
  is built from telomeric tips only), >= 80 kb of subtelomere (catalogue 'full' length) and no assembly gap (N) in the window.
De-duplication / relatedness: samples present in the catalogue (HPRC Release 2, HG002, pedigree) or in the Tan HPRCY2 set, or
  aliases of those samples (data/metadata/samples.tsv alternative_id, e.g. NA24385 = HG002), are removed;
  children of trios whose parent is in the same cohort are removed (1000 Genomes pedigree, work/tan2026/igsr_3202_ped.txt;
  APC sample 'apr-s', whose parents 'apr-m' and 'apr-f' are in the set, by its label).
Groups: each cohort overall; HGSVC3 by 1000 Genomes superpopulation; JASA split into Japanese (1000 Genomes JPT samples) and Saudi
  ('ksa' samples). HPRC comparator: Release 2 catalogue frequency overall (ALL) and, where the group has a matching 1000 Genomes
  superpopulation (CPC -> EAS, JASA Japanese -> EAS, HGSVC3 per superpopulation), that superpopulation
  (data/catalogue/allele_frequencies_by_superpop.tsv). APC and JASA Saudi have no matching superpopulation (ALL only).
Window: Tan records are classified with a window end drawn from the Release 2 subtelomere lengths of the arm (tan2026_05_compare);
  window_sensitivity in the summary repeats the pooled comparison with an arm-median clip and with unclipped 100 kb windows.
Structural alleles: tier A/B alleles (data/catalogue/alleles_structural.tsv); count = primary records of the arm assigned to the allele
  (catalogue compatibility rules, nearest-median tie-break, full-length
  assignment); denominator = all primary records of the arm (assigned or not). 95% Wilson CI;
  two-sided Fisher exact test against the HPRC counts; Benjamini-Hochberg q within each table kind.
Forms: arm>paralog family (or arm) of a foreign S20 call (catalogue revision-2 form); forms carried by >= 20 Release 2 ends; denominator = records of the arm
  with an S20 call (state not 'none'), as for Release 2.
Novel structures: primary non-HPRC records with no matching allele (full-length 'unassigned'), clustered per (arm, core structure)
  by complete linkage on the switch vector (Chebyshev, 10 kb; as catalogue_08_alleles). known_structure = the structure string
  exists among the arm's catalogue alleles (only the switch positions differ). candidate_allele = >= 2 unrelated samples.
Outputs: data/validation/tan2026_cohort_forms.tsv, data/validation/tan2026_novel_structures.tsv, work/tan2026/summary_cohorts.json
Usage  : python tan2026_06_cohorts.py
"""
import os, sys, json, numpy as np, pandas as pd
from scipy.stats import fisher_exact
from scipy.cluster.hierarchy import linkage, fcluster
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
import catalogue_common as C
import catalogue_08_alleles as C08
KEY_ALLELES = ('SA-20q-01', 'SA-3q-01', 'SA-6p-01', 'SA-2q-02', 'SA-8p-02')   # revision-2 forms 20q>PF04, 3q>PF03, 6p>PF03, 2q>16q, 8p>PF02

def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = ~np.isnan(p); pv = p[ok]
    if len(pv):
        o = np.argsort(pv); r = pv[o] * len(pv) / (np.arange(len(pv)) + 1); r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty(len(pv)); qq[o] = np.minimum(r, 1); q[ok] = qq
    return q

def cohort_records():
    D = pd.read_csv(f'{T.TW}/tan_classified.tsv', sep='\t', low_memory=False)
    E = pd.read_csv(f'{T.REPO}/data/catalogue/ends.tsv', sep='\t', low_memory=False, usecols=['set', 'sample'])
    M = pd.read_csv(f'{T.REPO}/data/metadata/samples.tsv', sep='\t')
    hprc = set(E['sample']) | set(D['sample'][D.cohort == 'HPRCY2'])   # catalogue samples (incl. HG002, pedigree) and Tan HPRCY2
    hprc |= set(M.alternative_id.dropna()[M['sample'].isin(hprc)])      # aliases, e.g. NA24385 = HG002
    X = D[D.cohort.isin(T.COHORTS_INDEP)].copy(); log = dict(n_records=int(len(X)), n_samples=int(X['sample'].nunique()))
    pedf = f'{T.TW}/igsr_3202_ped.txt'; ped = pd.read_csv(pedf, sep=' ') if os.path.exists(pedf) else pd.DataFrame(columns=['SampleID', 'FatherID', 'MotherID', 'Population', 'Superpopulation'])
    drop_h = sorted(set(X['sample']) & hprc)
    kids = set()
    for c, g in X.groupby('cohort'):
        s = set(g['sample']); m = ped[ped.SampleID.isin(s)]; kids |= set(m.SampleID[(m.FatherID.isin(s)) | (m.MotherID.isin(s))])
    if {'apr-s', 'apr-m', 'apr-f'} <= set(X['sample']): kids.add('apr-s')
    kids = sorted(kids - set(drop_h))
    X = X[~X['sample'].isin(set(drop_h) | set(kids))]
    sp = dict(zip(ped.SampleID, ped.Superpopulation)); pop = dict(zip(ped.SampleID, ped.Population))
    X['superpopulation'] = X['sample'].map(sp).fillna('')
    X['group'] = X.cohort
    X.loc[(X.cohort == 'JASA') & X['sample'].map(pop).eq('JPT'), 'group'] = 'JASA_Japanese'
    X.loc[(X.cohort == 'JASA') & X['sample'].str.lower().str.startswith('ksa'), 'group'] = 'JASA_Saudi'
    log.update(removed_hprc_overlap=drop_h, removed_trio_children=kids, n_samples_after=int(X['sample'].nunique()),
               jasa_unclassified=sorted(set(X['sample'][(X.cohort == 'JASA') & (X.group == 'JASA')])),
               hgsvc3_no_superpop=sorted(set(X['sample'][(X.cohort == 'HGSVC3') & (X.superpopulation == '')])))
    X['primary'] = X.has_telomere.astype(bool) & (X.sub_len >= C08.FULL_SUB) & X.core_structure.notna() & (X.n_gap_bases == 0)
    log['n_excluded_gap'] = int((X.n_gap_bases > 0).sum())
    return X, log

def groups(X):
    P = X[X.primary]
    yield 'ALL_INDEPENDENT', 'ALL', P
    for c in T.COHORTS_INDEP:
        g = P[P.cohort == c]; yield c, {'CPC': 'EAS'}.get(c, 'ALL'), g
    for c in ('JASA_Japanese', 'JASA_Saudi'): yield c, 'EAS' if c == 'JASA_Japanese' else 'ALL', P[P.group == c]
    for s, g in P[P.cohort == 'HGSVC3'].groupby('superpopulation'):
        if s: yield f'HGSVC3_{s}', s, g

def main():
    X, log = cohort_records()
    SA = pd.read_csv(f'{T.REPO}/data/catalogue/alleles_structural.tsv', sep='\t', low_memory=False)
    F = pd.read_csv(f'{T.REPO}/data/catalogue/allele_frequencies_by_superpop.tsv', sep='\t')
    F = F[F.allele_type == 'structural_allele'].set_index(['allele_id', 'superpopulation'])
    DEN = F.reset_index().groupby(['arm', 'superpopulation']).n_arm_ends_assessed.max().to_dict()
    E = pd.read_csv(f'{T.REPO}/data/catalogue/ends.tsv', sep='\t', low_memory=False)
    R2 = E[E.set == 'HPRC_R2']; SP = dict(zip(R2.end_id, R2.superpopulation))
    AB = SA[SA.tier.isin(['A', 'B'])]
    NAT = dict(zip(SA.allele_id, SA.nativity_class.fillna(''))) if 'nativity_class' in SA.columns else \
          dict(zip(SA.allele_id, np.where(SA.is_native.astype(bool), 'native', 'non_native')))
    rows = []
    for gname, sp, g in groups(X):
        narm = g.groupby('catalogue_arm').size(); cnt = g.groupby('structural_allele_id').size()
        for a in AB.itertuples():
            n = int(narm.get(a.arm, 0)); k = int(cnt.get(a.allele_id, 0))
            if n == 0: continue
            for ref in (['ALL'] + ([sp] if sp != 'ALL' else [])):
                if (a.allele_id, ref) in F.index:
                    K, N = int(F.at[(a.allele_id, ref), 'n_haplotype_ends']), int(F.at[(a.allele_id, ref), 'n_arm_ends_assessed'])
                else:
                    K = 0; N = int(DEN.get((a.arm, ref), 0))
                lo, hi = T.wilson(k, n); Lo, Hi = T.wilson(K, N)
                p = fisher_exact([[k, n - k], [K, N - K]])[1] if N else np.nan
                rows.append(dict(kind='structural_allele', group=gname, cohort=gname.split('_')[0] if gname != 'ALL_INDEPENDENT' else gname,
                                 hprc_reference=ref, arm=a.arm, item=a.allele_id, form=a.form if isinstance(a.form, str) else '', structure=a.structure,
                                 tier=a.tier, nativity_class=NAT.get(a.allele_id, ''), n=k, denominator=n, n_samples=int(g[g.structural_allele_id == a.allele_id]['sample'].nunique()),
                                 freq=round(k / n, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4), hprc_n=K, hprc_denominator=N,
                                 hprc_freq=round(K / N, 4) if N else np.nan, hprc_ci95_lo=round(Lo, 4), hprc_ci95_hi=round(Hi, 4), fisher_p=p))
    # forms
    R2f = R2[R2.S20_state.notna() & (R2.S20_state != 'none')]
    fc = R2f.form.fillna('').value_counts(); forms = [f for f in fc.index if f and fc[f] >= 20]
    for gname, sp, g in groups(X):
        gs = g[g.S20_state.notna() & (g.S20_state != 'none')]
        for f in forms:
            arm = f.split('>')[0]; n = int((gs.tan_arm == arm).sum()); k = int((gs.form == f).sum())
            if n == 0: continue
            for ref in (['ALL'] + ([sp] if sp != 'ALL' else [])):
                rr = R2f[R2f.home_arm == arm] if ref == 'ALL' else R2f[(R2f.home_arm == arm) & (R2f.superpopulation == ref)]
                K, N = int((rr.form == f).sum()), len(rr)
                lo, hi = T.wilson(k, n); Lo, Hi = T.wilson(K, N)
                rows.append(dict(kind='form', group=gname, cohort=gname.split('_')[0] if gname != 'ALL_INDEPENDENT' else gname, hprc_reference=ref,
                                 arm=arm, item=f, form=f, structure='', tier='', nativity_class='', n=k, denominator=n, n_samples=int(gs[gs.form == f]['sample'].nunique()),
                                 freq=round(k / n, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4), hprc_n=K, hprc_denominator=N,
                                 hprc_freq=round(K / N, 4) if N else np.nan, hprc_ci95_lo=round(Lo, 4), hprc_ci95_hi=round(Hi, 4),
                                 fisher_p=fisher_exact([[k, n - k], [K, N - K]])[1] if N else np.nan))
    R = pd.DataFrame(rows); R['bh_q'] = np.nan
    for kd in R.kind.unique():
        m = R.kind == kd; R.loc[m, 'bh_q'] = bh(R.loc[m, 'fisher_p'])
    R['fisher_p'] = R.fisher_p.map(lambda v: float(f'{v:.3g}')); R['bh_q'] = R.bh_q.map(lambda v: float(f'{v:.3g}'))
    C.write_tsv_atomic(R, f'{T.VAL}/tan2026_cohort_forms.tsv')
    N = novel(X, SA)
    C.write_tsv_atomic(N, f'{T.VAL}/tan2026_novel_structures.tsv')
    S = dict(log, **cohort_summary(X, R, N))
    S['window_sensitivity'] = window_sensitivity(X, SA, F)
    json.dump(S, open(f'{T.TW}/summary_cohorts.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps({k: S[k] for k in ('n_primary', 'replication_tierA_common')}, default=str))

def window_sensitivity(X, SA, F):
    """pooled independent cohorts vs HPRC ALL for tier A/B alleles under three window emulations of the Tan records:
    sampled (primary; Release 2 length distribution), arm-median clip, unclipped 100 kb"""
    P = X[X.primary]; AB = SA[SA.tier.isin(['A', 'B'])]; Fa = F.xs('ALL', level='superpopulation') if 'ALL' in F.index.get_level_values(1) else None
    out = {}
    for nm, col in (('sampled_R2_lengths', 'structural_allele_id'), ('arm_median_clip', 'structural_allele_id_medclip'), ('unclipped_100kb', 'structural_allele_id_100k')):
        if col not in P.columns: continue
        narm = P.groupby('catalogue_arm').size(); cnt = P.groupby(col).size(); rows = []
        for a in AB.itertuples():
            n = int(narm.get(a.arm, 0)); k = int(cnt.get(a.allele_id, 0))
            if n == 0 or Fa is None or a.allele_id not in Fa.index: continue
            K, Nn = int(Fa.at[a.allele_id, 'n_haplotype_ends']), int(Fa.at[a.allele_id, 'n_arm_ends_assessed'])
            rows.append(dict(native=bool(a.is_native), tier=a.tier, f=k / n, F=K / Nn if Nn else np.nan, p=fisher_exact([[k, n - k], [K, Nn - K]])[1],
                             common_A=a.tier == 'A' and Nn and K / Nn >= 0.05, k=k))
        Rr = pd.DataFrame(rows); Rr['q'] = bh(Rr.p); nn = Rr[~Rr.native]
        out[nm] = dict(n_tests=int(len(Rr)), n_q_lt_0_05=int((Rr.q < 0.05).sum()), non_native_pearson_r=round(float(np.corrcoef(nn.f, nn.F)[0, 1]), 4),
                       common_tierA_observed=f"{int((Rr[Rr.common_A.astype(bool)].k > 0).sum())}/{int(Rr.common_A.astype(bool).sum())}",
                       frac_unassigned=round(float((P[col].fillna('') == '').mean()), 4))
    return out

def novel(X, SA):
    P = X[X.primary & (X.sa_assignment == 'unassigned')]; known = set(zip(SA.arm, SA.structure)); out = []
    sw_by = {(r.arm, r.structure): [] for r in SA.itertuples()}
    for r in SA.itertuples(): sw_by[(r.arm, r.structure)].append((r.allele_id, T.ints(r.switches_median)))
    for (arm, st), h in P.groupby(['catalogue_arm', 'core_structure'], sort=True):
        ids = list(h.end_id); Xs = np.array([T.ints(s) for s in h.core_switches_sub.fillna('')], float).reshape(len(h), -1)
        lab = np.ones(len(h), int) if len(h) == 1 or Xs.shape[1] == 0 else fcluster(linkage(Xs, 'complete', metric='chebyshev'), t=C08.SW_TOL, criterion='distance')
        for l in sorted(set(lab)):
            ix = np.where(lab == l)[0]; g = h.iloc[ix]; Xi = Xs[ix]
            med = np.median(Xi, 0).round().astype(int).tolist() if Xi.shape[1] else []
            near = sorted(((max(abs(a - b) for a, b in zip(sw, med)), aid) for aid, sw in sw_by.get((arm, st), []) if len(sw) == len(med) and med)) or None
            out.append(dict(arm=arm, core_structure=st, n_core_blocks=st.count('|') + 1, switches_median=','.join(map(str, med)),
                            switches_min=','.join(map(str, Xi.min(0).astype(int))) if Xi.shape[1] else '',
                            switches_max=','.join(map(str, Xi.max(0).astype(int))) if Xi.shape[1] else '',
                            n_records=len(g), n_samples=int(g['sample'].nunique()), cohorts=','.join(sorted(g.cohort.unique())),
                            n_verkko=int((g.assembler == 'verkko').sum()), known_structure=(arm, st) in known,
                            nearest_known_allele=near[0][1] if near else '', nearest_known_max_switch_delta=int(near[0][0]) if near else np.nan,
                            contains_none_block='none' in st.split('|'), candidate_allele=g['sample'].nunique() >= 2,
                            example_record=sorted(g.end_id)[0]))
    N = pd.DataFrame(out).sort_values(['candidate_allele', 'n_samples', 'arm', 'core_structure'], ascending=[False, False, True, True])
    return N

def cohort_summary(X, R, N):
    P = X[X.primary]
    S = dict(n_records_after_dedup=int(len(X)), n_primary=int(len(P)), n_primary_samples=int(P['sample'].nunique()),
             primary_by_cohort={c: dict(records=int(len(g)), samples=int(g['sample'].nunique()), frac_with_telomere=round(float(X[X.cohort == c].has_telomere.mean()), 4))
                                for c, g in P.groupby('cohort')},
             assignment_by_cohort={c: g.sa_assignment.value_counts().to_dict() for c, g in P.groupby('cohort')},
             frac_assigned_by_cohort={c: round(float((g.sa_assignment == 'full').mean()), 4) for c, g in P.groupby('cohort')})
    A = R[(R.kind == 'structural_allele') & (R.hprc_reference == 'ALL')]
    S['n_sa_tests'] = int(len(R[R.kind == 'structural_allele'])); S['n_form_tests'] = int(len(R[R.kind == 'form']))
    S['n_sa_tests_q_lt_0.05'] = int((R[R.kind == 'structural_allele'].bh_q < 0.05).sum()); S['n_form_tests_q_lt_0.05'] = int((R[R.kind == 'form'].bh_q < 0.05).sum())
    # replication of common tier-A alleles (R2 freq >= 5%)
    com = A[(A.tier == 'A') & (A.hprc_freq >= 0.05) & (A.group == 'ALL_INDEPENDENT')]
    S['replication_tierA_common'] = dict(n_alleles=int(len(com)), n_observed=int((com.n > 0).sum()), n_observed_ge_3_samples=int((com.n_samples >= 3).sum()))
    for c in T.COHORTS_INDEP:
        cc = A[(A.tier == 'A') & (A.hprc_freq >= 0.05) & (A.group == c)]
        S['replication_tierA_common'][c] = dict(n_alleles=int(len(cc)), n_observed=int((cc.n > 0).sum()))
    ab = A[A.group == 'ALL_INDEPENDENT']
    S['tierAB_observed_by_nativity_class'] = {k: dict(n=int(len(g)), observed=int((g.n > 0).sum()), observed_ge_3_samples=int((g.n_samples >= 3).sum()))
                                              for k, g in ab.groupby(ab.nativity_class.fillna(''))}
    S['tierAB_observed_in_independent'] = {t: dict(n=int((ab.tier == t).sum()), observed=int((ab[ab.tier == t].n > 0).sum()),
                                               observed_ge_3_samples=int((ab[ab.tier == t].n_samples >= 3).sum())) for t in ('A', 'B')}
    i = ab[ab.denominator > 0]
    nn = i[~i.structure.eq('home')]
    S['freq_correlation_non_native'] = dict(n_alleles=int(len(nn)), pearson_r=round(float(np.corrcoef(nn.freq, nn.hprc_freq)[0, 1]), 4),
                                            spearman_rho=round(float(nn.freq.rank().corr(nn.hprc_freq.rank())), 4))
    S['freq_correlation_independent_vs_hprc'] = dict(n_alleles=int(len(i)), pearson_r=round(float(np.corrcoef(i.freq, i.hprc_freq)[0, 1]), 4) if len(i) > 2 else None,
                                                     spearman_rho=round(float(i.freq.rank().corr(i.hprc_freq.rank())), 4) if len(i) > 2 else None)
    key = {}
    for aid in KEY_ALLELES:
        rr = R[(R.item == aid)]
        if len(rr): key[aid] = rr[['group', 'hprc_reference', 'n', 'denominator', 'freq', 'ci95_lo', 'ci95_hi', 'hprc_freq', 'bh_q']].to_dict('records')
    S['key_alleles'] = key
    S['novel'] = dict(n_unassigned_primary=int((P.sa_assignment == 'unassigned').sum()),
                      frac_unassigned_primary=round(float((P.sa_assignment == 'unassigned').mean()), 4),
                      n_clusters=int(len(N)), n_candidate_alleles=int(N.candidate_allele.sum()),
                      n_candidate_new_structure=int((N.candidate_allele & ~N.known_structure).sum()),
                      n_candidate_known_structure_new_switches=int((N.candidate_allele & N.known_structure).sum()),
                      n_candidate_ge_3_samples=int((N.n_samples >= 3).sum()), n_candidate_without_none_block=int((N.candidate_allele & ~N.contains_none_block).sum()))
    return S

if __name__ == '__main__':
    main()
