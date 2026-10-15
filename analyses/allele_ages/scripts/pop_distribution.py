#!/usr/bin/env python3
"""Population distribution of foreign subtelomere forms across 1000 Genomes superpopulations.

Inputs  ($ATLAS_INPUT_DIR = analysis home): results/subtel_foreign_tips.tsv (11,998 chromosome-named Release 2 tips),
        data/metadata/samples.tsv (repo; built by build_samples_tsv.py).
Definitions
  haplotype = one assembly (asm). Denominator for a form on own arm A = haplotypes with a named tip on A (any t20_state),
  same definition as results/subtel_foreign_forms.tsv. Carrier = haplotype whose tip on A has t20_state == foreign and that form.
  Frequencies with Wilson 95% CIs per superpopulation.
Tests (forms with >= MIN_CARRIERS carrier haplotypes)
  * chi-square statistic on the 2 x 5 (carrier/non-carrier x superpopulation) table; p-value by permutation of superpopulation
    labels among SAMPLES (both haplotypes of a sample move together; N_PERM permutations, fixed seed), unstratified and
    stratified by HiFi input chemistry batch (labels permuted only within batch);
  * Fisher exact test AFR vs non-AFR (two-sided) with odds ratio;
  * Fst: haploid ANOVA estimator (Weir 1996) across 5 superpopulations; Hudson Fst (Bhatia et al. 2013) AFR vs non-AFR.
  Multiple testing: Benjamini-Hochberg across tested forms for each p-value family.
Junction clusters: for each base-resolved junction cluster (same form, |delta| <= 1 kb), carriers by superpopulation and whether
  it is shared between AFR and non-AFR samples.
Outputs (analyses/allele_ages/results/): form_freq_by_superpop.tsv (long), form_popdiff_tests.tsv, junction_cluster_superpop.tsv,
        haplotypes_by_superpop.tsv, pop_summary.json
"""
import os, json, numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.proportion import proportion_confint

IN = os.environ.get('ATLAS_INPUT_DIR')
if not IN:
    raise SystemExit('set ATLAS_INPUT_DIR to the analysis home (see analyses/allele_ages/README.md)')
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
OUT = os.path.join(REPO, 'analyses', 'allele_ages', 'results')
os.makedirs(OUT, exist_ok=True)
MIN_CARRIERS, N_PERM, SEED = 10, 10000, 20261015
SP = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

T = pd.read_csv(os.path.join(IN, 'results', 'subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
S = pd.read_csv(os.path.join(REPO, 'data', 'metadata', 'samples.tsv'), sep='\t')
T = T.merge(S[['sample', 'superpopulation', 'population']], on='sample', how='left')
ADMIXED_AFR = {'ACB', 'ASW'}  # African-descent populations with recent European admixture
T['afr_cont'] = (T.superpopulation == 'AFR') & ~T.population.isin(ADMIXED_AFR)
T['ooa'] = T.superpopulation.isin(['EUR', 'EAS', 'SAS'])
assert T.superpopulation.isin(SP).all(), 'unlabelled samples'
# chemistry batch per sample (most common hifi_input among its tips)
batch = T.groupby('sample').hifi_input.agg(lambda s: s.fillna('unknown').mode().iloc[0])
samp_sp = T.drop_duplicates('sample').set_index('sample').superpopulation

H = T.drop_duplicates('asm').groupby('superpopulation').agg(haplotypes=('asm', 'size'), samples=('sample', 'nunique'))
H['named_tips'] = T.groupby('superpopulation').size()
H['foreign_tips'] = T[T.t20_state == 'foreign'].groupby('superpopulation').size()
H['foreign_frac_of_named_tips'] = (H.foreign_tips / H.named_tips).round(4)
H.to_csv(os.path.join(OUT, 'haplotypes_by_superpop.tsv'), sep='\t')


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(k, n, alpha=0.05, method='wilson')
    return round(lo, 4), round(hi, 4)


def chi2stat(tab):
    # tab: 2 x K counts; returns Pearson chi-square (columns with zero total dropped)
    tab = tab[:, tab.sum(0) > 0]
    e = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.nansum((tab - e) ** 2 / e)


def fst_anova(k, n):
    # haploid ANOVA estimator of theta (Weir 1996, ch.5) for a biallelic indicator; k carriers, n haplotypes per pop
    k, n = np.asarray(k, float), np.asarray(n, float)
    m = n > 0; k, n = k[m], n[m]
    r, N = len(n), n.sum()
    if r < 2 or k.sum() in (0, N):
        return np.nan
    p = k / n; pbar = k.sum() / N
    msp = (n * (p - pbar) ** 2).sum() / (r - 1)
    msg = (k * (1 - p)).sum() / (N - r)  # sum n p(1-p) / (N-r)
    nc = (N - (n ** 2).sum() / N) / (r - 1)
    return (msp - msg) / (msp + (nc - 1) * msg)


def hudson(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num / den if den > 0 else np.nan


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); m = ~np.isnan(p)
    pv = p[m]; o = np.argsort(pv); n = len(pv)
    adj = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(adj, 1); q[m] = out
    return q


rng = np.random.default_rng(SEED)
F = T[T.t20_state == 'foreign']
long_rows, test_rows = [], []
for form, g in F.groupby('form'):
    arm = g.own.iloc[0]
    A = T[T.own == arm].drop_duplicates('asm')       # assessed haplotypes on that arm
    car = set(g.asm)
    A = A.assign(carrier=A.asm.isin(car).astype(int))
    for sp in SP:
        a = A[A.superpopulation == sp]; k, n = int(a.carrier.sum()), len(a)
        lo, hi = wilson(k, n)
        long_rows.append(dict(form=form, own_arm=arm, foreign_arm=g.foreign_arm.iloc[0], superpopulation=sp, carriers=k,
                              haplotypes_assessed=n, freq=round(k / n, 4) if n else np.nan, ci95_lo=lo, ci95_hi=hi))
    k_all, n_all = int(A.carrier.sum()), len(A)
    long_rows.append(dict(form=form, own_arm=arm, foreign_arm=g.foreign_arm.iloc[0], superpopulation='ALL', carriers=k_all,
                          haplotypes_assessed=n_all, freq=round(k_all / n_all, 4), ci95_lo=wilson(k_all, n_all)[0], ci95_hi=wilson(k_all, n_all)[1]))
    if k_all < MIN_CARRIERS:
        continue
    # per-sample carrier counts and assessed haplotype counts
    ps = A.groupby('sample').agg(k=('carrier', 'sum'), n=('asm', 'size'))
    lab = samp_sp.reindex(ps.index).values
    bat = batch.reindex(ps.index).values
    sp_idx = np.array([SP.index(x) for x in lab])
    kk, nn = ps.k.values, ps.n.values

    def table(idx):
        car_ = np.bincount(idx, weights=kk, minlength=5); tot = np.bincount(idx, weights=nn, minlength=5)
        return np.vstack([car_, tot - car_])
    obs_tab = table(sp_idx); obs = chi2stat(obs_tab)
    # unstratified permutation
    perm = np.array([chi2stat(table(rng.permutation(sp_idx))) for _ in range(N_PERM)])
    p_perm = (1 + (perm >= obs - 1e-9).sum()) / (N_PERM + 1)
    # stratified by chemistry batch
    pst = np.empty(N_PERM)
    groups = [np.where(bat == b)[0] for b in np.unique(bat)]
    for i in range(N_PERM):
        idx = sp_idx.copy()
        for gi in groups:
            idx[gi] = sp_idx[rng.permutation(gi)]
        pst[i] = chi2stat(table(idx))
    p_strat = (1 + (pst >= obs - 1e-9).sum()) / (N_PERM + 1)
    p_asym = chi2_contingency(obs_tab[:, obs_tab.sum(0) > 0])[1]
    kA, nA = obs_tab[0, 0], obs_tab[:, 0].sum(); kN, nN = obs_tab[0, 1:].sum(), obs_tab[:, 1:].sum()
    orr, pf = fisher_exact([[kA, nA - kA], [kN, nN - kN]])
    fr = obs_tab[0] / obs_tab.sum(0)
    test_rows.append(dict(form=form, own_arm=arm, foreign_arm=g.foreign_arm.iloc[0], carriers=k_all, haplotypes_assessed=n_all,
                          freq_all=round(k_all / n_all, 4), **{f'freq_{s}': round(fr[i], 4) for i, s in enumerate(SP)},
                          max_min_ratio=round(fr.max() / fr.min(), 3) if fr.min() > 0 else np.inf,
                          chi2=round(obs, 3), p_chi2_asymptotic=p_asym, p_perm=p_perm, p_perm_batch_stratified=p_strat,
                          fst_anova_5pop=round(fst_anova(obs_tab[0], obs_tab.sum(0)), 4),
                          fst_hudson_AFR_nonAFR=round(hudson(kA, nA, kN, nN), 4),
                          freq_nonAFR=round(kN / nN, 4), OR_AFR_vs_nonAFR=round(orr, 3), p_fisher_AFR_nonAFR=pf,
                          n_samples_with_carrier=int((ps.k > 0).sum()),
                          superpops_with_carrier=int((obs_tab[0] > 0).sum())))
L = pd.DataFrame(long_rows)
L.to_csv(os.path.join(OUT, 'form_freq_by_superpop.tsv'), sep='\t', index=False)
R = pd.DataFrame(test_rows).sort_values('p_perm')
for c in ['p_perm', 'p_perm_batch_stratified', 'p_fisher_AFR_nonAFR']:
    R['q_' + c[2:]] = bh(R[c].values)
R.to_csv(os.path.join(OUT, 'form_popdiff_tests.tsv'), sep='\t', index=False)

# junction clusters vs superpopulation
J = F[F.junction_cluster.notna()]
jr = []
for jc, g in J.groupby('junction_cluster'):
    s = g.drop_duplicates('sample')
    c = s.superpopulation.value_counts()
    jr.append(dict(junction_cluster=jc, form=g.form.iloc[0], junction_median=float(g.junction.median()),
                   junction_range=float(g.junction.max() - g.junction.min()), tips=len(g), samples=len(s),
                   **{f'samples_{x}': int(c.get(x, 0)) for x in SP},
                   superpops=int((c > 0).sum()), samples_AFR_continental=int(s.afr_cont.sum()), samples_EUR_EAS_SAS=int(s.ooa.sum()),
                   afr=int(c.get('AFR', 0)) > 0, nonafr=int(c.drop('AFR', errors='ignore').sum()) > 0))
JR = pd.DataFrame(jr)
JR['shared_AFR_nonAFR'] = JR.afr & JR.nonafr
JR['shared_AFRcont_EUR_EAS_SAS'] = (JR.samples_AFR_continental > 0) & (JR.samples_EUR_EAS_SAS > 0)
JR = JR.drop(columns=['afr', 'nonafr']).sort_values('samples', ascending=False)
JR.to_csv(os.path.join(OUT, 'junction_cluster_superpop.tsv'), sep='\t', index=False)
# expected sharing under random labelling: probability a cluster with s samples drawn from the 230 samples contains both
nA_s = int((samp_sp == 'AFR').sum()); nT = len(samp_sp)
from scipy.special import comb
def p_both(s):
    return 1 - comb(nT - nA_s, s) / comb(nT, s) - comb(nA_s, s) / comb(nT, s)
JR2 = JR[JR.samples >= 2]
ss = T.drop_duplicates('sample'); nC = int(ss.afr_cont.sum()); nO = int(ss.ooa.sum())
def p_both2(s):
    # P(a random s-sample set contains >=1 continental-African and >=1 EUR/EAS/SAS sample)
    other = nT - nC - nO
    return 1 - comb(nT - nC, s) / comb(nT, s) - comb(nT - nO, s) / comb(nT, s) + comb(other, s) / comb(nT, s)
summ = dict(
    haplotypes=int(T.asm.nunique()), samples=int(T['sample'].nunique()),
    haplotypes_by_superpop=H.haplotypes.to_dict(), samples_by_superpop=H.samples.to_dict(),
    foreign_frac_by_superpop=H.foreign_frac_of_named_tips.to_dict(),
    forms_tested=len(R), min_carriers=MIN_CARRIERS, n_perm=N_PERM, seed=SEED,
    forms_q05_perm=int((R.q_perm < 0.05).sum()), forms_q05_perm_batch=int((R.q_perm_batch_stratified < 0.05).sum()),
    forms_q05_fisher_AFR=int((R.q_fisher_AFR_nonAFR < 0.05).sum()),
    forms_in_all_5_superpops=int((R.superpops_with_carrier == 5).sum()),
    median_fst_anova=float(R.fst_anova_5pop.median()), max_fst_anova=float(R.fst_anova_5pop.max()),
    forms_freq_AFR_higher_q05=int(((R.q_fisher_AFR_nonAFR < 0.05) & (R.OR_AFR_vs_nonAFR > 1)).sum()),
    forms_freq_AFR_lower_q05=int(((R.q_fisher_AFR_nonAFR < 0.05) & (R.OR_AFR_vs_nonAFR < 1)).sum()),
    junction_clusters_ge2_samples=len(JR2), junction_clusters_ge2_shared_AFR_nonAFR=int(JR2.shared_AFR_nonAFR.sum()),
    junction_clusters_ge2_expected_shared_random=round(float(sum(p_both(s) for s in JR2.samples)), 2),
    samples_AFR_continental=nC, samples_EUR_EAS_SAS=nO,
    junction_clusters_ge2_shared_AFRcont_OOA=int(JR2.shared_AFRcont_EUR_EAS_SAS.sum()),
    junction_clusters_ge2_expected_shared_AFRcont_OOA_random=round(float(sum(p_both2(s) for s in JR2.samples)), 2),
    junction_clusters_ge5_shared_AFRcont_OOA=int(JR[JR.samples >= 5].shared_AFRcont_EUR_EAS_SAS.sum()),
    junction_clusters_ge10_samples=int((JR.samples >= 10).sum()),
    junction_clusters_ge10_shared_AFRcont_OOA=int(JR[JR.samples >= 10].shared_AFRcont_EUR_EAS_SAS.sum()),
    junction_cluster_samples_ge2_total=int(JR2.samples.sum()),
    junction_clusters_ge5_samples=int((JR.samples >= 5).sum()),
    junction_clusters_ge5_shared_AFR_nonAFR=int(JR[JR.samples >= 5].shared_AFR_nonAFR.sum()),
    junction_clusters_ge5_in_all5=int(JR[JR.samples >= 5].superpops.eq(5).sum()),
)
json.dump(summ, open(os.path.join(OUT, 'pop_summary.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))
