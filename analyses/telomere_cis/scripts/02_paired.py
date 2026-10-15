"""Step 2 - paired test of H1 (telomere length set in cis by the subtelomere allele) in HPRC Release 2 assemblies.

Unit: an individual heterozygous at an arm, i.e. both haplotype tips of the arm are named and one carries a foreign terminal
block (t20_state 'foreign', form 'R>D': recipient arm R carrying a block best matching CHM13 arm D) and the other the own
form (t20_state 'own'). 'tie'/'none' tips are excluded. Both haplotypes of a pair come from the same sample, the same
assembler run and the same HiFi input, so those factors cancel in the within-pair difference d = TL(foreign) - TL(own).
Tests
  T1 pooled Wilcoxon signed-rank on d (+ exact-style sign-flip permutation, 20,000 draws) and per common form (>= 10 pairs, BH);
  T2 within-arm label-shuffle control: foreign/own labels permuted among all foreign/own tips of the same arm (all samples),
     statistic recomputed on the pairs that are heterozygous after shuffling (5,000 permutations);
  T3 unpaired model over all foreign/own tips: rel_tl ~ arm + foreign (+ haplotype index), within-arm permutation p;
  T4 donor-likeness: regression of d on x = base(D) - base(R), where base = mean relative telomere length of own-form tips of
     the arm, leave-sample-out (the pair's own sample is removed from both baselines). Slope 1 = foreign-form end takes the
     donor arm's length; slope 0 = length follows the recipient chromosome. Repeated with the external Karimian et al. 2024
     arm values (x_K) which share no data with the assemblies.
  Sensitivity: hifiasm only, trio-phased only, Hi-C-phased only, both junk_len == 0, foreign tip with junction shared by >= 2
     samples (validated recurrent allele), haplotype-index check.
Power: minimum detectable mean difference (80 %, two-sided alpha 0.05) = 2.80 * SD(d) / sqrt(n).
Inputs : results/tip_telomere_lengths.tsv, results/arm_baseline.tsv (step 1), data/karimian2024_fig4e.tsv
Outputs: results/het_pairs.tsv, results/paired_tests.tsv, results/paired_by_form.tsv, results/paired_stats.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from tl_common import RES, DAT, SEED

rng = np.random.default_rng(SEED)
# allele definition: 't20' = published state (tip[s0:20 kb], telomere-length dependent window);
# 'S20' (primary), 'S10', 'S5' = fixed subtelomere-anchored windows from scripts/01b_fixed_window_states.py
TAG = sys.argv[1] if len(sys.argv) > 1 else 'S20'
T = pd.read_csv(os.path.join(RES, 'tip_telomere_lengths.tsv'), sep='\t', low_memory=False)
FW = pd.read_csv(os.path.join(RES, 'fixed_window_states.tsv'), sep='\t')
T = T.merge(FW[['q', 'S20_n_near', 'T20R_n_near'] + ([f'{TAG}_state', f'{TAG}_best_end'] if TAG != 't20' else [])], on='q', how='left')
if TAG != 't20':
    T['t20_state'] = T[f'{TAG}_state']; T['foreign_arm'] = np.where(T.t20_state == 'foreign', T[f'{TAG}_best_end'], np.nan)
    T['form'] = np.where(T.t20_state == 'foreign', T.arm + '>' + T.foreign_arm.astype(str), np.nan)
T['n_near'] = T.S20_n_near if TAG != 't20' else T.T20R_n_near
SUF = '' if TAG == 'S20' else f'_{TAG}'
T = T[T.t20_state.isin(['foreign', 'own'])].copy()
T['isf'] = (T.t20_state == 'foreign').astype(int)
T['phasing'] = np.where(T.hap.isin(['mat', 'pat']), 'trio', 'hic')
# ---- leave-sample-out own-form arm baselines (relative TL)
own = T[T.isf == 0]
arm_sum = own.groupby('arm').rel_tl.sum(); arm_n = own.groupby('arm').rel_tl.size()
so = own.groupby(['sample', 'arm']).rel_tl.agg(['sum', 'size'])
def base_lso(sample, arm):
    if arm not in arm_n: return np.nan
    s, n = arm_sum[arm], arm_n[arm]
    if (sample, arm) in so.index: s -= so.loc[(sample, arm), 'sum']; n -= so.loc[(sample, arm), 'size']
    return s / n if n >= 10 else np.nan
K = pd.read_csv(os.path.join(DAT, 'karimian2024_fig4e.tsv'), sep='\t').set_index('arm').rel_mean_kb_digitized * 1000
# ---- pairs
g = T.groupby(['sample', 'arm'])
rows = []
for (s, a), d in g:
    if len(d) != 2 or d.hapn.nunique() != 2: continue
    rows.append((s, a, d))
pairs = []
for s, a, d in rows:
    if d.isf.sum() != 1: continue
    f = d[d.isf == 1].iloc[0]; o = d[d.isf == 0].iloc[0]
    D = f.foreign_arm
    pairs.append(dict(sample=s, arm=a, form=f.form, donor=D, assembler=f.assembler, hifi_input=f.hifi_input, phasing=f.phasing,
        q_foreign=f.q, q_own=o.q, tl_foreign=f.telo_len, tl_own=o.telo_len, d=f.telo_len - o.telo_len,
        foreign_hapn=int(f.hapn), junk_foreign=f.junk_len, junk_own=o.junk_len,
        foreign_junction_samples=f.junction_cluster_samples, base_R=base_lso(s, a), base_D=base_lso(s, D),
        karimian_R=K.get(a, np.nan), karimian_D=K.get(D, np.nan)))
P = pd.DataFrame(pairs)
P['x'] = P.base_D - P.base_R; P['x_karimian'] = P.karimian_D - P.karimian_R
P.to_csv(os.path.join(RES, f'het_pairs{SUF}.tsv'), sep='\t', index=False, float_format='%.5g')

def signflip_p(d, n=20000):
    d = np.asarray(d, float); obs = abs(d.mean())
    sg = rng.choice([-1, 1], size=(n, len(d)))
    return float(((np.abs((sg * d).mean(1)) >= obs).sum() + 1) / (n + 1))
def boot_ci(d, n=5000):
    d = np.asarray(d, float); m = rng.choice(d, size=(n, len(d))).mean(1); return np.percentile(m, [2.5, 97.5])
def test(d, label):
    d = np.asarray(d, float); n = len(d)
    if n < 5: return dict(test=label, n=n)
    w = stats.wilcoxon(d); ci = boot_ci(d)
    return dict(test=label, n=n, mean_d=d.mean(), median_d=np.median(d), ci95_lo=ci[0], ci95_hi=ci[1], sd_d=d.std(ddof=1),
                frac_foreign_longer=(d > 0).mean(), wilcoxon_p=w.pvalue, signflip_p=signflip_p(d), mde80=2.80 * d.std(ddof=1) / np.sqrt(n))
R = [test(P.d, 'pooled_all_het_pairs')]
sub = {'hifiasm_only': P.assembler == 'hifiasm', 'trio_phased': P.phasing == 'trio', 'hic_phased': P.phasing == 'hic',
       'no_terminal_junk_both': (P.junk_foreign == 0) & (P.junk_own == 0),
       'validated_recurrent_junction': P.foreign_junction_samples.fillna(0) >= 2,
       'foreign_on_hap1': P.foreign_hapn == 1, 'foreign_on_hap2': P.foreign_hapn == 2}
for k, m in sub.items(): R.append(test(P.d[m], k))
for k, d in P.groupby('hifi_input'): R.append(test(d.d, f'hifi_input={k}'))
# ---- T2 within-arm label shuffle
tips = T[['sample', 'arm', 'hapn', 'telo_len', 'isf']].reset_index(drop=True)
pair_idx = []   # index pairs (i,j) of the two tips of each two-haplotype arm with foreign/own states
for (s, a), d in tips.groupby(['sample', 'arm']):
    if len(d) == 2 and d.hapn.nunique() == 2: pair_idx.append(d.index.values)
pair_idx = np.array(pair_idx); tl = tips.telo_len.values
arm_groups = [d.index.values for _, d in tips.groupby('arm')]
def stat(lab):
    a, b = lab[pair_idx[:, 0]], lab[pair_idx[:, 1]]; het = a != b
    dd = np.where(a[het] == 1, tl[pair_idx[het, 0]] - tl[pair_idx[het, 1]], tl[pair_idx[het, 1]] - tl[pair_idx[het, 0]])
    return dd.mean(), het.sum()
obs, nhet = stat(tips.isf.values)
null = np.empty(5000); lab = tips.isf.values.copy()
for i in range(len(null)):
    for ix in arm_groups: lab[ix] = rng.permutation(tips.isf.values[ix])
    null[i] = stat(lab)[0]
t2 = dict(observed_mean_d=float(obs), n_het_pairs=int(nhet), null_mean=float(null.mean()), null_sd=float(null.std()),
          perm_p_two_sided=float(((np.abs(null - null.mean()) >= abs(obs - null.mean())).sum() + 1) / (len(null) + 1)))
# ---- T3 unpaired: rel_tl ~ arm + foreign + hapn, within-arm permutation
def ols_coef(y, X): return np.linalg.lstsq(X, y, rcond=None)[0]
armd = pd.get_dummies(T.arm, dtype=float).values
X = np.column_stack([armd, T.isf.values, (T.hapn.values == 2).astype(float)])
y = T.rel_tl.values; b = ols_coef(y, X); bf = b[-2]; bh = b[-1]
resid = y - X @ b; se = np.sqrt(np.diag(resid @ resid / (len(y) - X.shape[1]) * np.linalg.pinv(X.T @ X)))
ix_groups = [np.where(T.arm.values == a)[0] for a in T.arm.unique()]
nullb = np.empty(2000); lab = T.isf.values.copy()
for i in range(len(nullb)):
    for ix in ix_groups: lab[ix] = rng.permutation(T.isf.values[ix])
    X[:, -2] = lab; nullb[i] = ols_coef(y, X)[-2]
t3 = dict(n_tips=int(len(T)), n_foreign=int(T.isf.sum()), beta_foreign_bp=float(bf), se_bp=float(se[-2]), beta_hap2_bp=float(bh),
          perm_p=float(((np.abs(nullb) >= abs(bf)).sum() + 1) / (len(nullb) + 1)), mde80_bp=float(2.80 * se[-2]))
# paralogy covariate: number of other CHM13 ends within 2 % of the best alignment score (S20 window; T20R for t20)
nn = np.minimum(T.n_near.fillna(0).values, 4)
X2 = np.column_stack([armd, T.isf.values, (T.hapn.values == 2).astype(float), nn]); b2 = ols_coef(y, X2)
r2 = y - X2 @ b2; se2 = np.sqrt(np.diag(r2 @ r2 / (len(y) - X2.shape[1]) * np.linalg.pinv(X2.T @ X2)))
t3.update(beta_foreign_adj_nnear_bp=float(b2[-3]), se_adj_nnear_bp=float(se2[-3]), beta_per_near_end_bp=float(b2[-1]), se_per_near_end_bp=float(se2[-1]))
own_only = T.isf.values == 0
X3 = np.column_stack([armd[own_only], nn[own_only]]); b3 = ols_coef(y[own_only], X3); r3 = y[own_only] - X3 @ b3
se3 = np.sqrt(np.diag(r3 @ r3 / (own_only.sum() - X3.shape[1]) * np.linalg.pinv(X3.T @ X3)))
t3.update(own_tips_beta_per_near_end_bp=float(b3[-1]), own_tips_se_per_near_end_bp=float(se3[-1]))
ob = T[T.isf == 0].groupby('arm').rel_tl.agg(['size', 'mean']).rename(columns={'size': 'n_own', 'mean': 'own_mean_rel_tl'})
ob.to_csv(os.path.join(RES, f'arm_own_baseline{SUF}.tsv'), sep='\t', float_format='%.5g')
# ---- T4 donor-likeness
def reg(x, y, label, groups=None):
    """OLS slope with cluster-robust (CR1) SE by allele form (tips/pairs of one form share x) and a form-level permutation
    (x values permuted among forms, 2,000 draws)."""
    import statsmodels.api as sm
    x = np.asarray(x, float); y = np.asarray(y, float); groups = np.asarray(groups if groups is not None else np.arange(len(x)))
    m = np.isfinite(x) & np.isfinite(y); x, y, groups = x[m], y[m], groups[m]
    X = sm.add_constant(x); fit = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(groups)[0]})
    b, se = fit.params[1], fit.bse[1]; ng = len(set(groups))
    ug, gi = np.unique(groups, return_inverse=True); xg = np.array([x[gi == k][0] for k in range(len(ug))])
    nul = []
    for _ in range(2000):
        xp = rng.permutation(xg)[gi]; nul.append(np.polyfit(xp, y, 1)[0])
    df = max(ng - 1, 1)
    return dict(test=label, n=int(m.sum()), n_forms=ng, slope=b, slope_se_cluster=se, slope_ci95=[b - 1.96 * se, b + 1.96 * se],
                intercept=fit.params[0], r=float(np.corrcoef(x, y)[0, 1]), p=2 * stats.t.sf(abs(b / se), df),
                perm_p_forms=float(((np.abs(nul) >= abs(b)).sum() + 1) / 2001), p_slope_eq_1=2 * stats.t.sf(abs((b - 1) / se), df), sd_x=float(np.std(x)))
t4 = [reg(P.x.values, P.d.values, 'd ~ base(D)-base(R) [HPRC own-form, leave-sample-out]', P.form.values),
      reg(P.x_karimian.values, P.d.values, 'd ~ Karimian(D)-Karimian(R)', P.form.values)]
# form-level (one point per form, mean d) to avoid weighting by common forms
fm = P.groupby('form').agg(n=('d', 'size'), mean_d=('d', 'mean'), x=('x', 'mean'), x_k=('x_karimian', 'mean')).reset_index()
fm5 = fm[fm.n >= 5]
t4.append(reg(fm5.x.values, fm5.mean_d.values, 'form-level mean d ~ x (forms with >= 5 pairs)', fm5.form.values))
# unpaired: every foreign tip of the allele definition
Fo = T[T.isf == 1].copy()
Fo['bR'] = [base_lso(s, a) for s, a in zip(Fo['sample'], Fo.arm)]; Fo['bD'] = [base_lso(s, a) for s, a in zip(Fo['sample'], Fo.foreign_arm)]
# multiple regression rel_tl ~ base(R) + base(D): H1 (donor) predicts beta_D = 1, beta_R = 0; recipient control predicts 0 / 1.
# (A difference-on-difference form would share the estimation error of base(R) between y and x and bias the slope upwards.)
import statsmodels.api as sm
fm_ = Fo.dropna(subset=['bR', 'bD']).copy()
Xm = sm.add_constant(fm_[['bR', 'bD']].values); ym = fm_.rel_tl.values; grp = pd.factorize(fm_.form)[0]
fit = sm.OLS(ym, Xm).fit(cov_type='cluster', cov_kwds={'groups': grp}); ng = len(set(grp))
t4.append(dict(test='all foreign tips: rel_tl ~ base(R) + base(D) [cluster-robust by form]', n=len(fm_), n_forms=ng, beta_R=fit.params[1], se_R=fit.bse[1],
               beta_D=fit.params[2], se_D=fit.bse[2], p_D=2 * stats.t.sf(abs(fit.params[2] / fit.bse[2]), ng - 1),
               p_R=2 * stats.t.sf(abs(fit.params[1] / fit.bse[1]), ng - 1)))
# Karimian version (external arm values)
fk = Fo.assign(kR=Fo.arm.map(K), kD=Fo.foreign_arm.map(K)).dropna(subset=['kR', 'kD'])
fitk = sm.OLS(fk.rel_tl.values, sm.add_constant(fk[['kR', 'kD']].values)).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(fk.form)[0]})
ngk = fk.form.nunique()
t4.append(dict(test='all foreign tips: rel_tl ~ Karimian(R) + Karimian(D) [cluster-robust by form]', n=len(fk), n_forms=ngk, beta_R=fitk.params[1], se_R=fitk.bse[1],
               beta_D=fitk.params[2], se_D=fitk.bse[2], p_D=2 * stats.t.sf(abs(fitk.params[2] / fitk.bse[2]), ngk - 1),
               p_R=2 * stats.t.sf(abs(fitk.params[1] / fitk.bse[1]), ngk - 1)))
# ---- per form
B = []
for f, d in P.groupby('form'):
    if len(d) < 10: continue
    r = test(d.d, f); r.update(form=f, donor=d.donor.iloc[0], recipient=d.arm.iloc[0], x_mean=d.x.mean(), x_karimian=d.x_karimian.mean()); B.append(r)
B = pd.DataFrame(B).sort_values('n', ascending=False)
p = B.wilcoxon_p.values; o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
B.loc[B.index[o], 'bh_q'] = np.minimum(q, 1)
B.to_csv(os.path.join(RES, f'paired_by_form{SUF}.tsv'), sep='\t', index=False, float_format='%.4g')
pd.DataFrame(R).to_csv(os.path.join(RES, f'paired_tests{SUF}.tsv'), sep='\t', index=False, float_format='%.4g')
# repeatability of assembled TL between the two haplotypes of an arm (own/own pairs, within sample) - biological + technical
oo = [(d.telo_len.values) for s, a, d in rows if d.isf.sum() == 0]
oo = np.array(oo); rep = stats.spearmanr(oo[:, 0], oo[:, 1]).statistic if len(oo) else np.nan
out = dict(allele_definition=TAG, n_two_haplotype_arms=len(rows), n_het_pairs=len(P), n_forms_with_pairs=int(P.form.nunique()),
           n_forms_tested_ge10=int(len(B)), n_forms_bh_q_lt_0_05=int((B.bh_q < 0.05).sum()), T2_within_arm_shuffle=t2, T3_unpaired=t3,
           T4_donor_likeness=t4, own_own_pairs=int(len(oo)), own_own_hap_spearman=float(rep),
           mean_abs_x_bp=float(np.nanmean(np.abs(P.x))), sd_x_bp=float(np.nanstd(P.x)))
json.dump(out, open(os.path.join(RES, f'paired_stats{SUF}.json'), 'w'), indent=1, default=float)
pd.set_option('display.width', 250)
print(pd.DataFrame(R).round(4).to_string(index=False)); print(B[['form', 'n', 'mean_d', 'median_d', 'wilcoxon_p', 'bh_q', 'x_mean', 'mde80']].round(3).to_string(index=False))
print(json.dumps(out, indent=1, default=lambda v: round(float(v), 4)))
