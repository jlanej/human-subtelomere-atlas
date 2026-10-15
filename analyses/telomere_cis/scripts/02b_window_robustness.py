"""Step 2b - is the negative foreign-vs-own telomere-length difference under boundary-anchored windows an artefact?

The paired difference d = TL(foreign) - TL(own) is negative under all boundary-anchored allele definitions (S5 < S10 < S20)
and positive only under the telomere-length-dependent published window (t20). The windows closest to the telomere carry
the largest effect, so the windows could be coupled to the telomere call. Checks:
  (a) boundary QC: tail_frac = fraction of the first 500 bp after the called boundary covered by CCCTAA/TTAGGG units;
      ITS = exact run of >= 4 units within [s0+500, s0+5 kb). Graded exclusions (strict: tail < 0.2 and no ITS; premature
      boundary: tail >= 0.5; ITS), a pair-level regression of d on the foreign-own differences in tail_frac and ITS bp
      (cluster-robust by sample), and a boundary-robust telomere length TL_ext (extension from the called boundary while
      500 bp windows, step 100 bp, keep >= 30 % unit coverage; capped at 30 kb) used to recompute d;
  (b) W5_15: new allele state on the window [s0+5 kb, s0+15 kb), which skips the first 5 kb next to the telomere
      (classification rules identical to scripts/01b_fixed_window_states.py);
  (c) concordant pairs: both tips have the same state under S5, S10 and S20 (and, stricter, also W5_15);
  (d) stratified paired tests by HiFi input, assembler and phasing, with an inverse-variance combined estimate across strata.
Heterogeneity across windows: cluster bootstrap over samples (2,000 draws; pairs of all windows recomputed per draw); 95 %
intervals of the differences between window estimates. The windows share individuals, so this is used instead of a
Cochran's Q test, which assumes independence.
Inputs : results/tip_telomere_lengths.tsv, results/fixed_window_states.tsv, 100 kb windows (as in 01b), CHM13 ends.
Outputs: work/telomere_cis/w515.{fa,paf}, results/window_qc.tsv, results/window_robustness.tsv, results/window_robustness.json
"""
import os, sys, re, json, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from tl_common import inp, RES, WORK, REPO, rc, SEED

rng = np.random.default_rng(SEED); THREADS = int(os.environ.get('ATLAS_THREADS', '3')); NEAR = 0.98
REF = inp('subtel', 'foreign', 'ref', 'chm13_ends500k.fa')
T = pd.read_csv(os.path.join(RES, 'tip_telomere_lengths.tsv'), sep='\t', low_memory=False)
T['s0'] = T.junk_len + T.telo_len
FW = pd.read_csv(os.path.join(RES, 'fixed_window_states.tsv'), sep='\t')
plan = pd.read_csv(os.path.join(REPO, 'work', 'windows', 'r2_plan.tsv'), sep='\t').set_index('q')
UNIT = re.compile('CCCTAA|TTAGGG'); ITS = re.compile('(?:CCCTAA){4,}|(?:TTAGGG){4,}')
def readfa(fn):
    name = None; buf = []
    for line in open(fn):
        if line[0] == '>':
            if name: yield name, ''.join(buf)
            name = line[1:].strip().split()[0]; buf = []
        else: buf.append(line.strip())
    if name: yield name, ''.join(buf)
qfa = os.path.join(WORK, 'w515.fa'); paf = os.path.join(WORK, 'w515.paf'); qc = []
with open(qfa, 'w') as fo:
    for asm, g in T.groupby('asm'):
        src = {}
        for fn in (os.path.join(REPO, 'work', 'windows', 'fetch', f'{asm}.fa'), inp('subtel', 'foreign', 'r2_100k', f'{asm}.fa')):
            if os.path.exists(fn):
                for n, s in readfa(fn): src.setdefault(n, s)
        for r in g.itertuples():
            s = src[plan.at[r.q, 'region']].upper(); s = s if r.end == 'L' else rc(s); s0 = int(r.s0)
            tail = s[s0:s0 + 500]; its = [m for m in ITS.finditer(s, s0 + 500, s0 + 5000)]
            e = s0
            while e + 500 <= len(s) and e - s0 < 30000 and 6 * len(UNIT.findall(s, e, e + 500)) / 500 >= 0.3: e += 100
            qc.append(dict(q=r.q, boundary_tail_frac=round(6 * len(UNIT.findall(tail)) / max(1, len(tail)), 3),
                           its_runs_5kb=len(its), its_bp_5kb=sum(m.end() - m.start() for m in its), tl_ext=e - int(r.junk_len), ext_bp=e - s0))
            sub = s[s0 + 5000:s0 + 15000]
            if len(sub) >= 9000: fo.write(f'>{r.q}@W515\n{sub}\n')
subprocess.run(f'minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 -t {THREADS} {REF} {qfa} > {paf} 2> {paf}.log', shell=True, check=True)
best = collections.defaultdict(dict)
for line in open(paf):
    p = line.split('\t')
    if int(p[3]) - int(p[2]) < 0.5 * int(p[1]): continue
    AS = int(next(t[5:] for t in p[12:] if t.startswith('AS:i:')))
    if AS > best[p[0]].get(p[5], -1): best[p[0]][p[5]] = AS
Q = pd.DataFrame(qc).merge(T[['q', 'arm']], on='q')
st, be = [], []
for q, a in zip(Q.q, Q.arm):
    d = best.get(f'{q}@W515')
    if not d: st.append('none'); be.append(None); continue
    b = max(d, key=d.get); oAS = d.get(f'chr{a}', 0)
    st.append('own' if b == f'chr{a}' else ('tie' if oAS >= NEAR * d[b] else 'foreign')); be.append(b.replace('chr', ''))
Q['W515_state'] = st; Q['W515_best_end'] = be
Q['uncertain_boundary'] = (Q.boundary_tail_frac >= 0.2) | (Q.its_runs_5kb > 0)
Q.drop(columns='arm').to_csv(os.path.join(RES, 'window_qc.tsv'), sep='\t', index=False)
# ---------------- pairs
X = T.merge(FW[['q', 'S5_state', 'S10_state', 'S20_state', 'T20R_state']], on='q').merge(Q.drop(columns='arm'), on='q')
X['phasing'] = np.where(X.hap.isin(['mat', 'pat']), 'trio', 'hic')
X['conc3'] = (X.S5_state == X.S10_state) & (X.S10_state == X.S20_state)
X['conc4'] = X.conc3 & (X.W515_state == X.S20_state)
two = X.groupby(['sample', 'arm']).filter(lambda d: len(d) == 2 and d.hapn.nunique() == 2).sort_values(['sample', 'arm', 'hapn'])
A = two.iloc[0::2].reset_index(drop=True); B = two.iloc[1::2].reset_index(drop=True)
assert (A['sample'].values == B['sample'].values).all() and (A.arm.values == B.arm.values).all()
def pairs(state, mask=None, col='telo_len'):
    """d (foreign - own) for heterozygous foreign/own pairs under a state column; mask = extra per-pair filter"""
    sa, sb = A[state].values, B[state].values
    het = ((sa == 'foreign') & (sb == 'own')) | ((sa == 'own') & (sb == 'foreign'))
    if mask is not None: het &= mask
    d = np.where(sa == 'foreign', A[col].values - B[col].values, B[col].values - A[col].values)
    return het, d
def fdiff(state, col):
    sa = A[state].values; return np.where(sa == 'foreign', A[col].values - B[col].values, B[col].values - A[col].values)
def summarize(label, het, d, extra=None):
    x = d[het]; n = len(x); r = dict(test=label, n=n)
    if n >= 5:
        smp = A['sample'].values[het]; us = np.unique(smp); idx = {s: np.where(smp == s)[0] for s in us}
        bs = [x[np.concatenate([idx[s] for s in rng.choice(us, len(us))])].mean() for _ in range(2000)]
        r.update(mean_d=x.mean(), median_d=np.median(x), ci95_lo=np.percentile(bs, 2.5), ci95_hi=np.percentile(bs, 97.5),
                 se_boot=np.std(bs), wilcoxon_p=stats.wilcoxon(x).pvalue, frac_foreign_longer=(x > 0).mean(), mde80=2.8 * x.std(ddof=1) / np.sqrt(n))
    if extra: r.update(extra)
    return r
clean = ~(A.uncertain_boundary.values | B.uncertain_boundary.values)
nopre = ~((A.boundary_tail_frac.values >= 0.5) | (B.boundary_tail_frac.values >= 0.5))
noits = ~((A.its_runs_5kb.values > 0) | (B.its_runs_5kb.values > 0))
c3 = A.conc3.values & B.conc3.values; c4 = A.conc4.values & B.conc4.values
rows = []
for w in ['T20R', 'S20', 'S10', 'S5', 'W515']:
    h, d = pairs(f'{w}_state'); rows.append(summarize(f'{w}', h, d, dict(window=w, subset='all')))
h, d = pairs('T20R_state', col='tl_ext'); rows.append(summarize('T20R | boundary-robust TL_ext', h, d, dict(window='T20R', subset='tl_ext')))
for w in ['S20', 'S10', 'S5', 'W515']:
    h, d = pairs(f'{w}_state', clean); rows.append(summarize(f'{w} | boundary-clean (tail<0.2, no ITS)', h, d, dict(window=w, subset='boundary_clean')))
    h, d = pairs(f'{w}_state', nopre); rows.append(summarize(f'{w} | no premature boundary (tail<0.5)', h, d, dict(window=w, subset='no_premature')))
    h, d = pairs(f'{w}_state', noits); rows.append(summarize(f'{w} | no ITS in first 5 kb', h, d, dict(window=w, subset='no_its')))
    h, d = pairs(f'{w}_state', col='tl_ext'); rows.append(summarize(f'{w} | boundary-robust TL_ext', h, d, dict(window=w, subset='tl_ext')))
h, d = pairs('S20_state', c3); rows.append(summarize('concordant S5=S10=S20', h, d, dict(window='S5=S10=S20', subset='concordant3')))
h, d = pairs('S20_state', c4); rows.append(summarize('concordant S5=S10=S20=W5_15', h, d, dict(window='S5=S10=S20=W5_15', subset='concordant4')))
h, d = pairs('S20_state', c3 & nopre); rows.append(summarize('concordant S5=S10=S20 | no premature boundary', h, d, dict(window='S5=S10=S20', subset='concordant3_nopre')))
h, d = pairs('S20_state', c3, col='tl_ext'); rows.append(summarize('concordant S5=S10=S20 | TL_ext', h, d, dict(window='S5=S10=S20', subset='concordant3_tl_ext')))
# pair-level adjustment for boundary features (cluster-robust by sample)
import statsmodels.api as sm
adj = {}
for w in ['S20', 'S10', 'S5', 'W515']:
    h, d = pairs(f'{w}_state')
    Z = np.column_stack([np.ones(h.sum()), fdiff(f'{w}_state', 'boundary_tail_frac')[h], fdiff(f'{w}_state', 'its_bp_5kb')[h] / 100])
    f = sm.OLS(d[h], Z).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(A['sample'].values[h])[0]})
    adj[w] = dict(n=int(h.sum()), adjusted_mean_d=f.params[0], se=f.bse[0], p=f.pvalues[0], beta_dtail_per_unit=f.params[1], se_dtail=f.bse[1],
                  beta_dits_per_100bp=f.params[2], se_dits=f.bse[2], mean_dtail=float(Z[:, 1].mean()), mean_dits_bp=float(Z[:, 2].mean() * 100))
# (d) strata for the concordant set and for S20
strata_out = {}
for setname, m0 in [('concordant3', c3), ('S20', None)]:
    for col in ['hifi_input', 'assembler', 'phasing']:
        ests = []
        for k in sorted(A[col].unique()):
            m = (A[col].values == k) & (m0 if m0 is not None else True)
            h, d = pairs('S20_state', m); r = summarize(f'{setname} | {col}={k}', h, d, dict(window='S20', subset=f'{setname}_{col}'))
            rows.append(r)
            if r['n'] >= 30: ests.append((r['mean_d'], r['se_boot']))
        e = np.array(ests); w_ = 1 / e[:, 1] ** 2; comb = (w_ * e[:, 0]).sum() / w_.sum(); se = 1 / np.sqrt(w_.sum())
        Qh = (w_ * (e[:, 0] - comb) ** 2).sum(); dfh = len(e) - 1
        strata_out[f'{setname}|{col}'] = dict(combined_mean_d=comb, se=se, ci95=[comb - 1.96 * se, comb + 1.96 * se], p=2 * stats.norm.sf(abs(comb / se)),
                                             n_strata=len(e), cochran_Q=Qh, Q_p=float(stats.chi2.sf(Qh, dfh)) if dfh > 0 else None,
                                             I2=max(0.0, (Qh - dfh) / Qh) if Qh > 0 else 0.0)
R = pd.DataFrame(rows); R.to_csv(os.path.join(RES, 'window_robustness.tsv'), sep='\t', index=False, float_format='%.4g')
# ---------------- heterogeneity across windows (cluster bootstrap over samples)
wins = ['S20', 'S10', 'S5', 'W515']; H = {w: pairs(f'{w}_state') for w in wins}
smp = A['sample'].values; us = np.unique(smp); idx = {s: np.where(smp == s)[0] for s in us}
obs = {w: H[w][1][H[w][0]].mean() for w in wins}; boot = {w: [] for w in wins}
for _ in range(2000):
    ii = np.concatenate([idx[s] for s in rng.choice(us, len(us))])
    for w in wins:
        h, d = H[w][0][ii], H[w][1][ii]; boot[w].append(d[h].mean())
het = {}
for i, a in enumerate(wins):
    for b in wins[i + 1:]:
        diff = np.array(boot[a]) - np.array(boot[b])
        het[f'{a}-{b}'] = dict(diff=obs[a] - obs[b], ci95=list(np.percentile(diff, [2.5, 97.5])), p_boot=float(2 * min((diff <= 0).mean(), (diff >= 0).mean())))
# overlap of foreign calls between windows (tip level)
ov = {f'{a}&{b}': int(((X[f'{a}_state'] == 'foreign') & (X[f'{b}_state'] == 'foreign')).sum()) for a in ['S5', 'S10', 'W515'] for b in ['S20']}
# does boundary QC / ITS associate with telomere length and with S5 foreign state? (arm-adjusted, relative TL)
import statsmodels.formula.api as smf
mq = smf.ols('rel_tl ~ C(arm) + uncertain_boundary', data=X).fit()
fs = pd.crosstab(X.uncertain_boundary, X.S5_state == 'foreign', normalize='index')
out = dict(n_tips=len(X), n_uncertain_boundary=int(X.uncertain_boundary.sum()), n_its_5kb=int((X.its_runs_5kb > 0).sum()),
           n_tail_ge_0_2=int((X.boundary_tail_frac >= 0.2).sum()),
           uncertain_boundary_rel_tl_effect_bp=dict(beta=float(mq.params['uncertain_boundary[T.True]']), se=float(mq.bse['uncertain_boundary[T.True]'])),
           frac_S5_foreign_by_uncertain_boundary={str(k): float(v) for k, v in fs[True].items()},
           W515_states=X.W515_state.value_counts().to_dict(), foreign_overlap_with_S20=ov,
           n_conc3_tips=int(X.conc3.sum()), n_conc4_tips=int(X.conc4.sum()),
           ext_bp_median=float(X.ext_bp.median()), frac_ext_gt_0=float((X.ext_bp > 0).mean()),
           rel_tl_vs_tl_ext_spearman=float(stats.spearmanr(X.telo_len, X.tl_ext).statistic),
           boundary_adjusted=adj, window_heterogeneity_bootstrap=het, stratified=strata_out)
# ---------------- TL_ext re-analysis of the arm baseline, donor-likeness and per-form tests (S20 states)
from tl_common import DAT
X['rel_tl_ext'] = X.tl_ext - X.groupby('sample').tl_ext.transform('mean')
K = pd.read_csv(os.path.join(DAT, 'karimian2024_fig4e.tsv'), sep='\t').set_index('arm').rel_mean_kb_digitized * 1000
am = X.groupby('arm').rel_tl_ext.mean(); jk = am.to_frame('m').join(K.rename('k'), how='inner')
out['karimian_vs_hprc_rel_tl_ext'] = dict(n=len(jk), rho=float(stats.spearmanr(jk.k, jk.m).statistic), p=float(stats.spearmanr(jk.k, jk.m).pvalue))
own = X[X.S20_state == 'own']; asum = own.groupby('arm').rel_tl_ext.sum(); an = own.groupby('arm').rel_tl_ext.size()
so = own.groupby(['sample', 'arm']).rel_tl_ext.agg(['sum', 'size'])
def blso(smp_, arm):
    if arm not in an: return np.nan
    s_, n_ = asum[arm], an[arm]
    if (smp_, arm) in so.index: s_ -= so.loc[(smp_, arm), 'sum']; n_ -= so.loc[(smp_, arm), 'size']
    return s_ / n_ if n_ >= 10 else np.nan
X2 = X.merge(FW[['q', 'S20_best_end']], on='q')
A2 = A.merge(FW[['q', 'S20_best_end']], on='q', how='left'); B2 = B.merge(FW[['q', 'S20_best_end']], on='q', how='left')
h, dx = pairs('S20_state', col='tl_ext')
fa = np.where(A.S20_state.values == 'foreign', True, False)
donor = np.where(fa, A2.S20_best_end.values, B2.S20_best_end.values); arm_ = A.arm.values; smp_ = A['sample'].values
xv = np.array([blso(s_, D) - blso(s_, a_) if hh else np.nan for s_, D, a_, hh in zip(smp_, donor, arm_, h)])
form = np.array([f'{a_}>{D}' for a_, D in zip(arm_, donor)])
m = h & np.isfinite(xv)
f = sm.OLS(dx[m], sm.add_constant(xv[m])).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(form[m])[0]})
out['donor_likeness_tl_ext_S20_pairs'] = dict(n=int(m.sum()), n_forms=int(len(set(form[m]))), slope=float(f.params[1]), se=float(f.bse[1]),
                                             p=float(2 * stats.t.sf(abs(f.params[1] / f.bse[1]), len(set(form[m])) - 1)))
Fo = X2[X2.S20_state == 'foreign'].copy()
Fo['bR'] = [blso(a, b) for a, b in zip(Fo['sample'], Fo.arm)]; Fo['bD'] = [blso(a, b) for a, b in zip(Fo['sample'], Fo.S20_best_end)]
Fo = Fo.dropna(subset=['bR', 'bD']); Fo['form'] = Fo.arm + '>' + Fo.S20_best_end
fu = sm.OLS(Fo.rel_tl_ext.values, sm.add_constant(Fo[['bR', 'bD']].values)).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(Fo.form)[0]})
out['donor_likeness_tl_ext_S20_unpaired'] = dict(n=len(Fo), n_forms=int(Fo.form.nunique()), beta_R=float(fu.params[1]), se_R=float(fu.bse[1]),
                                                beta_D=float(fu.params[2]), se_D=float(fu.bse[2]))
pf = []
for fm_ in sorted(set(form[h])):
    mm = h & (form == fm_)
    if mm.sum() >= 10:
        pf.append(dict(form=fm_, n=int(mm.sum()), mean_d_telo_len=float(pairs('S20_state')[1][mm].mean()), mean_d_tl_ext=float(dx[mm].mean()),
                       wilcoxon_p_tl_ext=float(stats.wilcoxon(dx[mm]).pvalue)))
PF = pd.DataFrame(pf); p_ = PF.wilcoxon_p_tl_ext.values; o_ = np.argsort(p_); q_ = p_[o_] * len(p_) / np.arange(1, len(p_) + 1)
PF.loc[PF.index[o_], 'bh_q_tl_ext'] = np.minimum(np.minimum.accumulate(q_[::-1])[::-1], 1)
PF.sort_values('wilcoxon_p_tl_ext').to_csv(os.path.join(RES, 'paired_by_form_tl_ext.tsv'), sep='\t', index=False, float_format='%.4g')
out['forms_tl_ext'] = dict(n_forms=len(PF), n_bh_q_lt_0_05=int((PF.bh_q_tl_ext < 0.05).sum()),
                           top=PF.sort_values('wilcoxon_p_tl_ext').head(3)[['form', 'n', 'mean_d_tl_ext', 'wilcoxon_p_tl_ext', 'bh_q_tl_ext']].to_dict('records'))
json.dump(out, open(os.path.join(RES, 'window_robustness.json'), 'w'), indent=1, default=float)
pd.set_option('display.width', 250)
print(R[['test', 'n', 'mean_d', 'ci95_lo', 'ci95_hi', 'wilcoxon_p', 'mde80']].round(3).to_string(index=False))
print(json.dumps(out, indent=1, default=lambda v: round(float(v), 4)))
