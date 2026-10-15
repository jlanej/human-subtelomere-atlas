"""Step 3 tests - HG002 per-read telomere length: platform concordance, truncation, arm ranking vs the population baseline,
and the allele effect between haplotypes (heterozygous foreign/own arms) and versus donor arms.

Primary measure: per-end median of complete, MAPQ >= 20 ONT UL reads (>= 5 reads); PacBio Revio HiFi as a replicate.
Allele states: S20 (telomere-length-independent, primary) and the published t20 state, from results/hg002_fixed_window_states.tsv.
Inputs : results/hg002_end_telomere.tsv, results/hg002_fixed_window_states.tsv, work/telomere_cis/hg002_reads.tsv.gz (step 3),
         results/arm_baseline.tsv, results/arm_own_baseline.tsv (S20 own-form population baseline), data/karimian2024_fig4e.tsv
Outputs: results/hg002_ends_summary.tsv, results/hg002_het_arms.tsv, results/hg002_stats.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from tl_common import RES, DAT, WORK, ACRO_P, SEED

rng = np.random.default_rng(SEED); MINR = 5
E = pd.read_csv(os.path.join(RES, 'hg002_end_telomere.tsv'), sep='\t')
S = pd.read_csv(os.path.join(RES, 'hg002_fixed_window_states.tsv'), sep='\t')
Rd = pd.read_csv(os.path.join(WORK, 'hg002_reads.tsv.gz'), sep='\t')
Rd = Rd[(Rd.status == 'complete') & (Rd.mapq >= 20)]
Bl = pd.read_csv(os.path.join(RES, 'arm_baseline.tsv'), sep='\t').set_index('arm')
Bo = pd.read_csv(os.path.join(RES, 'arm_own_baseline.tsv'), sep='\t').set_index('arm')
K = pd.read_csv(os.path.join(DAT, 'karimian2024_fig4e.tsv'), sep='\t').set_index('arm').rel_mean_kb_digitized * 1000
W = E.pivot_table(index='end', columns='platform', values='median_tl')
N = E.pivot_table(index='end', columns='platform', values='n_complete_mq20')
W = W.where(N >= MINR)
X = S.set_index('end').join(W).join(N.add_prefix('n_'))
X['ref_s0'] = X.s0
for p in ['ont_ul_dorado', 'revio_3cells', 'hifi_dcv1.1']: X[f'rel_{p}'] = X[p] - X[p].mean()
X.to_csv(os.path.join(RES, 'hg002_ends_summary.tsv'), sep='\t', float_format='%.5g')
out = {}
def sp(a, b):
    m = a.notna() & b.notna(); r = stats.spearmanr(a[m], b[m]); return dict(n=int(m.sum()), rho=round(float(r.statistic), 4), p=float(r.pvalue))
out['n_ends_ge5_reads'] = {p: int(X[p].notna().sum()) for p in ['ont_ul_dorado', 'revio_3cells', 'hifi_dcv1.1']}
out['median_of_end_medians_bp'] = {p: float(X[p].median()) for p in ['ont_ul_dorado', 'revio_3cells', 'hifi_dcv1.1']}
out['n_complete_reads_mq20'] = Rd.groupby('platform').size().to_dict()
out['concordance'] = {'ont_vs_revio': sp(X.ont_ul_dorado, X.revio_3cells), 'ont_vs_hifi_dc': sp(X.ont_ul_dorado, X['hifi_dcv1.1']),
                      'revio_vs_hifi_dc': sp(X.revio_3cells, X['hifi_dcv1.1']), 'ont_vs_assembly_v1.1_telo_len': sp(X.ont_ul_dorado, X.ref_telo_len)}
# truncation: per-read TL vs read length within platform; fraction of reads where TL >= 50 % of read length
tr = {}
for p, g in Rd.groupby('platform'):
    r = stats.spearmanr(g.tl, g.read_len)
    tr[p] = dict(n=len(g), median_read_len=float(g.read_len.median()), spearman_tl_vs_read_len=round(float(r.statistic), 4), p=float(r.pvalue),
                 frac_tl_ge_half_read=round(float((g.tl >= 0.5 * g.read_len).mean()), 4),
                 within_end_spearman_tl_vs_len=round(float(np.nanmedian([stats.spearmanr(x.tl, x.read_len).statistic for _, x in g.groupby('end') if len(x) >= 8])), 4))
out['truncation'] = tr
# arm ranking: average of the two haplotype ends (relative ONT), non-acrocentric autosomal arms
X['arm'] = X.arm.astype(str)
A = X[~X.arm.isin(ACRO_P) & ~X.arm.str.startswith(('X', 'Y'))].groupby('arm')[['rel_ont_ul_dorado', 'rel_revio_3cells']].mean()
A = A.join(Bl[['mean_rel_tl']]).join(Bo[['own_mean_rel_tl']]).join(K.rename('karimian'))
out['arm_ranking'] = {'ont_vs_hprc_all': sp(A.rel_ont_ul_dorado, A.mean_rel_tl), 'ont_vs_hprc_own_S20': sp(A.rel_ont_ul_dorado, A.own_mean_rel_tl),
                      'ont_vs_karimian': sp(A.rel_ont_ul_dorado, A.karimian), 'revio_vs_hprc_all': sp(A.rel_revio_3cells, A.mean_rel_tl),
                      'revio_vs_karimian': sp(A.rel_revio_3cells, A.karimian)}
# per-end (haplotype-resolved) vs population own baseline of its own arm (own-state ends only)
o = X[X.S20_state == 'own']
out['own_ends_vs_hprc_own_baseline'] = sp(o.rel_ont_ul_dorado, o.arm.map(Bo.own_mean_rel_tl))
# ---- heterozygous arms
H = []
for tag in ['S20', 't20']:
    st = X[f'{tag}_state']; be = X[f'{tag}_best_end'].astype(str).str.replace('chr', '')
    for arm, g in X.groupby('arm'):
        if len(g) != 2: continue
        sts = list(st[g.index])
        if sorted(sts) != ['foreign', 'own']: continue
        fe = g.index[st[g.index] == 'foreign'][0]; oe = g.index[st[g.index] == 'own'][0]
        D = be[fe]
        rec = dict(allele_definition=tag, arm=arm, foreign_end=fe, own_end=oe, donor=D)
        for p in ['ont_ul_dorado', 'revio_3cells']:
            rf = Rd[(Rd.end == fe) & (Rd.platform == p)].tl.values; ro = Rd[(Rd.end == oe) & (Rd.platform == p)].tl.values
            rec[f'n_f_{p}'] = len(rf); rec[f'n_o_{p}'] = len(ro)
            if len(rf) >= MINR and len(ro) >= MINR:
                bs = [np.median(rng.choice(rf, len(rf))) - np.median(rng.choice(ro, len(ro))) for _ in range(2000)]
                rec.update({f'd_{p}': float(np.median(rf) - np.median(ro)), f'ci_lo_{p}': float(np.percentile(bs, 2.5)), f'ci_hi_{p}': float(np.percentile(bs, 97.5)),
                            f'mwu_p_{p}': float(stats.mannwhitneyu(rf, ro).pvalue)})
        rec['x_hprc'] = Bo.own_mean_rel_tl.get(D, np.nan) - Bo.own_mean_rel_tl.get(arm, np.nan)
        rec['x_karimian'] = K.get(D, np.nan) - K.get(arm, np.nan)
        H.append(rec)
H = pd.DataFrame(H); H.to_csv(os.path.join(RES, 'hg002_het_arms.tsv'), sep='\t', index=False, float_format='%.4g')
het = {}
for tag, g in H.groupby('allele_definition'):
    for p in ['ont_ul_dorado', 'revio_3cells']:
        d = g[f'd_{p}'].dropna()
        r = dict(n_arms=int(len(d)), arms=';'.join(g.loc[d.index, 'arm']), median_d=float(d.median()) if len(d) else None,
                 n_foreign_longer=int((d > 0).sum()), n_nominal_p_lt_0_05=int((g.loc[d.index, f'mwu_p_{p}'] < 0.05).sum()))
        if len(d) >= 3:
            r['sign_test_p'] = float(stats.binomtest(int((d > 0).sum()), len(d)).pvalue)
            r['wilcoxon_p'] = float(stats.wilcoxon(d).pvalue) if len(d) >= 5 else None
            m = g.loc[d.index]; ok = m.x_hprc.notna()
            if ok.sum() >= 4:
                lr = stats.linregress(m.x_hprc[ok], d[ok]); r.update(donor_slope_hprc=float(lr.slope), donor_slope_se=float(lr.stderr), donor_slope_p=float(lr.pvalue))
        het[f'{tag}|{p}'] = r
out['het_arms'] = het
# donor-likeness over all HG002 foreign ends (S20): rel ONT ~ base(R) + base(D)
F = X[(X.S20_state == 'foreign')].copy(); F['bR'] = F.arm.map(Bo.own_mean_rel_tl); F['bD'] = F.S20_best_end.astype(str).map(Bo.own_mean_rel_tl)
F = F.dropna(subset=['bR', 'bD', 'rel_ont_ul_dorado'])
if len(F) >= 6:
    Xm = np.column_stack([np.ones(len(F)), F.bR, F.bD]); y = F.rel_ont_ul_dorado.values; b = np.linalg.lstsq(Xm, y, rcond=None)[0]
    r = y - Xm @ b; se = np.sqrt(np.diag(r @ r / (len(y) - 3) * np.linalg.inv(Xm.T @ Xm)))
    out['hg002_foreign_ends_rel_tl_on_base'] = dict(n=len(F), beta_R=float(b[1]), se_R=float(se[1]), beta_D=float(b[2]), se_D=float(se[2]),
                                                   p_D=float(2 * stats.t.sf(abs(b[2] / se[2]), len(y) - 3)), p_R=float(2 * stats.t.sf(abs(b[1] / se[1]), len(y) - 3)))
# HG002 end-level foreign vs own (all ends, arm not controlled) as context
fo = X[X.S20_state.isin(['foreign', 'own'])]
out['hg002_ends_foreign_vs_own_mwu'] = dict(n_foreign=int((fo.S20_state == 'foreign').sum()), n_own=int((fo.S20_state == 'own').sum()),
    median_rel_foreign=float(fo[fo.S20_state == 'foreign'].rel_ont_ul_dorado.median()), median_rel_own=float(fo[fo.S20_state == 'own'].rel_ont_ul_dorado.median()),
    p=float(stats.mannwhitneyu(fo[fo.S20_state == 'foreign'].rel_ont_ul_dorado.dropna(), fo[fo.S20_state == 'own'].rel_ont_ul_dorado.dropna()).pvalue))
json.dump(out, open(os.path.join(RES, 'hg002_stats.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float)); print(H.round(3).to_string(index=False))
