"""Step 1 - arm-specific assembled telomere length baseline (HPRC Release 2 named tips) with truncation diagnostics, and
comparison with published arm-specific telomere length (Karimian et al. 2024, doi:10.1126/science.ado0431, Fig. 4e).

Inputs (ATLAS_INPUT_DIR): results/subtel_foreign_tips.tsv (11,998 named tips; telo_len = assembled terminal (CCCTAA)n block,
        caller pipeline/catalogue_telomere.py on 20 kb tip windows), subtel/foreign/r2_100k/<asm>.fa (100 kb windows, used to
        re-measure telomeres that reach the 20 kb window edge), meta/explorer_seqdata.csv (HiFi file N50/coverage),
        results/single_read_zone_by_tip.tsv + gfa/<asm>.Alines.tsv (hifiasm read layouts, 10 haplotypes),
        data/karimian2024_fig4e.tsv (digitized, scripts/00_karimian_digitize.py).
Outputs: results/tip_telomere_lengths.tsv, results/arm_baseline.tsv, results/sample_covariates.tsv,
        results/truncation_diagnostics.tsv, results/layout_outer_read.tsv, results/karimian_comparison.tsv,
        results/baseline_stats.json
Relative telomere length (rel_tl) = tip telo_len minus the mean telo_len of all named tips of that sample (both haplotypes),
as in Karimian et al. (deviation from the individual's grand mean); it removes sample-level effects (cell-line passage, age,
read-length distribution of the sequencing run).
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from tl_common import inp, RES, DAT, rc, find_block, arm_order, ACRO_P, load_tips, SEED

F = load_tips()
V = pd.read_csv(inp('results', 'subtel_foreign_tips.tsv'), sep='\t', usecols=['q', 'vk_verkko_telo_len', 'vk_ident', 'vk_verkko_hit'], low_memory=False)
F = F.merge(V, on='q', how='left')
# ---- telomeres reaching the 20 kb window edge: re-measure on 100 kb windows where available
F['censored20'] = (F.junk_len + F.telo_len) >= 19990
F['telo_len_raw20'] = F.telo_len
F['telo_source'] = 'w20'
need = F[F.censored20]
for asm, g in need.groupby('asm'):
    fa = inp('subtel', 'foreign', 'r2_100k', f'{asm}.fa')
    if not os.path.exists(fa): continue
    seqs = {}; name = None; buf = []
    for line in open(fa):
        if line[0] == '>':
            if name: seqs[name] = ''.join(buf)
            name = line[1:].strip(); buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    for r in g.itertuples():
        key = f'{r.contig}:1-100000' if r.end == 'L' else f'{r.contig}:{r.contig_len - 99999}-{r.contig_len}'
        if key not in seqs: continue
        s = seqs[key].upper(); s = s if r.end == 'L' else rc(s)
        b = find_block(s)
        if b and b['junk_len'] == r.junk_len:
            F.loc[r.Index, ['telo_len', 'telo_source']] = [b['telo_len'], 'w100']
F['censored'] = F.censored20 & (F.telo_source == 'w20')
# ---- sample covariates: coverage-weighted HiFi N50 (WGS PacBio files)
S = pd.read_csv(inp('meta', 'explorer_seqdata.csv'), low_memory=False)
P = S[(S.platform == 'PACBIO_SMRT') & (S.library_strategy == 'WGS')].copy()
P['n50'] = pd.to_numeric(P.n50, errors='coerce'); P['coverage'] = pd.to_numeric(P.coverage, errors='coerce')
P = P.dropna(subset=['n50', 'coverage'])
cov = P.groupby('sample_id').apply(lambda d: pd.Series(dict(hifi_n50=np.average(d.n50, weights=d.coverage), hifi_cov=d.coverage.sum(),
        revio_frac=d.coverage[d.instrument_model == 'Revio'].sum() / d.coverage.sum())), include_groups=False)
O = S[S.platform == 'OXFORD_NANOPORE'].copy(); O['c100'] = pd.to_numeric(O.coverage_over_100kb, errors='coerce')
cov['ont_cov_100kb'] = O.groupby('sample_id').c100.sum(min_count=1)
cov = cov.reset_index().rename(columns={'sample_id': 'sample'})
smp = F.groupby('sample').agg(n_tips=('q', 'size'), mean_tl=('telo_len', 'mean'), median_tl=('telo_len', 'median'),
        assembler=('assembler', 'first'), hifi_input=('hifi_input', 'first')).reset_index().merge(cov, on='sample', how='left')
smp.to_csv(os.path.join(RES, 'sample_covariates.tsv'), sep='\t', index=False, float_format='%.4g')
F = F.merge(smp[['sample', 'mean_tl', 'hifi_n50']], on='sample', how='left')
F['rel_tl'] = F.telo_len - F.mean_tl
F.drop(columns=['mean_tl']).to_csv(os.path.join(RES, 'tip_telomere_lengths.tsv'), sep='\t', index=False, float_format='%.5g')
# ---- arm baseline
def summ(d):
    return pd.Series(dict(n=len(d), n_samples=d['sample'].nunique(), median_tl=d.telo_len.median(), mean_tl=d.telo_len.mean(),
        q25=d.telo_len.quantile(.25), q75=d.telo_len.quantile(.75), mean_rel_tl=d.rel_tl.mean(), se_rel_tl=d.rel_tl.std() / np.sqrt(len(d)),
        frac_censored=d.censored.mean(), frac_foreign=(d.t20_state == 'foreign').mean()))
B = F.groupby('arm').apply(summ, include_groups=False)
Bo = F[F.t20_state == 'own'].groupby('arm').apply(summ, include_groups=False)[['n', 'median_tl', 'mean_rel_tl', 'se_rel_tl']].add_suffix('_own')
B = B.join(Bo); B = B.loc[sorted(B.index, key=arm_order)]
B['rank_rel_all'] = B.mean_rel_tl.rank(); B['rank_rel_own'] = B.mean_rel_tl_own.rank()
# split-half reproducibility of the arm ranking (random halves of samples)
rng = np.random.default_rng(SEED); sams = np.array(sorted(F['sample'].unique())); rhos = []
for _ in range(200):
    h = set(rng.choice(sams, len(sams) // 2, replace=False))
    a = F[F['sample'].isin(h)].groupby('arm').rel_tl.mean(); b = F[~F['sample'].isin(h)].groupby('arm').rel_tl.mean()
    j = a.index.intersection(b.index); j = [x for x in j if B.n[x] >= 40]
    rhos.append(stats.spearmanr(a[j], b[j]).statistic)
B.to_csv(os.path.join(RES, 'arm_baseline.tsv'), sep='\t', float_format='%.4g')
# arm effect size (one-way ANOVA on rel_tl, arms with >= 40 tips)
big = B.index[B.n >= 40]; g = [F.loc[F.arm == a, 'rel_tl'].values for a in big]
fst = stats.f_oneway(*g); allv = np.concatenate(g); grand = allv.mean()
eta2 = sum(len(x) * (x.mean() - grand) ** 2 for x in g) / ((allv - grand) ** 2).sum()
# ---- truncation diagnostics
D = []
def add(name, value, note=''): D.append(dict(diagnostic=name, value=value, note=note))
add('named_tips', len(F)); add('samples', F['sample'].nunique())
add('tips_reaching_20kb_window', int(F.censored20.sum())); add('remeasured_on_100kb', int((F.telo_source == 'w100').sum()))
add('still_censored_at_20kb', int(F.censored.sum()))
add('median_telo_len_bp', float(F.telo_len.median())); add('frac_telo_gt_15kb', round(float((F.telo_len > 15000).mean()), 4))
for k, d in F.groupby('hifi_input'): add(f'median_telo_len_bp[hifi_input={k}]', float(d.telo_len.median()), f'n={len(d)}')
kw = stats.kruskal(*[d.telo_len for _, d in F.groupby('hifi_input')]); add('kruskal_p_telo_len_by_hifi_input', float(kw.pvalue))
for k, d in F.groupby('assembler'): add(f'median_telo_len_bp[assembler={k}]', float(d.telo_len.median()), f'n={len(d)}')
s2 = smp.dropna(subset=['hifi_n50'])
r = stats.spearmanr(s2.hifi_n50, s2.mean_tl); add('spearman_sample_mean_tl_vs_hifi_n50', round(float(r.statistic), 4), f'p={r.pvalue:.3g}; n={len(s2)} samples')
r = stats.spearmanr(s2.hifi_cov, s2.mean_tl); add('spearman_sample_mean_tl_vs_hifi_cov', round(float(r.statistic), 4), f'p={r.pvalue:.3g}')
s3 = s2.dropna(subset=['ont_cov_100kb']); r = stats.spearmanr(s3.ont_cov_100kb, s3.mean_tl)
add('spearman_sample_mean_tl_vs_ont_cov_100kb', round(float(r.statistic), 4), f'p={r.pvalue:.3g}; n={len(s3)}')
# linear model: sample mean TL ~ N50 + input type
import numpy.linalg as la
X = pd.get_dummies(s2.hifi_input, drop_first=True, dtype=float); X.insert(0, 'n50_kb', s2.hifi_n50 / 1000); X.insert(0, 'const', 1.0)
beta, *_ = la.lstsq(X.values, s2.mean_tl.values, rcond=None); resid = s2.mean_tl.values - X.values @ beta
se = np.sqrt(np.diag(resid @ resid / (len(s2) - X.shape[1]) * la.inv(X.values.T @ X.values)))
add('lm_sample_mean_tl_per_kb_hifi_n50_bp', round(float(beta[1]), 1), f'SE {se[1]:.1f}; adjusted for hifi_input')
# hifiasm vs Verkko (year-2 trio) telomere length on the same haplotype end
vv = F[(F.vk_verkko_hit == True) & F.vk_verkko_telo_len.notna() & (F.vk_ident >= 0.995)] if 'vk_verkko_hit' in F else F.iloc[:0]
if len(vv):
    w = stats.wilcoxon(vv.vk_verkko_telo_len, vv.telo_len)
    add('verkko_vs_hifiasm_same_end_n', len(vv)); add('verkko_minus_hifiasm_median_bp', float((vv.vk_verkko_telo_len - vv.telo_len).median()), f'Wilcoxon p={w.pvalue:.3g}')
    add('spearman_verkko_vs_hifiasm_telo_len', round(float(stats.spearmanr(vv.vk_verkko_telo_len, vv.telo_len).statistic), 4))
# hifiasm read layouts: length of the outermost read that carries the contig tip
Z = pd.read_csv(inp('results', 'single_read_zone_by_tip.tsv'), sep='\t'); rows = []
for asm, g2 in Z.groupby('asm'):
    fn = inp('gfa', f'{asm}.Alines.tsv')
    if not os.path.exists(fn): continue
    Al = pd.read_csv(fn, sep='\t', header=None, usecols=[1, 2, 4, 5, 6], names=['ctg', 'cpos', 'read', 'rs', 're'])
    span = (Al.re - Al.rs).groupby([Al.ctg, Al.read]).max()
    for r in g2.itertuples():
        rows.append(dict(asm=asm, contig=r.contig, end=r.end, telo_len=r.telo_len, junk_len=r.junk_len, one_read_zone=r.zone,
                         outer_read_len=span.get((r.contig, r.outer_read), np.nan)))
L = pd.DataFrame(rows).dropna(subset=['outer_read_len']); L['telo_frac_of_outer_read'] = L.telo_len / L.outer_read_len
L.to_csv(os.path.join(RES, 'layout_outer_read.tsv'), sep='\t', index=False, float_format='%.4g')
r = stats.spearmanr(L.telo_len, L.outer_read_len)
add('layout_tips', len(L), f'{L.asm.nunique()} hifiasm haplotypes with read layouts')
add('layout_median_outer_read_len_bp', float(L.outer_read_len.median()))
add('layout_spearman_telo_len_vs_outer_read_len', round(float(r.statistic), 4), f'p={r.pvalue:.3g}')
add('layout_frac_telo_ge_0.8_of_outer_read', round(float((L.telo_frac_of_outer_read >= 0.8).mean()), 4))
add('layout_frac_telomere_inside_one_read_zone', round(float((L.one_read_zone >= L.telo_len + L.junk_len).mean()), 4), 'whole telomere supported by a single read')
Dd = pd.DataFrame(D); Dd.to_csv(os.path.join(RES, 'truncation_diagnostics.tsv'), sep='\t', index=False)
# ---- comparison with Karimian et al. 2024
K = pd.read_csv(os.path.join(DAT, 'karimian2024_fig4e.tsv'), sep='\t').set_index('arm')
C = K.join(B[['n', 'mean_rel_tl', 'se_rel_tl', 'median_tl', 'mean_rel_tl_own', 'n_own']], how='inner')
C.to_csv(os.path.join(RES, 'karimian_comparison.tsv'), sep='\t', float_format='%.4g')
res = dict(n_arms_compared=len(C))
for col in ['mean_rel_tl', 'mean_rel_tl_own']:
    c = C.dropna(subset=[col])
    sp = stats.spearmanr(c.rel_mean_kb_digitized, c[col]); pe = stats.pearsonr(c.rel_mean_kb_digitized, c[col] / 1000)
    sl = stats.linregress(c.rel_mean_kb_digitized, c[col] / 1000)
    res[col] = dict(n=len(c), spearman_rho=round(float(sp.statistic), 4), spearman_p=float(sp.pvalue), pearson_r=round(float(pe.statistic), 4),
                    pearson_p=float(pe.pvalue), slope_hprc_kb_per_karimian_kb=round(float(sl.slope), 3))
kshort = list(K.index[:3]); klong = list(K.index[-3:])
res['hprc_rank_of_karimian_shortest3'] = {a: int(C.mean_rel_tl.rank()[a]) for a in kshort}
res['hprc_rank_of_karimian_longest3'] = {a: int(C.mean_rel_tl.rank()[a]) for a in klong}
res['arm_anova'] = dict(arms=int(len(big)), F=round(float(fst.statistic), 2), p=float(fst.pvalue), eta2=round(float(eta2), 4))
res['split_half_spearman_median'] = round(float(np.median(rhos)), 4); res['split_half_spearman_5_95'] = [round(float(x), 4) for x in np.percentile(rhos, [5, 95])]
res['sd_rel_tl_bp'] = round(float(F.rel_tl.std()), 1)
json.dump(res, open(os.path.join(RES, 'baseline_stats.json'), 'w'), indent=1)
print(Dd.to_string(index=False)); print(json.dumps(res, indent=1))
