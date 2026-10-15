"""Figures and robustness summary for the telomere_cis track.
Outputs: figures/fig1_arm_baseline.{png,pdf}, figures/fig2_allele_test.{png,pdf}, figures/fig3_hg002_reads.{png,pdf},
         results/robustness_summary.json (arm ranking vs Karimian within HiFi input groups and with S20 own-form baseline)
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tl_common import RES, FIG, DAT, WORK, arm_order, ACRO_P

plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'font.family': 'DejaVu Sans', 'axes.titlelocation': 'left'})
C_OWN, C_FOR, C_GREY, C_FOCAL = '#4C72B0', '#DD8452', '#9a9a9a', '#222222'
PLAT = {'ont_ul_dorado': ('ONT UL', '#4C72B0'), 'revio_3cells': ('HiFi Revio', '#55A868'), 'hifi_dcv1.1': ('HiFi DeepConsensus', '#8172B2')}
def letter(ax, s): ax.text(-0.13, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
def save(fig, name):
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(FIG, f'{name}.{ext}'), dpi=300, bbox_inches='tight', metadata={'Software': None} if ext == 'png' else {'Creator': None, 'Producer': None, 'CreationDate': None})
    plt.close(fig)

T = pd.read_csv(os.path.join(RES, 'tip_telomere_lengths.tsv'), sep='\t', low_memory=False)
B = pd.read_csv(os.path.join(RES, 'arm_baseline.tsv'), sep='\t').set_index('arm')
K = pd.read_csv(os.path.join(DAT, 'karimian2024_fig4e.tsv'), sep='\t').set_index('arm')
Bo = pd.read_csv(os.path.join(RES, 'arm_own_baseline.tsv'), sep='\t').set_index('arm')
L = pd.read_csv(os.path.join(RES, 'layout_outer_read.tsv'), sep='\t')
bs = json.load(open(os.path.join(RES, 'baseline_stats.json')))
# ---------------- robustness numbers
rob = {}
j = K.join(Bo, how='inner'); r = stats.spearmanr(j.rel_mean_kb_digitized, j.own_mean_rel_tl)
rob['karimian_vs_hprc_own_S20'] = dict(n=len(j), rho=round(float(r.statistic), 4), p=float(r.pvalue))
for k, d in T.groupby('hifi_input'):
    m = d.groupby('arm').rel_tl.agg(['mean', 'size']); m = m[m['size'] >= 20].join(K, how='inner')
    if len(m) >= 10:
        r = stats.spearmanr(m.rel_mean_kb_digitized, m['mean']); rob[f'karimian_vs_hprc[hifi_input={k}]'] = dict(n_arms=len(m), rho=round(float(r.statistic), 4), p=float(r.pvalue))
m = T[T.assembler == 'hifiasm'].groupby('arm').rel_tl.mean().to_frame('m').join(K, how='inner')
rob['karimian_vs_hprc[hifiasm_only]'] = dict(n=len(m), rho=round(float(stats.spearmanr(m.rel_mean_kb_digitized, m.m).statistic), 4))
# window artefact: relative TL ~ arm + state (own reference) under the published-window replica (T20R) and S20
import statsmodels.formula.api as smf
FWs = pd.read_csv(os.path.join(RES, 'fixed_window_states.tsv'), sep='\t')
TS = T.merge(FWs[['q', 'T20R_state', 'S20_state']], on='q')
for t in ['T20R', 'S20']:
    m = smf.ols(f'rel_tl ~ C(arm) + C({t}_state, Treatment("own"))', data=TS).fit()
    rob[f'state_effects_bp[{t}]'] = {k.split('[T.')[1].rstrip(']'): dict(beta=round(float(m.params[k]), 1), se=round(float(m.bse[k]), 1),
                                     n=int((TS[f'{t}_state'] == k.split('[T.')[1].rstrip(']')).sum())) for k in m.params.index if f'{t}_state' in k}
json.dump(rob, open(os.path.join(RES, 'robustness_summary.json'), 'w'), indent=1)

# ---------------- Figure 1: baseline
fig = plt.figure(figsize=(7.2, 6.0)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.62, wspace=0.55)
ax = fig.add_subplot(gs[0, :]); letter(ax, 'a')
b = B[B.n >= 20].sort_values('mean_rel_tl'); x = np.arange(len(b))
cols = [C_FOCAL if a in K.index else C_GREY for a in b.index]
ax.axhline(0, color='#cccccc', lw=0.6, zorder=0)
ax.errorbar(x, b.mean_rel_tl / 1000, yerr=1.96 * b.se_rel_tl / 1000, fmt='none', ecolor='#bbbbbb', lw=0.8, zorder=1)
ax.scatter(x, b.mean_rel_tl / 1000, c=cols, s=12, zorder=2)
ax.set_xticks(x); ax.set_xticklabels(b.index, rotation=90); ax.set_xlim(-0.8, len(b) - 0.2)
ax.set_ylabel('Assembled telomere length\nrelative to sample mean (kb)')
ax.set_title(f'Chromosome arms differ in assembled telomere length (HPRC Release 2, {int(B.n.sum()):,} tips, {T["sample"].nunique()} samples; mean ± 95% CI)')
ax.text(0.01, 0.93, 'black: arm in Karimian et al. 2024 Fig. 4e;  grey: not reported there', transform=ax.transAxes, fontsize=6, color='#555555')
ax = fig.add_subplot(gs[1, 0]); letter(ax, 'b')
c = K.join(B[['mean_rel_tl']], how='inner')
ax.scatter(c.rel_mean_kb_digitized, c.mean_rel_tl / 1000, s=10, color=C_FOCAL)
for a in ['17p', '20q', '12p', '3p', '12q', '4q']:
    ax.annotate(a, (c.rel_mean_kb_digitized[a], c.mean_rel_tl[a] / 1000), xytext=(3, 2), textcoords='offset points', fontsize=6)
rr = bs['mean_rel_tl']; ax.set_title(f'Arm ranking matches nanopore\nprofiling (ρ = {rr["spearman_rho"]:.2f}, {c.shape[0]} arms)')
ax.set_xlabel('Karimian et al. 2024 relative\nmean length (kb, digitized)'); ax.set_ylabel('HPRC relative length (kb)')
lr = stats.linregress(c.rel_mean_kb_digitized, c.mean_rel_tl / 1000); xx = np.array(ax.get_xlim()); ax.plot(xx, lr.intercept + lr.slope * xx, color=C_GREY, lw=0.8)
ax.margins(0.06)
ax = fig.add_subplot(gs[1, 1]); letter(ax, 'c')
order = ['Revio only', 'DeepConsensus only', 'DeepConsensus + Revio', 'other/standard CCS', 'unknown']
data = [T.loc[T.hifi_input == k, 'telo_len'] / 1000 for k in order]
bp = ax.boxplot(data, widths=0.6, showfliers=False, patch_artist=True, medianprops=dict(color='k', lw=1))
for p_ in bp['boxes']: p_.set_facecolor('#dddddd'); p_.set_edgecolor('#555555')
ax.set_xticks(range(1, 6)); ax.set_xticklabels(['Revio\nonly', 'DC\nonly', 'DC +\nRevio', 'other\nCCS', 'unknown'])
ax.set_xticklabels([f'{l}\n{len(d)/1000:.1f}k' for l, d in zip(['Revio\nonly', 'DC\nonly', 'DC +\nRevio', 'other\nCCS', 'unk-\nnown'], data)], fontsize=5.5)
ax.set_ylabel('Assembled telomere length (kb)'); ax.set_title('Absolute length depends\non HiFi input batch')
ax = fig.add_subplot(gs[1, 2]); letter(ax, 'd')
ax.scatter(L.outer_read_len / 1000, L.telo_len / 1000, s=4, alpha=0.5, color=C_FOCAL, lw=0)
xx = np.array([0, L.outer_read_len.max() / 1000 * 1.05]); ax.plot(xx, xx, color=C_GREY, lw=0.8, ls='--'); ax.text(xx[1] * 0.62, xx[1] * 0.70, 'telomere = read', fontsize=6, color='#555555', rotation=0)
ax.set_xlim(0, xx[1]); ax.set_ylim(0, xx[1])
rho = stats.spearmanr(L.telo_len, L.outer_read_len).statistic
ax.set_xlabel('Tip-carrying read length (kb)'); ax.set_ylabel('Assembled telomere length (kb)')
ax.set_title(f'Telomeres lie well below\ntip-read length (ρ = {rho:.2f})')
save(fig, 'fig1_arm_baseline')

# ---------------- Figure 2: allele test
fig = plt.figure(figsize=(7.2, 8.4)); gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.7], hspace=0.42, wspace=0.95)
ax = fig.add_subplot(gs[0, 0]); letter(ax, 'a')
WR = pd.read_csv(os.path.join(RES, 'window_robustness.tsv'), sep='\t').set_index('test')
defs = [('T20R', 'T20R | boundary-robust TL_ext', 'published tip[s0, 20 kb)\n(TL-dependent)'), ('S20', 'S20 | boundary-robust TL_ext', 'S20: [s0, s0+20 kb)'),
        ('S10', 'S10 | boundary-robust TL_ext', 'S10: [s0, s0+10 kb)'), ('S5', 'S5 | boundary-robust TL_ext', 'S5: [s0, s0+5 kb)'),
        ('W515', 'W515 | boundary-robust TL_ext', 'W5-15: [s0+5, s0+15 kb)'),
        ('concordant S5=S10=S20', 'concordant S5=S10=S20 | TL_ext', 'same state at\nS5, S10 and S20')]
for i_, (raw, ext, lab) in enumerate(defs):
    for key, off, mk, fc in [(raw, -0.14, 'o', None), (ext, 0.14, 'o', 'white')]:
        r = WR.loc[key]; col = C_FOR if raw == 'T20R' else C_FOCAL
        ax.errorbar(r.mean_d / 1000, -i_ + off, xerr=[[(r.mean_d - r.ci95_lo) / 1000], [(r.ci95_hi - r.mean_d) / 1000]], fmt=mk, color=col, ms=3.5, lw=0.9,
                    mfc=fc if fc else col)
    rr = WR.loc[raw]; re_ = WR.loc[ext]
    ax.text(1.02, -i_, f'n={int(rr.n)}\np={rr.wilcoxon_p:.1g} / {re_.wilcoxon_p:.1g}', fontsize=5.5, va='center', transform=ax.get_yaxis_transform())
ax.axvline(0, color='#cccccc', lw=0.6, zorder=0); ax.set_yticks([-i_ for i_ in range(len(defs))]); ax.set_yticklabels([d[2] for d in defs], fontsize=5.5)
ax.set_ylim(-len(defs) + 0.5, 0.5); ax.set_xlabel('Foreign − own haplotype\ntelomere length (kb)')
ax.set_title('Negative shift disappears with a\nboundary-robust telomere length', loc='left')
from matplotlib.lines import Line2D as _L
ax.legend(handles=[_L([], [], marker='o', ls='', color=C_FOCAL, label='called TL'), _L([], [], marker='o', ls='', color=C_FOCAL, mfc='white', label='TL_ext')],
          frameon=False, fontsize=5.5, loc='upper left', bbox_to_anchor=(-0.05, -0.30), ncol=2, handletextpad=0.2)
ax = fig.add_subplot(gs[1, 0:2]); letter(ax, 'b')
F2 = pd.read_csv(os.path.join(RES, 'paired_by_form.tsv'), sep='\t').sort_values('mean_d')
y = np.arange(len(F2))
ax.errorbar(F2.mean_d / 1000, y, xerr=[(F2.mean_d - F2.ci95_lo) / 1000, (F2.ci95_hi - F2.mean_d) / 1000], fmt='o', ms=3, color=C_FOCAL, ecolor='#bbbbbb', lw=0.8)
sig = F2.bh_q < 0.05; ax.scatter(F2.mean_d[sig] / 1000, y[sig], s=22, facecolor='none', edgecolor=C_FOR, lw=1, zorder=3, label='BH q < 0.05')
ax.set_yticks(y); ax.set_yticklabels([f'{f} (n={n})' for f, n in zip(F2.form, F2.n)], fontsize=6); ax.axvline(0, color='#cccccc', lw=0.6, zorder=0)
ax.set_xlabel('Foreign − own haplotype telomere length (kb, mean, 95% bootstrap CI)')
ax.set_title('Per-form paired differences (S20; forms with ≥ 10 heterozygous individuals)')
if sig.any(): ax.legend(frameon=False, loc='lower right')
ax.margins(y=0.02)
ax = fig.add_subplot(gs[0, 1]); letter(ax, 'c')
P = pd.read_csv(os.path.join(RES, 'het_pairs.tsv'), sep='\t')
fm = P.groupby('form').agg(n=('d', 'size'), d=('d', 'mean'), x=('x', 'mean')).dropna(); fm = fm[fm.n >= 5]
ax.scatter(fm.x / 1000, fm.d / 1000, s=np.clip(fm.n, 5, 60) * 0.8, color=C_FOCAL, alpha=0.7, lw=0)
js = json.load(open(os.path.join(RES, 'paired_stats.json')))
t4 = [t for t in js['T4_donor_likeness'] if t['test'].startswith('form-level')][0]
xx = np.array([fm.x.min(), fm.x.max()]) / 1000
ax.plot(xx, t4['intercept'] / 1000 + t4['slope'] * xx, color=C_FOCAL, lw=1)
ax.plot(xx, xx, color=C_FOR, lw=0.8, ls='--'); ax.text(xx[1], xx[1], ' slope 1\n (donor)', fontsize=6, color=C_FOR, va='center')
ax.axhline(0, color='#cccccc', lw=0.6, zorder=0)
ax.set_xlabel('Donor − recipient\narm baseline (kb)'); ax.set_ylabel('Foreign − own (kb, form mean)')
ax.set_title(f'Weak donor-likeness, called TL\n(slope {t4["slope"]:.2f} ± {t4["slope_se_cluster"]:.2f}, {t4["n"]} forms)')
ax = fig.add_subplot(gs[1, 2]); letter(ax, 'e')
FW = pd.read_csv(os.path.join(RES, 'fixed_window_states.tsv'), sep='\t')
TT = T.merge(FW, on='q')
for t, col, lab, off in [('T20R', C_FOR, 'published window', -0.12), ('S20', C_FOCAL, 'S20 window', 0.12)]:
    o = TT[TT[f'{t}_state'] == 'own'].copy(); o['nn'] = np.minimum(o[f'{t}_n_near'], 3)
    g = o.groupby('nn').rel_tl.agg(['mean', 'sem', 'size'])
    ax.errorbar(g.index + off, g['mean'] / 1000, yerr=1.96 * g['sem'] / 1000, fmt='o-', color=col, ms=3, lw=1, label=lab)
ax.axhline(0, color='#cccccc', lw=0.6, zorder=0); ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(['0', '1', '2', '≥3'])
ax.set_xlabel('Other CHM13 ends within\n2% of best score'); ax.set_ylabel('Relative telomere length (kb)')
ax.set_title('Paralogy–length link is a\nwindow artefact (own tips)'); ax.legend(frameon=False, loc='upper left')
ax = fig.add_subplot(gs[0, 2]); letter(ax, 'd')
labs = []; vals = []; errs = []
for t, suf in [('t20', '_t20'), ('S20', ''), ('S10', '_S10'), ('S5', '_S5')]:
    jj = json.load(open(os.path.join(RES, f'paired_stats{suf}.json')))
    for tt in jj['T4_donor_likeness']:
        if tt['test'].startswith('all foreign tips: rel_tl ~ base'):
            vals.append((tt['beta_R'], tt['beta_D'])); errs.append((tt['se_R'], tt['se_D'])); labs.append(t)
wr = json.load(open(os.path.join(RES, 'window_robustness.json')))['donor_likeness_tl_ext_S20_unpaired']
vals.append((wr['beta_R'], wr['beta_D'])); errs.append((wr['se_R'], wr['se_D'])); labs.append('S20\nTL_ext')
vals = np.array(vals); errs = np.array(errs); yy = np.arange(len(labs))
ax.errorbar(vals[:, 1], yy + 0.12, xerr=1.96 * errs[:, 1], fmt='o', color=C_FOR, ms=4, lw=1, label='donor arm (β_D)')
ax.errorbar(vals[:, 0], yy - 0.12, xerr=1.96 * errs[:, 0], fmt='s', color=C_OWN, ms=4, lw=1, label='recipient arm (β_R)')
ax.axvline(0, color='#cccccc', lw=0.6, zorder=0); ax.axvline(1, color='#cccccc', lw=0.6, ls='--', zorder=0)
ax.set_yticks(yy); ax.set_yticklabels(labs); ax.invert_yaxis(); ax.set_xlabel('Weight on arm baseline\n(95% CI)')
ax.set_ylim(len(labs) - 0.4, -1.4)
ax.set_title('Foreign ends: donor vs\nrecipient baseline'); ax.legend(frameon=False, loc='upper center', fontsize=6, ncol=1)
save(fig, 'fig2_allele_test')

# ---------------- Figure 3: HG002
X = pd.read_csv(os.path.join(RES, 'hg002_ends_summary.tsv'), sep='\t').set_index('end')
Rd = pd.read_csv(os.path.join(WORK, 'hg002_reads.tsv.gz'), sep='\t'); Rd = Rd[(Rd.status == 'complete') & (Rd.mapq >= 20)]
Hh = pd.read_csv(os.path.join(RES, 'hg002_het_arms.tsv'), sep='\t'); Hh = Hh[Hh.allele_definition == 'S20']
hs = json.load(open(os.path.join(RES, 'hg002_stats.json')))
fig = plt.figure(figsize=(7.2, 5.8)); gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.6)
ax = fig.add_subplot(gs[0, :]); letter(ax, 'a')
o = Rd[Rd.platform == 'ont_ul_dorado']; med = o.groupby('end').tl.median(); cnt = o.groupby('end').size()
ends = med[cnt >= 5].sort_values().index; xi = {e: i for i, e in enumerate(ends)}
st = X.S20_state.reindex(ends)
rng = np.random.default_rng(1)
for e in ends:
    v = o.loc[o.end == e, 'tl'] / 1000; colr = C_FOR if st[e] == 'foreign' else (C_OWN if st[e] == 'own' else C_GREY)
    ax.scatter(xi[e] + rng.uniform(-0.25, 0.25, len(v)), v, s=1.5, color=colr, alpha=0.35, lw=0)
    ax.plot([xi[e] - 0.35, xi[e] + 0.35], [med[e] / 1000] * 2, color='k', lw=0.9)
ax.set_xticks(range(len(ends))); ax.set_xticklabels([e.split('|')[1].split('_')[0].replace('chr', '') + e.split('|')[2] + ('-M' if 'MATERNAL' in e else '-P') for e in ends], rotation=90, fontsize=4.5)
ax.set_ylabel('Per-read telomere length (kb)'); ax.set_xlim(-1, len(ends))
ax.set_ylim(0, np.percentile(o.tl, 99.7) / 1000)
ax.set_title(f'HG002 ONT ultra-long reads: end-specific telomere length ({len(o):,} complete reads, {len(ends)} ends with ≥ 5 reads; black = median)')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=C_FOR, label='foreign (S20)'), Line2D([], [], marker='o', ls='', color=C_OWN, label='own'),
                   Line2D([], [], marker='o', ls='', color=C_GREY, label='tie')], frameon=False, loc='upper left', ncol=3)
ax = fig.add_subplot(gs[1, 0]); letter(ax, 'b')
ax.scatter(X.ont_ul_dorado / 1000, X['hifi_dcv1.1'] / 1000, s=8, color=PLAT['hifi_dcv1.1'][1], label='HiFi DeepConsensus')
ax.scatter(X.ont_ul_dorado / 1000, X.revio_3cells / 1000, s=8, color=PLAT['revio_3cells'][1], marker='s', label='HiFi Revio')
cc = hs['concordance']
ax.set_xlabel('ONT UL end median (kb)'); ax.set_ylabel('HiFi end median (kb)')
ax.set_title(f'End profiles differ by DNA\nbatch (ρ {cc["ont_vs_hifi_dc"]["rho"]:.2f} DC, {cc["ont_vs_revio"]["rho"]:.2f} Revio)')
ax.legend(frameon=False, loc='upper left', fontsize=6); ax.margins(0.05)
ax = fig.add_subplot(gs[1, 1]); letter(ax, 'c')
k = 0; ticks = []
for r in Hh.itertuples():
    for e, col, off in [(r.own_end, C_OWN, -0.18), (r.foreign_end, C_FOR, 0.18)]:
        v = o.loc[o.end == e, 'tl'] / 1000
        ax.scatter(k + off + rng.uniform(-0.08, 0.08, len(v)), v, s=3, color=col, alpha=0.6, lw=0)
        if len(v): ax.plot([k + off - 0.12, k + off + 0.12], [v.median()] * 2, color='k', lw=1)
    ticks.append(f'{r.arm}\n>{r.donor}'); k += 1
ax.set_xticks(range(len(ticks))); ax.set_xticklabels(ticks, fontsize=6)
ax.set_ylabel('Per-read telomere length (kb, ONT)')
het = hs['het_arms']['S20|ont_ul_dorado']
ax.set_title(f'Het arms: foreign longer\nat {het["n_foreign_longer"]}/{het["n_arms"]} (sign p = {het["sign_test_p"]:.2g})')
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=C_OWN, label='own haplotype'), Line2D([], [], marker='o', ls='', color=C_FOR, label='foreign haplotype')], frameon=False, fontsize=6, loc='upper right')
ax = fig.add_subplot(gs[1, 2]); letter(ax, 'd')
X['arm'] = X.arm.astype(str)
A = X[~X.arm.isin(ACRO_P) & ~X.arm.str.startswith(('X', 'Y'))].groupby('arm').rel_ont_ul_dorado.mean().to_frame('hg').join(B.mean_rel_tl)
ax.scatter(A.mean_rel_tl / 1000, A.hg / 1000, s=8, color=C_FOCAL)
ar = hs['arm_ranking']['ont_vs_hprc_all']
ax.set_xlabel('HPRC arm baseline\n(kb, relative)'); ax.set_ylabel('HG002 ONT arm mean\n(kb, relative)')
ax.set_title(f'No match to population\narm ranking (ρ = {ar["rho"]:.2f}, {ar["n"]} arms)'); ax.margins(0.05)
save(fig, 'fig3_hg002_reads')
print(json.dumps(rob, indent=1))
