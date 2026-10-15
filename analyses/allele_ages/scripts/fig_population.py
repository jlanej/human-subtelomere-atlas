#!/usr/bin/env python3
"""Figure: population distribution of foreign subtelomere forms (results of pop_distribution.py).
(a) carrier frequency by superpopulation for the 30 forms with most carriers (Wilson CIs in form_freq_by_superpop.tsv);
(b) fraction of junction clusters shared between continental-African and EUR/EAS/SAS samples, observed vs expected
    under random sample labelling, by cluster size. Output: figures/fig_population.{png,pdf}"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import comb

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'sans-serif',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
SP = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
R = pd.read_csv(os.path.join(RES, 'form_popdiff_tests.tsv'), sep='\t').sort_values('carriers', ascending=False).head(30)
J = pd.read_csv(os.path.join(RES, 'junction_cluster_superpop.tsv'), sep='\t')
H = pd.read_csv(os.path.join(RES, 'haplotypes_by_superpop.tsv'), sep='\t').set_index('superpopulation')
import json
summ = json.load(open(os.path.join(RES, 'pop_summary.json')))

fig = plt.figure(figsize=(7.0, 5.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.55)
ax = fig.add_subplot(gs[0])
M = R[[f'freq_{s}' for s in SP]].values
im = ax.imshow(M, cmap='viridis', vmin=0, vmax=0.8, aspect='auto')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=5.5, color='white' if M[i, j] < 0.45 else 'black')
labels = [f"{f}{' *' if q < 0.05 else ''}" for f, q in zip(R.form, R.q_perm_batch_stratified)]
ax.set_yticks(range(len(R))); ax.set_yticklabels(labels)
ax.set_xticks(range(5)); ax.set_xticklabels([f'{s}\n({int(H.loc[s, "haplotypes"])})' for s in SP])
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax2 = ax.twinx(); ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(range(len(R)))
ax2.set_yticklabels([f'{x:.2f}' for x in R.fst_anova_5pop]); ax2.tick_params(length=0)
for s in ax2.spines.values():
    s.set_visible(False)
ax2.set_ylabel('$F_{ST}$ (5 superpopulations)')
ax.set_title('Common forms are present in all superpopulations\nbut many differ in frequency', loc='left')
cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.09)
cb.set_label('carrier frequency (haplotypes)\n* superpopulation difference q < 0.05 (chemistry-stratified permutation)')

bx = fig.add_subplot(gs[1])
nC, nO, nT = summ['samples_AFR_continental'], summ['samples_EUR_EAS_SAS'], summ['samples']
def p_both(s):
    return 1 - comb(nT - nC, s) / comb(nT, s) - comb(nT - nO, s) / comb(nT, s) + comb(nT - nC - nO, s) / comb(nT, s)
bins = [(2, 2, '2'), (3, 4, '3-4'), (5, 9, '5-9'), (10, 999, '≥10')]
obs, exp, ns = [], [], []
for lo, hi, _ in bins:
    g = J[(J.samples >= lo) & (J.samples <= hi)]
    obs.append(g.shared_AFRcont_EUR_EAS_SAS.mean()); exp.append(np.mean([p_both(s) for s in g.samples])); ns.append(len(g))
x = np.arange(len(bins))
bx.bar(x - 0.18, obs, 0.36, color='#2b6c8f', label='observed')
bx.bar(x + 0.18, exp, 0.36, color='#bdbdbd', label='random labelling')
for xi, o, n in zip(x, obs, ns):
    bx.text(xi - 0.18, o + 0.02, f'{o:.2f}', ha='center', fontsize=6)
bx.set_xticks(x); bx.set_xticklabels([f'{b[2]}\n(n={n})' for b, n in zip(bins, ns)])
bx.set_xlabel('samples sharing the junction cluster'); bx.set_ylabel('fraction of clusters in both continental\nAfrican and EUR/EAS/SAS samples')
bx.set_ylim(0, 1.08); bx.legend(frameon=False, loc='upper left')
bx.set_title('Most large junction clusters span the\nout-of-Africa split (below random expectation)', loc='left')
for a_, l in ((ax, 'a'), (bx, 'b')):
    a_.text(-0.02 if a_ is bx else -0.28, 1.08, l, transform=a_.transAxes, fontweight='bold', fontsize=10, va='bottom', ha='right')
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FIG, f'fig_population.{ext}'), dpi=300, bbox_inches='tight')
print('saved')
