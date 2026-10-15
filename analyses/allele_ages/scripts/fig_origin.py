#!/usr/bin/env python3
"""Figure: single origin and relative age of foreign forms (results of origin_age.py, origin_summary.py, ape_forms.py,
ape_junction_check.py). Output: figures/fig_origin.{png,pdf}"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                     'legend.fontsize': 6.5, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
V = pd.read_csv(os.path.join(RES, 'origin_verdicts.tsv'), sep='\t'); V = V[V.status == 'ok']
M = pd.read_csv(os.path.join(RES, 'origin_members.tsv'), sep='\t')
P = pd.read_csv(os.path.join(RES, 'ape_junction_probes.tsv'), sep='\t')
D = pd.read_csv(os.path.join(RES, 'ape_junction_hits_diagnosis.tsv'), sep='\t')
VC = {'consistent with single origin': '#1b7837', 'clustered, low diversity, no carrier-fixed SNVs': '#7fbf7b',
      'clustered, carriers as diverse as natives': '#9970ab', 'not supported': '#bababa'}
GC = {'carrier': '#d95f02', 'native': '#1f78b4', 'other_recipient': '#8c8c8c'}
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, hspace=0.85, wspace=0.55)
# a
ax = fig.add_subplot(gs[0, 0])
for v, c in VC.items():
    g = V[V.single_origin == v]
    ax.scatter(g.pi_native, g.pi_carrier, s=np.sqrt(g.carriers_used) * 4, color=c, edgecolor='k', lw=0.3, label=f'{v} ({len(g)})', zorder=3)
lim = (0.1, 12); ax.plot(lim, lim, color='k', lw=0.6, ls='--', zorder=1); ax.text(0.3, 0.37, 'equal diversity', fontsize=6, rotation=38, ha='center')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(0.08, 12)
ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(['0.1', '0.3', '1', '3', '10']); ax.set_yticks([0.1, 0.3, 1, 3, 10]); ax.set_yticklabels(['0.1', '0.3', '1', '3', '10'])
for f in ['20q>3q', '6p>11p', '1p>11p', '3q>11p', '8p>1p', '16p>9q']:
    r = V[V.form == f].iloc[0]
    ax.annotate(f, (r.pi_native, r.pi_carrier), xytext=(4, -2), textcoords='offset points', fontsize=6)
ax.set_xlabel('native donor-arm diversity, $\\pi$ (SNVs/kb)'); ax.set_ylabel('carrier diversity, $\\pi$ (SNVs/kb)')
ax.set_title(f'Carriers of {int((V.pi_carrier < V.pi_native).sum())} of {len(V)} forms are less diverse\nthan native donor alleles', loc='left')
ax.legend(frameon=False, loc='upper left', fontsize=5.5, handletextpad=0.2, ncol=2, bbox_to_anchor=(-0.12, -0.17), columnspacing=0.6)
# b
bx = fig.add_subplot(gs[0, 1])
forms = ['20q>3q', '3q>11p', '6p>11p', '1p>11p']
rng = np.random.default_rng(1)
for i, f in enumerate(forms):
    g = M[M.form == f]
    for j, grp in enumerate(['carrier', 'native', 'other_recipient']):
        y = g[g.group == grp].dist_to_carrier_consensus_per_kb.values
        x = i + (j - 1) * 0.27 + rng.uniform(-0.08, 0.08, len(y))
        bx.scatter(x, y + 0.02, s=3, color=GC[grp], alpha=0.6, lw=0, label=grp.replace('_', ' ') + 's' if i == 0 else None)
        if len(y):
            bx.plot([i + (j - 1) * 0.27 - 0.1, i + (j - 1) * 0.27 + 0.1], [np.median(y) + 0.02] * 2, color='k', lw=1)
bx.set_yscale('log'); bx.set_xticks(range(len(forms))); bx.set_xticklabels(forms)
bx.set_yticks([0.02, 0.1, 1, 10]); bx.set_yticklabels(['0', '0.1', '1', '10'])
bx.set_ylabel('distance to carrier consensus\n(SNVs/kb, +0.02)')
bx.set_title('Carriers sit on one tight lineage; donor lineage\nsurvives in natives (20q>3q, 3q>11p) or not (6p>11p)', loc='left')
bx.legend(frameon=False, loc='upper center', ncol=3, markerscale=3, bbox_to_anchor=(0.5, 1.0), columnspacing=0.8)
bx.set_ylim(0.015, 40)
# c
cx = fig.add_subplot(gs[1, 0])
W = V[V.single_origin != 'not supported'].sort_values('age_kyr_median_mu')
y = np.arange(len(W))
cx.scatter(W.age_kyr_median_mu.clip(lower=5), y, color=[VC[v] for v in W.single_origin], edgecolor='k', lw=0.3, s=18, zorder=3)
cx.scatter(W.age_kyr_mu, y, marker='|', color='k', s=25, lw=0.8, zorder=2, label='mean-based ($\\rho$) age')
for yi, (a, b) in enumerate(zip(W.age_kyr_median_mu.clip(lower=5), W.age_kyr_mu)):
    cx.plot([a, b], [yi, yi], color='#999999', lw=0.5, zorder=1)
cx.set_yticks(y); cx.set_yticklabels([f + (' (0; at 5)' if a == 0 else '') for f, a in zip(W.form, W.age_kyr_median_mu)], fontsize=5.5); cx.set_xscale('log')
cx.set_xticks([10, 100, 1000, 10000]); cx.set_xticklabels(['10', '100', '1k', '10k'])
cx.set_xlabel('clock age, kyr ($\\mu$ = 1.25e-8 /bp/gen, 29 y/gen)')
n_mid = int(((W.age_kyr_median_mu >= 200) & (W.age_kyr_median_mu <= 1500)).sum())
cx.set_title(f'Median-based ages of {n_mid} of {len(W)} clustered forms\nfall in 0.2-1.5 Myr; mean-based ages are higher', loc='left')
cx.legend(frameon=False, loc='lower right', fontsize=6)
cx.text(0.99, 0.12, 'dots: median-based age', transform=cx.transAxes, ha='right', fontsize=6)
# d
dx = fig.add_subplot(gs[1, 1])
sp = ['chimpanzee', 'bonobo', 'gorilla', 'sumatran_orangutan', 'bornean_orangutan']
n_tests = P.groupby('species').size().reindex(sp)
n_dup = D.groupby('species').size().reindex(sp).fillna(0)
x = np.arange(len(sp))
dx.bar(x, n_tests - n_dup, color='#d9d9d9', edgecolor='#707070', lw=0.4, label='probe not found across switch')
dx.bar(x, n_dup, bottom=n_tests - n_dup, color='#9970ab', label='found, but native probes align equally\n(shared duplicon)')
dx.set_xticks(x); dx.set_xticklabels(['chimpanzee', 'bonobo', 'gorilla', 'Sumatran\norangutan', 'Bornean\norangutan'], fontsize=5.5)
dx.set_ylabel(f'recurrent human junctions tested (n={P.junction_cluster.nunique()})')
dx.set_title('No recurrent human junction is present\nin a great-ape T2T genome', loc='left')
dx.legend(frameon=False, loc='upper left', fontsize=6, borderaxespad=0)
dx.set_ylim(0, 100); dx.set_yticks([0, 20, 40, 60, 70])
for a_, l in ((ax, 'a'), (bx, 'b'), (cx, 'c'), (dx, 'd')):
    a_.text(-0.26, 1.12, l, transform=a_.transAxes, fontweight='bold', fontsize=10, va='bottom')
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FIG, f'fig_origin.{ext}'), dpi=300, bbox_inches='tight')
print('saved')
