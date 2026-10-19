"""Figure: does every stretch of every molecule align better to the end it was assigned to?
From report/pitch_data/trace_windows.tsv.gz and trace_halves.tsv (website/tour_traces.py).

  figures/pitch/fig_alignment_contrast.(png|pdf)
    a  window-by-window contrast along single molecules of four claims (2 kb windows, each aligned on its own to the
       CHM13 donor end D and to the CHM13 end P of the block after the switch; contrast = d_P - d_D, differences per
       1,000 aligned letters; positive = closer to the donor end)
    b  every molecule side (before the switch interval; first 20 kb after it): median divergence to the end the side
       was assigned to against the median divergence to the other end, informative windows only

Usage: python report/pitch_alignment_contrast.py
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'report', 'pitch_data'); OUT = os.path.join(REPO, 'figures', 'pitch')
DONOR, RECIP, TELO, SWITCH, GREY = '#D55E00', '#0072B2', '#222222', '#F0E442', '#9A9A9A'
CAP, FLOOR, AFTER = 60.0, 0.5, 20_000
TRACE_CLAIMS = [('HG002_2_chr3_MATERNAL_R', 'HG002 3q: chr11p-like block, then chr3q'),
                ('HG01258__20q__SA-20q-01', 'HG01258 20q: chr3q-like block, then chr20q'),
                ('HG002_2_chr9_MATERNAL_L', 'HG002 9p: chr15q-like block, then chr9p'),
                ('HG002_2_chr17_MATERNAL_L', 'HG002 17p: indel-rich chr11p-like block (window test fails)')]
FAM = lambda p: 'PacBio HiFi' if p.startswith('HiFi') else ('Nanopore R9' if 'R9' in p else 'Nanopore R10')
MARK = {'PacBio HiFi': 'o', 'Nanopore R10': 's', 'Nanopore R9': '^'}


def style():
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6,
                         'ytick.labelsize': 6, 'legend.fontsize': 7, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'axes.titleweight': 'normal', 'axes.titlelocation': 'left', 'legend.frameon': False})


def main():
    style()
    W = pd.read_csv(os.path.join(DATA, 'trace_windows.tsv.gz'), sep='\t')
    H = pd.read_csv(os.path.join(DATA, 'trace_halves.tsv'), sep='\t'); H = H[H.after_span == AFTER].copy()
    J = {c['slug']: c for c in json.load(open(os.path.join(REPO, 'docs', 'evidence', 'tour', 'data', 'tour_data.json')))['claims']}
    fig = plt.figure(figsize=(7.2, 9.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.4, 3.6], hspace=0.28, left=0.11, right=0.87, top=0.95, bottom=0.07)
    ga = gs[0].subgridspec(len(TRACE_CLAIMS), 1, hspace=0.75)
    axs = [fig.add_subplot(ga[i]) for i in range(len(TRACE_CLAIMS))]
    for k, (ax, (slug, lab)) in enumerate(zip(axs, TRACE_CLAIMS)):
        c = J[slug]; g = W[W.slug == slug]
        ax.axhspan(0, CAP * 1.08, color=DONOR, alpha=0.07, lw=0); ax.axhspan(-CAP * 1.08, 0, color=RECIP, alpha=0.07, lw=0)
        for rid, m in g.groupby('read_id'):
            ax.axvspan(m.switch_lo.iloc[0] / 1e3, m.switch_hi.iloc[0] / 1e3, color=SWITCH, alpha=0.35, lw=0)
        for rid, m in g.groupby('read_id'):
            x = (m.win_start + 1000) / 1e3; ax.plot(x, m.delta, color=TELO, lw=0.6, alpha=0.75)
        ax.axhline(0, color=GREY, lw=0.5)
        ax.set_xlim(0, 100); ax.set_ylim(-CAP * 1.08, CAP * 1.08); ax.set_yticks([-60, 0, 60]); ax.set_yticklabels(['\u2264\u221260', '0', '\u226560'])
        ax.set_title(f'{lab} ({g.read_id.nunique()} molecules)', fontsize=7)
        ax.text(100.5, CAP * 0.55, f'closer to\n{c["D"]}', color=DONOR, fontsize=6, va='center', ha='left')
        ax.text(100.5, -CAP * 0.55, f'closer to\n{c["P"]}', color=RECIP, fontsize=6, va='center', ha='left')
        if k < len(TRACE_CLAIMS) - 1: ax.set_xticklabels([])
    axs[-1].set_xlabel('position along the molecule, from the telomere (kb)')
    fig.text(0.025, (axs[0].get_position().y1 + axs[-1].get_position().y0) / 2,
             'contrast: divergence to the other end minus divergence to the donor end\n'
             '(differences per 1,000 aligned letters; \u00b160 = capped, or aligns to one end only)', rotation=90,
             va='center', ha='center', fontsize=6)
    axs[0].text(-0.13, 1.42, 'a', transform=axs[0].transAxes, fontsize=10, fontweight='bold', va='bottom')
    axs[0].text(0, 1.42, 'Along each molecule: closer to the donor end before the switch (yellow), to the other end after it',
                transform=axs[0].transAxes, fontsize=8, va='bottom')
    # b
    gb = gs[1].subgridspec(1, 2, width_ratios=[0.66, 0.34], wspace=0.05)
    ax = fig.add_subplot(gb[0])
    H['fam'] = H.platform.map(FAM); H['p17'] = H.slug.str.contains('chr17_|17p', regex=True)
    H = H[H.n_informative > 0]
    H['own'] = np.where(H.side == 'before', H.med_d_D, H.med_d_P); H['oth'] = np.where(H.side == 'before', H.med_d_P, H.med_d_D)
    P = H.dropna(subset=['own', 'oth']).copy()
    Wi = W.copy(); Wi['cls'] = np.where(Wi.delta.isna() | ((Wi.d_D.fillna(1e9) > 100) & (Wi.d_P.fillna(1e9) > 100)), 'poor',
                                        np.where(Wi.delta.abs() < 2, 'shared', 'informative'))
    Wi = Wi[(Wi.cls == 'informative') & ((Wi.side == 'before') | ((Wi.side == 'after') & (Wi.win_end <= Wi.switch_hi + AFTER)))]
    Wi['own'] = np.where(Wi.side == 'before', Wi.d_D, Wi.d_P)
    own_med = Wi.groupby(['slug', 'read_id', 'side']).own.median()
    O = H[H.oth.isna()].copy(); O['own'] = [own_med.get((r.slug, r.read_id, r.side), np.nan) for r in O.itertuples()]
    O = O.dropna(subset=['own']); O['own'] = O.own.clip(lower=FLOOR)
    BAND = 900
    for v in ('own', 'oth'): P[v] = P[v].clip(lower=FLOOR)
    for (side, fam, p17), g in P.groupby(['side', 'fam', 'p17']):
        col = DONOR if side == 'before' else RECIP
        ax.scatter(g.own, g.oth, s=16, marker=MARK[fam], facecolor='none' if p17 else col, edgecolor=col, lw=0.7,
                   alpha=0.85, zorder=3)
    for (side, fam, p17), g in O.groupby(['side', 'fam', 'p17']):
        col = DONOR if side == 'before' else RECIP
        jit = BAND * np.exp(np.random.default_rng(1).uniform(-0.12, 0.12, len(g)))
        ax.scatter(g.own, jit, s=16, marker=MARK[fam], facecolor='none' if p17 else col, edgecolor=col, lw=0.7, alpha=0.85, zorder=3)
    ax.axhspan(BAND * 0.8, BAND * 1.25, color=GREY, alpha=0.12, lw=0)
    ax.text(0.75, BAND * 1.32, 'no alignment to the other end', fontsize=6, color=TELO, va='bottom')
    lim = (0.7, 400); ax.plot(lim, lim, color=GREY, lw=0.7, zorder=1)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(2.5, BAND * 1.9)
    tk = [1, 10, 100]; ax.set_xticks(tk); ax.set_yticks(tk); ax.set_xticklabels(['1', '10', '100']); ax.set_yticklabels(['1', '10', '100'])
    ax.set_xlabel('divergence to the end this side was assigned to\n(differences per 1,000 aligned letters, median)')
    ax.set_ylabel('divergence to the other end')
    ax.text(0.8, 300, 'above the line: closer to\nthe assigned end', fontsize=6, color=TELO, va='top')
    ax.text(25, 3.0, 'below: closer to the other end', fontsize=6, color=TELO, va='bottom')
    hv = H[~H.p17 & (H.verdict != 'not assessable')]; ab, nall = int((hv.verdict == 'as claimed').sum()), len(hv)
    ax.text(-0.17, 1.03, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
    ax.set_title(f'Molecule sides outside the 17p claims: {ab} of {nall} closer to the end they were assigned to', fontsize=8, loc='left')
    hs = [Line2D([], [], marker='o', ls='', color=DONOR, label='before the switch (assigned: donor end)'),
          Line2D([], [], marker='o', ls='', color=RECIP, label='after the switch (assigned: other end)')]
    hs += [Line2D([], [], marker=MARK[f], ls='', color=TELO, mfc='none' if False else TELO, label=f) for f in MARK]
    hs += [Line2D([], [], marker='o', ls='', color=TELO, mfc='none', label='17p claims (indel-rich donor block)')]
    lg = fig.add_subplot(gb[1]); lg.axis('off'); [s.set_visible(False) for s in lg.spines.values()]
    lg.legend(handles=hs, loc='upper left', fontsize=6, handletextpad=0.3, borderaxespad=0)
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUT, f'fig_alignment_contrast.{ext}'))
    s = P.assign(ok=P.oth > P.own).groupby(['p17', 'side']).ok.agg(['sum', 'size'])
    print(s.to_string()); print('one-end band', O.groupby(['p17', 'side']).size().to_dict(), 'not plotted', int(len(H) - len(P) - len(O)))
    print('headline', ab, 'of', nall)
    return fig


if __name__ == '__main__':
    main()
