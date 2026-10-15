"""Figure for the junctions track (H2): one-way transfer at homologous breakpoints.
Inputs : results/reciprocal_summary.json, results/junction_features.tsv, results/junction_null_summary.tsv,
         results/exemplars.tsv, results/exemplar_<id>_track.tsv, ATLAS_INPUT_DIR/results/subtel_foreign_tips.tsv
Outputs: figures/junctions_mechanism.png (300 dpi) and .pdf, figures/junctions_mechanism_caption.md
Usage  : python scripts/figure.py
"""
import os, sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 8,
                     'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
DONOR, RECIP, NULLC, OBS = '#D55E00', '#0072B2', '#9E9E9E', '#333333'
STATUS = [('native', 'native donor end', '#4D4D4D'), ('replaced_other', 'donor end foreign to a 3rd arm', '#A6A6A6'),
          ('has_recipient_window', 'recipient-like window inside', '#CC79A7'), ('reciprocal', 'reciprocal (donor carries recipient end)', '#E69F00')]

def letter(ax, s): ax.text(-0.13, 1.09, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')

def main():
    S = json.load(open(os.path.join(J.RES, 'reciprocal_summary.json')))
    X = pd.read_csv(os.path.join(J.RES, 'junction_features.tsv'), sep='\t')
    NS = pd.read_csv(os.path.join(J.RES, 'junction_null_summary.tsv'), sep='\t')
    E = pd.read_csv(os.path.join(J.RES, 'exemplars.tsv'), sep='\t')
    T = pd.read_csv(J.inp('results/subtel_foreign_tips.tsv'), sep='\t', low_memory=False, usecols=['t20_state', 'switch_res', 'sw_own_dist'])
    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.8], hspace=0.62, wspace=0.42, left=0.1, right=0.94, top=0.95, bottom=0.07)
    # a: donor-end status
    ax = fig.add_subplot(gs[0, 0]); P = S['pooled']
    for i, grp in enumerate(['carriers', 'noncarriers']):
        left = 0
        for k, lab, col in STATUS:
            f = P[grp][k] / P[grp]['assessed']; ax.barh(i, f, left=left, color=col, height=0.6, edgecolor='white', linewidth=0.4); left += f
    ax.set_yticks([0, 1]); ax.set_yticklabels([f"carriers\n(n = {P['carriers']['assessed']:,})", f"non-carriers\n(n = {P['noncarriers']['assessed']:,})"])
    ax.set_ylim(2.5, -0.5); ax.set_xlim(0, 1); ax.set_xlabel('fraction of same-haplotype donor ends')
    ax.set_title('Donor ends of carriers match non-carriers', loc='left')
    ax.text(1.0, 1.48, f"χ² p = {P['chi2_donor_status']['p']:.2f}", ha='right', va='center', fontsize=7)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for k, l, c in STATUS]
    ax.legend(h, [l for k, l, c in STATUS], loc='lower left', bbox_to_anchor=(0.0, 0.0), frameon=False, ncol=2, fontsize=6, handlelength=1, borderaxespad=0.2, columnspacing=0.8)
    letter(ax, 'a')
    # b: replaced own terminal sequence
    ax = fig.add_subplot(gs[1, 0])
    v = T[(T.t20_state == 'foreign') & (T.switch_res == 'base')].sw_own_dist.dropna().values / 1000
    bins = np.logspace(np.log10(0.5), np.log10(300), 36)
    ax.hist(v, bins=bins, color=RECIP, alpha=0.85, label=f'all base-resolved tips (n = {len(v):,})')
    ax.hist(X.sw_own_dist.values / 1000, bins=bins, color='white', edgecolor=OBS, linewidth=0.8, histtype='step', label=f'one per junction cluster (n = {len(X)})')
    ax.set_xscale('log'); ax.set_xticks([1, 3, 10, 30, 100, 300]); ax.set_xticklabels(['1', '3', '10', '30', '100', '300'])
    ax.set_xlabel('recipient sequence replaced (kb from telomere)'); ax.set_ylabel('tips')
    ax.axvline(np.median(v), ymax=0.76, color=OBS, lw=0.6, ls='--'); ax.text(np.median(v) * 0.93, ax.get_ylim()[1] * 0.97, f'median {np.median(v):.1f} kb', fontsize=6, ha='right')
    ax.set_title('Terminal 1-200 kb of the arm is replaced', loc='left'); ax.legend(frameon=False, loc='upper right', fontsize=6)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3); letter(ax, 'b')
    # c: parental identity at the breakpoint
    ax = fig.add_subplot(gs[0, 1])
    cls = [('0', lambda m: m == 0), ('1-4', lambda m: (m >= 1) & (m <= 4)), ('5-25', lambda m: (m >= 5) & (m <= 25)), ('26-200', lambda m: (m >= 26) & (m <= 200)), ('>200', lambda m: m > 200)]
    obs = [f(X.mh_parental).mean() for _, f in cls]
    nullf = json.load(open(os.path.join(J.RES, 'junction_summary.json')))['mh_classes_null_frac']
    nk = ['blunt_0', 'micro_1_4', 'micro_5_25', 'homology_26_200', 'homology_gt200']; nul = [nullf.get(k, 0.0) for k in nk]
    xx = np.arange(len(cls)); ax.bar(xx - 0.2, obs, 0.38, color=OBS, label=f'junctions (n = {len(X)})'); ax.bar(xx + 0.2, nul, 0.38, color=NULLC, label='random parental pairs')
    for xi, n_ in zip(xx + 0.2, nul):
        if n_ < 0.005: ax.plot([xi], [0.004], marker='_', color=NULLC, ms=6, mew=1.2)
    ax.set_xticks(xx); ax.set_xticklabels([c for c, _ in cls]); ax.set_xlabel('identical parental sequence at breakpoint (bp)'); ax.set_ylabel('fraction')
    ax.set_title('Breakpoints lie in donor-recipient homology', loc='left'); ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.36, 1.0), fontsize=6)
    ax.set_ylim(0, max(max(obs), max(nul)) * 1.15); letter(ax, 'c')
    # d: ITS
    ax = fig.add_subplot(gs[1, 1])
    rows = [('all', 'its_run_c_any', 'carrier, all'), ('far', 'its_run_c_any', 'carrier, ≥ 5 kb'), ('all', 'its_run_P_any', 'parents, all'), ('far', 'its_run_P_any', 'parents, ≥ 5 kb')]
    for i, (sub, k, lab) in enumerate(rows):
        r = NS[(NS.subset == sub) & (NS.feature == k)].iloc[0]
        ax.barh(i - 0.18, r.obs_mean, 0.34, color=OBS); ax.barh(i + 0.18, r.null_mean, 0.34, color=NULLC)
        ax.text(max(r.obs_mean, r.null_mean) + 0.005, i, f"p = {r.perm_p_greater:.2g}\nn = {int(r.junctions)}", va='center', fontsize=6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([l for _, _, l in rows]); ax.invert_yaxis(); ax.set_xlim(0, 0.24)
    ax.set_xlabel('fraction with ≥ 3 tandem TTAGGG-like\ncopies within ±1 kb')
    ax.set_title('ITS excess is confined to the first 5 kb', loc='left')
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=OBS), plt.Rectangle((0, 0), 1, 1, color=NULLC)], ['junctions', 'random positions'], frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.04), fontsize=6)
    letter(ax, 'd')
    # e: exemplar identity tracks
    sub = gs[2, :].subgridspec(1, 3, wspace=0.12)
    names = {'E1_20q-3q': '20q>3q, HPRC R2 (121 samples)', 'E2_3q-11p': '3q>11p, HG002 maternal chr3', 'E3_16p-9q': '16p>9q, HG002 paternal chr16'}
    for i, r in enumerate(E.itertuples()):
        ax = fig.add_subplot(sub[i]); tk = pd.read_csv(os.path.join(J.RES, f'exemplar_{r.exemplar}_track.tsv'), sep='\t')
        ax.plot(tk.rel / 1000, tk.ident_donor_100bp, color=DONOR, lw=0.8); ax.plot(tk.rel / 1000, tk.ident_recipient_100bp, color=RECIP, lw=0.8)
        ax.axvline(0, color=OBS, lw=0.5, ls=':'); ax.set_ylim(0.15, 1.04); ax.set_xlim(-5, 5); ax.set_yticks([0.5, 1.0])
        ds = tk[tk.site == 'donor_specific'].rel / 1000; rs = tk[tk.site == 'recipient_specific'].rel / 1000
        ax.vlines(ds, 0.17, 0.23, color=DONOR, lw=0.3); ax.vlines(rs, 0.25, 0.31, color=RECIP, lw=0.3)
        ax.set_title(f"{names[r.exemplar]}\n{int(r.mh_parental)} bp parental identity at switch", loc='left', fontsize=7)
        ax.text(-4.8, 0.36, f'donor {r.donor_end[3:]}', color=DONOR, fontsize=6); ax.text(4.8, 0.36, f'recipient {r.recipient_end[3:]}', color=RECIP, fontsize=6, ha='right')
        if i == 2: ax.text(5.2, 0.24, 'diagnostic\nsites', fontsize=6, va='center', ha='left')
        if i == 0: ax.set_ylabel('identity to CHM13\n(100 bp window)'); letter(ax, 'e')
        else: ax.set_yticklabels([])
        if i == 1: ax.set_xlabel('position relative to junction (kb; telomere to the left)')
    out = os.path.join(J.FIG, 'junctions_mechanism')
    fig.savefig(out + '.png', dpi=300); fig.savefig(out + '.pdf')
    r = fig.canvas.get_renderer(); import matplotlib as mpl
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    print('text overlaps:', ov)
    cap = f"""**Foreign terminal blocks arise by one-way transfer joined within donor-recipient homology.**
(a) Status of the donor arm's telomeric end in the same haplotype assembly, pooled over the {S['forms']} forms with >= 10 carrier tips (HPRC Release 2):
carriers of form R>D (n = {P['carriers']['assessed']:,} assessed donor ends; {P['carriers']['absent']:,} carriers had no assembled named donor end) versus haplotypes whose
R end is assessed but does not carry R>D (n = {P['noncarriers']['assessed']:,}); chi-square test on the four categories.
(b) Distance from the CHM13 recipient telomere at base-resolved switches (= own terminal sequence replaced); filled: all base-resolved tips, outline: one representative per junction cluster.
(c) Length of sequence identical between the two CHM13 parental ends at the breakpoint (forward + backward extension from the junction coordinates) for {len(X)} junction-cluster
representatives versus {NS[NS.feature == 'mh_parental'].null_draws.iloc[0]} random position pairs per junction in the same two parental ends (5-200 kb from the telomere).
(d) Interstitial telomere-like arrays (>= 3 in-phase TTAGGG/CCCTAA copies with <= 1 mismatch each) within +-1 kb of the junction in the carrier and in either parent, versus random
positions in the same sequences; '>= 5 kb' = junctions >= 5 kb from the telomere in carrier and both parents; one-sided permutation p.
(e) 100 bp rolling identity of three exemplar carriers to the donor (orange) and recipient (blue) CHM13 ends around the junction (single best
alignment on each side); ticks below: diagnostic sites where the carrier matches only the donor (orange) or only the recipient (blue, includes bases where the donor does not align).
Read support: E2 and E3 are crossed by {int(E.ont_ul_dorado_crossing_10kb.iloc[1])} and {int(E.ont_ul_dorado_crossing_10kb.iloc[2])} ONT UL and {int(E.revio_3cells_crossing_10kb.iloc[1])} and {int(E.revio_3cells_crossing_10kb.iloc[2])} PacBio Revio primary alignments (switch +-5 kb, none clipped at the switch).
"""
    open(out + '_caption.md', 'w').write(cap)

if __name__ == '__main__':
    main()
