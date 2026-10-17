"""Figures of the one-page pitch (report/PITCH.md), from report/pitch_data/ (report/pitch_data.py) and the catalogue.

  figures/pitch/fig_pitch_molecule.(png|pdf)   one 513 kb molecule: telomere, chromosome-11-like DNA, chromosome 3
  figures/pitch/fig_pitch_common.(png|pdf)     the most common foreign terminal alleles and which were seen in single molecules

Usage: python report/pitch_figures.py
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'report', 'pitch_data')
OUT = os.path.join(REPO, 'figures', 'pitch')
DONOR, RECIP, TELO, SWITCH, GREY = '#D55E00', '#0072B2', '#222222', '#F0E442', '#9A9A9A'
HERO_SLUG, HERO_READ = 'HG002_2_chr3_MATERNAL_R', '3788f60a-2f28-4b28-917b-f1ecae3645e3'
# claims whose distal side is free of long tandem repeats, so the 2 kb window comparison is interpretable
# (17p: tandem-repeat-rich 11p/17p block; NA20346 2q: switch at 7-8 kb leaves too few windows)
SLOPE_CLAIMS = [('HG002_2_chr3_MATERNAL_R', 'HG002 3q\n(11p-like end)'),
                ('HG002_2_chr9_MATERNAL_L', 'HG002 9p\n(15q-like end)'),
                ('HG004__HG002_2_chr9_MATERNAL_L', 'HG004 9p\n(HG002\'s mother)'),
                ('HG002_2_chr19_MATERNAL_R', 'HG002 19q\n(22q-like end)'),
                ('HG01258__20q__SA-20q-01', 'HG01258 20q\n(3q-like end)')]


def style():
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6,
                         'ytick.labelsize': 6, 'legend.fontsize': 7, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'axes.titleweight': 'normal', 'axes.titlelocation': 'left', 'legend.frameon': False})


def letter(ax, s, dx=-0.06, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='right')


def idx(dD, dR):
    """closeness index from differences per kb: +1 = matches only the donor end, -1 = only the recipient end"""
    return (dR - dD) / (dR + dD) if (dR + dD) > 0 else np.nan


def read_summary(P):
    rows = []
    for (slug, rid), g in P.groupby(['slug', 'read_id']):
        t, lo, hi = g.telo_len.iloc[0], g.switch_lo.iloc[0], g.switch_hi.iloc[0]
        b = g[(g.window_mid > t + 1500) & (g.window_mid < lo - 1000)]
        a = g[(g.window_mid > hi + 1000) & (g.window_mid < hi + 15000)]
        rows.append(dict(slug=slug, read_id=rid, platform=g.platform.iloc[0],
                         before=idx(b.diff_per_kb_D.median(), b.diff_per_kb_R.median()),
                         after=idx(a.diff_per_kb_D.median(), a.diff_per_kb_R.median())))
    return pd.DataFrame(rows)


def fig_molecule():
    P = pd.read_csv(os.path.join(DATA, 'profiles.tsv'), sep='\t')
    S = pd.read_csv(os.path.join(DATA, 'sites.tsv'), sep='\t')
    G = pd.read_csv(os.path.join(DATA, 'segments.tsv'), sep='\t')
    p = P[P.read_id == HERO_READ]
    s = S[S.read_id == HERO_READ]
    g = G[G.read_id == HERO_READ].set_index('segment')
    telo, lo, hi, L = p.telo_len.iloc[0] / 1e3, p.switch_lo.iloc[0] / 1e3, p.switch_hi.iloc[0] / 1e3, p.read_len.iloc[0] / 1e3
    XMAX, chr3_len = 80, 201_105_948
    d0, d1 = g.loc['distal', 'pos'] / 1e3, g.loc['distal', 'end'] / 1e3          # chr11 p-arm: tip frame = forward strand
    r_off = (chr3_len - g.loc['proximal', 'end']) / 1e3                         # chr3 q-arm: kb from the CHM13 tip
    pre, post = s[s.read_pos < lo * 1e3], s[s.read_pos > hi * 1e3]
    nD_pre, n_pre = int((pre.call == 'D').sum()), len(pre)
    nP_post, n_post = int((post.call == 'P').sum()), len(post)

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(4, 5, height_ratios=[1.75, 0.5, 0.95, 1.45], hspace=0.75, wspace=0.08,
                          left=0.10, right=0.985, top=0.955, bottom=0.06)
    # ---------------- a: the molecule against the two reference ends, to scale
    ax = fig.add_subplot(gs[0, :])
    y11, yM, y3, h = 2.3, 1.15, 0.0, 0.30
    ax.add_patch(Rectangle((0, y11 - h / 2), XMAX, h, color=DONOR, alpha=0.85, lw=0))
    ax.add_patch(Rectangle((0, y3 - h / 2), r_off, h, color=RECIP, alpha=0.22, lw=0))
    ax.add_patch(Rectangle((r_off, y3 - h / 2), XMAX - r_off, h, color=RECIP, alpha=0.85, lw=0))
    for y in (y11, y3): ax.add_patch(Rectangle((0, y - h / 2), 0.6, h, color=TELO, lw=0))
    ax.add_patch(Rectangle((0, yM - h / 2), telo, h, color=TELO, lw=0))
    ax.add_patch(Rectangle((telo, yM - h / 2), lo - telo, h, color=DONOR, lw=0))
    ax.add_patch(Rectangle((lo, yM - h / 2), hi - lo, h, facecolor=SWITCH, edgecolor='#8a7a00', lw=0.6))
    ax.add_patch(Rectangle((hi, yM - h / 2), XMAX - hi, h, color=RECIP, lw=0))
    ax.add_patch(Polygon([(telo, yM + h / 2), (lo, yM + h / 2), (d1, y11 - h / 2), (d0, y11 - h / 2)],
                         color=DONOR, alpha=0.15, lw=0))
    ax.add_patch(Polygon([(hi, yM - h / 2), (XMAX, yM - h / 2), (XMAX, y3 + h / 2), (r_off, y3 + h / 2)],
                         color=RECIP, alpha=0.15, lw=0))
    lab = dict(ha='left', va='bottom', fontsize=7)
    ax.text(0, y11 + h / 2 + 0.06, 'Reference genome (CHM13): the end of chromosome 11', color=DONOR, **lab)
    ax.text(0, yM + h / 2 + 0.06, f'One DNA molecule from HG002 (mother\'s chromosome 3), {L:.0f} kb long',
            color=TELO, fontweight='bold', **lab)
    ax.text(0, y3 + h / 2 + 0.06, 'Reference genome (CHM13): the end of chromosome 3', color=RECIP, **lab)
    ax.text((telo + lo) / 2, yM, 'chromosome-11-like', ha='center', va='center', color='white', fontsize=7)
    ax.text((hi + XMAX) / 2, yM, f'chromosome 3 (continues {L - XMAX:.0f} kb further)', ha='center', va='center',
            color='white', fontsize=7)
    ax.text(r_off / 2, y3, 'not in this molecule', ha='center', va='center', color=RECIP, fontsize=6)
    ax.annotate('telomere', (telo / 2, yM - h / 2), (2.2, yM - 0.48), fontsize=6, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.6, color=TELO))
    ax.annotate('switch', ((lo + hi) / 2, yM - h / 2), ((lo + hi) / 2 - 6, yM - 0.48), fontsize=6, ha='center',
                va='center', arrowprops=dict(arrowstyle='-', lw=0.6, color='#8a7a00'))
    ax.text(d1 + 1.5, (yM + y11) / 2 + 0.05, f'this part maps uniquely to chromosome 11 (MAPQ 60)', fontsize=6,
            color=DONOR, va='center')
    ax.text(XMAX - 0.5, (yM + y3) / 2 - 0.08, f'this part maps uniquely to chromosome 3 (MAPQ 60, '
            f'{g.loc["proximal", "aligned_bp"] / 1e3:.0f} kb aligned)', fontsize=6, color=RECIP, va='center', ha='right')
    ax.set_xlim(0, XMAX); ax.set_ylim(-0.3, 2.75); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('distance from the telomeric end (kb)')
    ax.set_title('Cut this molecule at the switch: one part maps only to chromosome 11, the rest only to chromosome 3')
    letter(ax, 'a', dx=-0.015)

    # ---------------- b: every diagnostic letter along the molecule
    axt = fig.add_subplot(gs[1, :])
    for call, col, y in (('D', DONOR, 1), ('P', RECIP, 0)):
        axt.vlines(s[s.call == call].read_pos / 1e3, y - 0.36, y + 0.36, color=col, lw=0.7)
    axt.axvspan(lo, hi, color=SWITCH, alpha=0.6, lw=0)
    axt.set_xlim(0, XMAX); axt.set_ylim(-0.6, 1.6); axt.set_yticks([0, 1])
    axt.set_yticklabels(['chr3 letter', 'chr11 letter']); axt.tick_params(axis='y', length=0)
    axt.text(lo / 2, 1.5, f'{nD_pre} of {n_pre} sites: chr11 letter', color=DONOR, fontsize=6, ha='center', va='bottom')
    axt.text(post.read_pos.max() / 1e3 + 1.5, 0, f'{nP_post} of {n_post} sites: chr3 letter', color=RECIP, fontsize=6,
             ha='left', va='center')
    axt.set_xlabel('distance from the telomeric end (kb)', labelpad=1)
    axt.set_title('At every site where the two chromosomes differ, which letter does the molecule carry?', pad=10)
    letter(axt, 'b', dx=-0.015, dy=1.22)

    # ---------------- c: the letters themselves
    axl = fig.add_subplot(gs[2, :])
    k = 6
    pick = pd.concat([pre.iloc[np.linspace(0, len(pre) - 1, k).round().astype(int)],
                      post.iloc[np.linspace(0, len(post) - 1, k).round().astype(int)]])
    for i, r in enumerate(pick.itertuples()):
        xx = i + (0.8 if i >= k else 0)
        axl.text(xx, 2, r.donor_letter, ha='center', va='center', color=DONOR, fontsize=8, family='monospace')
        axl.text(xx, 0, r.recipient_letter, ha='center', va='center', color=RECIP, fontsize=8, family='monospace')
        mc = DONOR if r.call == 'D' else RECIP if r.call == 'P' else GREY
        axl.add_patch(Rectangle((xx - 0.36, 1 - 0.40), 0.72, 0.80, color=mc, lw=0))
        axl.text(xx, 1, r.read_letter, ha='center', va='center', color='white', fontsize=8, family='monospace',
                 fontweight='bold')
        axl.text(xx, -0.9, f'{r.read_pos / 1e3:.1f} kb', ha='center', va='center', fontsize=6, color='#555555')
    axl.axvline(k - 0.1, color='#8a7a00', lw=0.8, ls=(0, (2, 2)), ymin=0.08, ymax=0.95)
    axl.text(k - 0.1, 2.75, 'switch', ha='center', va='bottom', fontsize=6, color='#8a7a00')
    axl.text(-0.7, 2, 'chromosome 11 says', ha='right', va='center', color=DONOR, fontsize=7)
    axl.text(-0.7, 1, 'the molecule says', ha='right', va='center', color=TELO, fontsize=7, fontweight='bold')
    axl.text(-0.7, 0, 'chromosome 3 says', ha='right', va='center', color=RECIP, fontsize=7)
    axl.set_xlim(-3.6, 2 * k + 0.4); axl.set_ylim(-1.3, 3.2); axl.axis('off')
    axl.set_title(f'The letters, at {k} sites evenly spaced on each side of the switch', x=0.0)
    letter(axl, 'c', dx=-0.015, dy=0.96)

    # ---------------- d: the same flip in every interpretable showcase molecule
    R = read_summary(P)
    axes = []
    for j, (slug, lab_) in enumerate(SLOPE_CLAIMS):
        axs = fig.add_subplot(gs[3, j], sharey=axes[0] if axes else None); axes.append(axs)
        rr = R[R.slug == slug]
        for r in rr.itertuples():
            flip = r.before > 0 and r.after < 0
            axs.plot([0, 1], [r.before, r.after], '-o', color=TELO if flip else GREY, lw=0.9, ms=2.6, alpha=0.9)
        axs.axhline(0, color=GREY, lw=0.6)
        axs.set_xlim(-0.35, 1.35); axs.set_xticks([0, 1]); axs.set_xticklabels(['before\nswitch', 'after\nswitch'])
        n_flip = int(((rr.before > 0) & (rr.after < 0)).sum())
        axs.set_title(f'{lab_}\n{n_flip} of {len(rr)} molecules flip', fontsize=7, loc='center')
        if j: axs.tick_params(axis='y', labelleft=False)
    axes[0].set_ylabel('closer to donor (+1)\nor recipient (−1)')
    axes[0].set_ylim(-1.08, 1.08); axes[0].set_yticks([-1, 0, 1])
    letter(axes[0], 'd', dx=-0.42, dy=1.3)
    tot = sum(int(((R[R.slug == sl].before > 0) & (R[R.slug == sl].after < 0)).sum()) for sl, _ in SLOPE_CLAIMS)
    ntot = int(R.slug.isin([sl for sl, _ in SLOPE_CLAIMS]).sum())
    bb = axes[0].get_position()
    fig.text(bb.x0, bb.y1 + 0.075, f'Not a one-off: {tot} of {ntot} molecules (5 claims, 3 people) turn from donor-like to '
             f'recipient-like at their switch', fontsize=8, ha='left', va='bottom')
    return fig, dict(telo=telo, lo=lo, hi=hi, L=L, d0=d0, d1=d1, r_off=r_off, nD_pre=nD_pre, n_pre=n_pre,
                     nP_post=nP_post, n_post=n_post, summary=R, pick=pick)


def fig_common():
    A = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'alleles_structural.tsv'), sep='\t')
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False)
    F = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'paralog_families.tsv'), sep='\t')
    fam = dict(zip(F.family, F.members.str.replace(',', '/')))
    fam['PF01'] = '4p/acrocentric'                    # PF01 = 4p, 14p, 15p, 21p, 22p
    d = C[C.definitive_verdict.astype(str).str.startswith('definitive')]
    dd = d.groupby('allele_id').individual.nunique()
    top = A[(A.nativity_class == 'foreign') & A.form.notna()].sort_values('r2_freq', ascending=False).head(12).copy()
    top['people_mol'] = top.allele_id.map(dd).fillna(0).astype(int)

    def nice(form):
        arm, src_ = form.split('>')
        return f'{arm} end carries {fam.get(src_, src_)}-like DNA'
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    plt.subplots_adjust(left=0.37, right=0.70, top=0.84, bottom=0.13)
    y = np.arange(len(top))[::-1]
    for yy, r in zip(y, top.itertuples()):
        seen = r.people_mol > 0
        ax.hlines(yy, 0, r.r2_freq * 100, color=TELO if seen else GREY, lw=1.4)
        ax.plot(r.r2_freq * 100, yy, 'o', color=TELO if seen else 'white', mec=TELO if seen else GREY, ms=5)
        ax.text(-1.5, yy, f'{nice(r.form)} ({r.allele_id})', ha='right', va='center', fontsize=7)
        ax.text(r.r2_freq * 100 + 1.2, yy, f'{r.n_samples} people', va='center', fontsize=6, color='#555555')
        if seen:
            ax.text(64, yy, f'single-molecule evidence in {r.people_mol} {"person" if r.people_mol == 1 else "people"}',
                    va='center', fontsize=6, color=TELO, clip_on=False)
    ax.set_xlim(0, 60); ax.set_ylim(-0.7, len(top) - 0.3); ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('share of assembled ends of that chromosome arm (%)')
    fig.text(0.02, 0.95, 'These are common: the 12 most frequent chromosome ends that carry another chromosome\'s DNA',
             fontsize=8, ha='left', va='top')
    fig.text(0.02, 0.90, '\u25cf definitive single-molecule evidence in at least one person    '
             '\u25cb not yet shown in single molecules (not assessable or not confirmed)', fontsize=6, ha='left',
             va='top', color='#555555')
    return fig, top


def main():
    style(); os.makedirs(OUT, exist_ok=True)
    f1, info = fig_molecule()
    for ext in ('png', 'pdf'): f1.savefig(os.path.join(OUT, f'fig_pitch_molecule.{ext}'))
    f2, top = fig_common()
    for ext in ('png', 'pdf'): f2.savefig(os.path.join(OUT, f'fig_pitch_common.{ext}'))
    print({k: v for k, v in info.items() if k != 'summary'})
    print(top[['allele_id', 'form', 'r2_freq', 'n_samples', 'people_mol']].to_string(index=False))


if __name__ == '__main__':
    main()
