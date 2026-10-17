#!/usr/bin/env python3
"""Figures of the read-evidence resource (evidence/figures/; site copies under docs/evidence/fig/).

  extract   per-unit calls and anchor alignments of every v1-conclusive read, from the pipeline outputs
            (work/evidence/{hg002,trio}/stage.pkl, work/evidence/population/units.json.gz + genome.sam)
            -> work/evidence/release_build/units_conclusive.json.gz
  showcase  (a) read-through figure per showcase exemplar + combined panel figure
  summary   (b) claims by cohort and tier, control specificity, dataset independence, trio transmission
  claims    (c) per-claim small multiples of all definitive reads (one PNG + SVG per definitive claim)
Usage: python evidence/scripts/evidence_figures.py extract|showcase|summary|claims|all
"""
import os, sys, json, gzip, pickle, collections, argparse
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_lib as L  # noqa: E402

W = os.path.join(REPO, 'work', 'evidence'); B = os.path.join(W, 'release_build')
FIG = os.path.join(REPO, 'evidence', 'figures')
UNITS = os.path.join(B, 'units_conclusive.json.gz')
CD, CP, CO, CA, CT = '#D55E00', '#0072B2', '#9a9a9a', '#009E73', '#3a3a3a'     # donor, proximal/recipient, other, anchor, telomere
CI_BAND = '#f2d675'


def style():
    import matplotlib as mpl
    mpl.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'DejaVu Sans',
                         'axes.titlelocation': 'left', 'savefig.dpi': 300, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                         'legend.frameon': False})


# --------------------------------------------------------------------------------------------------------------- extract
def _best_anchor(alns, R, th_min_mapq=60, th_min_bp=10000, th_max=2_000_000):
    """the anchor alignment chosen by readevidence.finalize (same rule)"""
    chrom, arm = R[:-1], R[-1]; want = '+' if arm == 'p' else '-'; best = None
    for a in alns or []:
        if a['chrom'] != chrom or a['strand'] != want or a['mapq'] < th_min_mapq or a['aln_len'] < th_min_bp: continue
        dist = a['ref_start'] if arm == 'p' else a['chrom_len'] - a['ref_end']
        if dist > th_max: continue
        if best is None or a['aln_len'] > best['aln_len']: best = a
    return best


def _sam_alns(sam, names):
    import pysam
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(sam) as f:
        L_ = dict(zip(f.references, f.lengths))
        for r in f:
            if r.is_unmapped or r.query_name not in names: continue
            ql = r.infer_read_length(); qs, qe = r.query_alignment_start, r.query_alignment_end
            if r.is_reverse: qs, qe = ql - qe, ql - qs
            out[r.query_name].append(dict(chrom=r.reference_name, chrom_len=L_[r.reference_name], strand='-' if r.is_reverse else '+',
                                          ref_start=r.reference_start, ref_end=r.reference_end, q_start=qs, q_end=qe,
                                          mapq=r.mapping_quality, aln_len=qe - qs))
    return out


def _pack_units(u):
    if isinstance(u, dict): return [u['read_pos'], u['allele'][0], u['kind'][0], u['d_pos'], u['p_pos']]
    return list(u)


def extract(a):
    import pandas as pd
    Rr = pd.read_csv(os.path.join(REPO, 'evidence', 'reads.tsv.gz'), sep='\t', low_memory=False)
    conc = Rr[(Rr.v1_verdict == 'conclusive') & ~Rr.claim_id.str.startswith('control|')]
    want = {(c, d, r) for c, d, r in zip(conc.claim_id, conc.dataset.astype(str), conc.read_id)}
    out = {}
    for name, prefix in (('hg002', ''), ('trio', 'trio')):
        S = pickle.load(open(os.path.join(W, name, 'stage.pkl'), 'rb')); Gm = S['Gm']
        for j, rows, _ in S['res']:
            for r, u in rows:
                cid = r['claim_id'] if name == 'hg002' else f"{r.get('sample') or j.get('sample')}|{r['claim_id']}"
                k = (cid, str(r.get('dataset')), r['read'])
                if k not in want or '|'.join(k) in out: continue
                if name == 'hg002' and r.get('role') != 'carrier': continue
                R = conc[(conc.claim_id == cid) & (conc.read_id == r['read'])].iloc[0]
                best = _best_anchor(Gm.get(r['_a']), _recip(cid, R))
                out['|'.join(k)] = dict(units=[_pack_units(x) for x in u], read_len=r['read_len'], s0=r['s0_read'],
                                        telo_off=r['telo_offset'], telo_len=r['telo_len'], anchor_from=r['anchor_from'],
                                        anchor_read=[r['anchor_from'] + best['q_start'], r['anchor_from'] + best['q_end']] if best else None)
        print(name, len(out), file=sys.stderr)
    # population
    P = pd.read_csv(os.path.join(W, 'population', 'reads_full.tsv'), sep='\t', low_memory=False,
                    usecols=['claim_id', 'role', 'sample', 'read', 'job', 'verdict', 'anchor_from', 'R', 'read_len', 's0_read', 'telo_offset', 'telo_len'])
    P = P[(P.role == 'carrier') & (P.verdict == 'conclusive')]
    U = json.load(gzip.open(os.path.join(W, 'population', 'units.json.gz'), 'rt'))
    keys = {f"A|{s}|{n}|{int(af)}": (c, s, n) for c, s, n, af in zip(P.claim_id, P['sample'], P.read, P.anchor_from)}
    G = _sam_alns(os.path.join(W, 'population', 'genome.sam'), set(keys))
    for r in P.itertuples():
        k = (r.claim_id, f'1KG_ONT_{r.sample}', r.read)
        if k not in want: continue
        best = _best_anchor(G.get(f"A|{r.sample}|{r.read}|{int(r.anchor_from)}"), r.R)
        out['|'.join(k)] = dict(units=U[f'{r.job}|{r.read}'], read_len=int(r.read_len), s0=int(r.s0_read), telo_off=r.telo_offset,
                                telo_len=int(r.telo_len), anchor_from=int(r.anchor_from),
                                anchor_read=[int(r.anchor_from) + best['q_start'], int(r.anchor_from) + best['q_end']] if best else None)
    miss = [k for k in want if '|'.join(k) not in out]
    print(f'units for {len(out)} of {len(want)} conclusive rows; missing {len(miss)}', file=sys.stderr)
    with gzip.open(UNITS, 'wt') as f: json.dump(out, f, default=lambda o: o.item() if hasattr(o, 'item') else str(o))


_CL = {}


def _recip(cid, R):
    import pandas as pd
    if not _CL:
        C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t')
        _CL.update(dict(zip(C.claim_id, C.recipient_R)))
    return _CL[cid]


# --------------------------------------------------------------------------------------------------------------- helpers
def _load():
    import pandas as pd
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False)
    Rr = pd.read_csv(os.path.join(REPO, 'evidence', 'reads.tsv.gz'), sep='\t', low_memory=False)
    U = json.load(gzip.open(UNITS, 'rt'))
    fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(os.path.join(os.environ.get('ATLAS_INPUT_DIR', ''), 'ref', 'chm13v2.0.fa.gz.fai'))}
    return C, Rr, U, fai


def _label(c):
    if c['cohort'] == 'HG002_Q100': return f"HG002 {c['arm']} {c['haplotype']} ({c['allele_id']})"
    if c['cohort'] == 'HG002_trio': return f"{c['individual']} {c['arm']} ({c['haplotype'].split(' (')[0].replace('_', ' ')}, {c['allele_id']})"
    return f"{c['individual']} {c['arm']} {c['allele_id']}"


def _end_label(e):
    return f"CHM13 {e[3:]}"


def _units_split(u):
    import numpy as np
    a = np.array([[x[0], {'D': 1, 'P': 0}.get(x[1], 2), 1 if x[2] == 's' else 0, x[3], x[4]] for x in u], float).reshape(-1, 5)
    return a


def _fmt_coord(chrom, x):
    return f"{chrom}:{int(x):,}"


# --------------------------------------------------------------------------------------------------------------- (a)
def read_through(ax, c, r, u, fai, compact=False):
    """one read, telomere -> anchor: units coloured donor vs proximal, switch interval, anchor, CHM13 coordinates below"""
    import numpy as np
    from matplotlib.patches import Rectangle, Polygon
    A = _units_split(u['units']); kb = 1e-3
    L_ = u['read_len']; s0 = u['s0']; t0 = u['telo_off'] or 0; tl = u['telo_len'] or 0
    lo, hi = r['switch_ci_read_lo'], r['switch_ci_read_hi']; sw = r['switch_read']
    an = u['anchor_read']
    full = (an[1] if an else L_)
    last_unit = float(A[:, 0].max()) if len(A) else hi
    xmax = min(full, max(last_unit, hi) + 15000) * kb          # crop: units and switch plus 15 kb of anchor
    yR, yS, yB, yRef1, yRef2 = 3.0, 2.55, 2.25, 1.15, 0.35
    ax.axvspan(lo * kb, hi * kb, ymin=0.0, ymax=1, color=CI_BAND, alpha=0.55, lw=0, zorder=0)
    ax.plot([0, L_ * kb], [yR, yR], color='#bdbdbd', lw=1.2, solid_capstyle='butt', zorder=1)
    ax.add_patch(Rectangle((t0 * kb, yR - 0.09), tl * kb, 0.18, color=CT, lw=0, zorder=2))
    # units: SNV row and private-bin row
    for kind, y in ((1, yS), (0, yB)):
        for lab, col in ((2, CO), (0, CP), (1, CD)):
            m = (A[:, 2] == kind) & (A[:, 1] == lab)
            if m.any(): ax.vlines(A[m, 0] * kb, y - 0.13, y + 0.13, color=col, lw=0.6 if kind else 1.6, zorder=3, rasterized=True)
    ax.axvline(sw * kb, color='#555555', lw=0.7, ls='--', zorder=4)
    if an:
        ax.add_patch(Rectangle((an[0] * kb, yR - 0.06), (an[1] - an[0]) * kb, 0.12, color=CA, lw=0, zorder=2))
    # CHM13 reference bars in tip-frame coordinates: donor D (units distal to the switch) and P (proximal units + anchor)
    D, P, R = c['donor_D'], c['proximal_P'], c['recipient_R']
    dmask = (A[:, 0] < lo) & (A[:, 1] == 1) & (A[:, 3] >= 0)
    pmask = (A[:, 0] > hi) & (A[:, 1] == 0) & (A[:, 4] >= 0)
    rows = []
    if dmask.any(): rows.append((D, CD, yRef1, A[dmask, 0].min(), A[dmask, 0].max(), A[dmask, 3].min(), A[dmask, 3].max()))
    pa = A[pmask]
    if len(pa):
        q0, q1, x0, x1 = pa[:, 0].min(), pa[:, 0].max(), pa[:, 4].min(), pa[:, 4].max()
        if an and P == R and isinstance(r.get('anchor_chm13'), str) and ':' in r['anchor_chm13']:
            cs, ce = (int(v) for v in r['anchor_chm13'].split(':')[1].split('-'))
            te = ce if R.endswith('p') else fai[R[:-1]] - cs + 1          # tip-frame end of the anchor on R
            q1, x1 = an[1], max(x1, float(te))
        rows.append((P, CP, yRef2, q0, q1, x0, x1))
    for end, col, y, q0, q1, x0, x1 in rows:
        # the bar is drawn in CHM13 tip-frame coordinates when the block is near the CHM13 telomere; a block internal to the
        # CHM13 end (e.g. 192-227 kb from the 11p telomere) is drawn shifted under the read and labelled with its true position
        sh = 0.0 if abs(q0 - x0) <= 10000 else q0 - x0
        b0, b1 = (x0 + sh) * kb, (x1 + sh) * kb
        ax.add_patch(Rectangle((b0, y - 0.08), b1 - b0, 0.16, facecolor='white', edgecolor=col, lw=1.0, zorder=2))
        ax.add_patch(Polygon([(q0 * kb, yB - 0.18), (q1 * kb, yB - 0.18), (b1, y + 0.1), (b0, y + 0.1)], closed=True,
                             facecolor=col, alpha=0.13, lw=0, zorder=1))
        chrom = end[:-1]
        c0, c1 = sorted([L.chm13_coord(fai, end, x0), L.chm13_coord(fai, end, x1)])
        note = f" ({x0 / 1e3:.0f}-{x1 / 1e3:.0f} kb from the CHM13 {end[3:]} telomere; drawn under the read)" if sh else ''
        ax.text(max(b0, 0), y - 0.16, f"{_end_label(end)}  {_fmt_coord(chrom, c0)}-{int(c1):,}{note}", ha='left', va='top', fontsize=6, color=col)
    # labels
    if an and full * kb > xmax:
        ax.annotate('', xy=(xmax, yR), xytext=(xmax - 0.04 * xmax, yR), arrowprops=dict(arrowstyle='->', color=CA, lw=1))
    if not compact:
        ax.text(t0 * kb, yR + 0.17, f'telomere {tl / 1e3:.1f} kb', fontsize=6, color=CT, ha='left', va='bottom')
        ax.text(sw * kb, yR + 0.42, f"switch, 95% interval {lo / 1e3:.1f}-{hi / 1e3:.1f} kb", fontsize=6, color='#555555',
                ha='center', va='bottom')
        if an:
            ax.text(xmax, yR + 0.17, f"anchor {r['anchor_chm13']}, MAPQ {int(r['anchor_mapq'])}, {int(r['anchor_unique_bp']) / 1e3:.0f} kb "
                    f"CHM13-unique; read continues to {full / 1e3:.0f} kb", fontsize=6, color=CA, ha='right', va='bottom')
        ax.text(-0.01, yS, 'SNV sites', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6)
        ax.text(-0.01, yB, 'private bins', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6)
        ax.text(-0.01, yR, 'read', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6)
    ax.set_xlim(-0.02 * xmax, xmax * 1.01); ax.set_ylim(-0.25, yR + 0.75)
    ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('read position from the telomeric end (kb)')


def _purity_note(r):
    return (f"donor segment {int(r['seg_D_units'])} units, {100 * r['seg_D_purity']:.0f}% donor alleles; proximal segment "
            f"{int(r['seg_P_units'])} units, {100 * r['seg_P_purity']:.0f}% recipient alleles; log10 LR {min(r['log10_lr_vs_allP'], r['log10_lr_vs_allD']):.0f}")


def showcase(a):
    import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    style(); C, Rr, U, fai = _load(); Ci = C.set_index('claim_id')
    S = pd.read_csv(os.path.join(REPO, 'evidence', 'showcase', 'showcase_claims.tsv'), sep='\t')
    os.makedirs(os.path.join(FIG, 'showcase'), exist_ok=True)
    ex = []
    for s in S.itertuples():
        c = Ci.loc[s.claim_id]; c = dict(c, claim_id=s.claim_id)
        sr = pd.read_csv(os.path.join(REPO, 'evidence', 'showcase', s.slug, f'{s.slug}.reads.tsv'), sep='\t')
        r = sr.iloc[0].to_dict(); k = f"{s.claim_id}|{r['dataset']}|{r['read_id']}"
        ex.append((c, r, U[k], s.slug))
    leg = [Line2D([], [], color=CD, lw=2, label='donor allele (CHM13 D)'), Line2D([], [], color=CP, lw=2, label='recipient / proximal allele (CHM13 P)'),
           Line2D([], [], color=CO, lw=2, label='neither'), Line2D([], [], color=CT, lw=4, label='telomere'),
           Line2D([], [], color=CA, lw=4, label='unique anchor'), Line2D([], [], color=CI_BAND, lw=6, label='switch 95% interval')]
    for c, r, u, slug in ex:
        fig, ax = plt.subplots(figsize=(7.2, 2.9))
        read_through(ax, c, r, u, fai)
        ax.set_title(f"{_label(c)}: {_end_label(c['donor_D'])} alleles distal to the switch, {_end_label(c['proximal_P'])} proximal", loc='left')
        fig.text(0.01, 0.125, f"read {r['read_id']} ({r['dataset']})", fontsize=6, color='#444444')
        fig.text(0.01, 0.085, _purity_note(r), fontsize=6, color='#444444')
        fig.legend(handles=leg, loc='lower center', ncol=6, fontsize=6, bbox_to_anchor=(0.5, 0.0))
        fig.subplots_adjust(left=0.1, right=0.99, top=0.9, bottom=0.33)
        for ext in ('png', 'svg'): fig.savefig(os.path.join(FIG, 'showcase', f'{slug}.readthrough.{ext}'), dpi=200)
        plt.close(fig)
    # combined panel figure (report / overview page)
    n = len(ex); fig, axs = plt.subplots(n, 1, figsize=(7.2, 1.55 * n + 0.5))
    for i, (ax, (c, r, u, slug)) in enumerate(zip(axs, ex)):
        read_through(ax, c, r, u, fai, compact=True)
        ax.set_title(f"{'abcdefghij'[i]}  {_label(c)}: {_end_label(c['donor_D'])} to {_end_label(c['proximal_P'])}, switch "
                     f"{r['switch_ci_read_lo'] / 1e3:.1f}-{r['switch_ci_read_hi'] / 1e3:.1f} kb, anchor {int(r['anchor_unique_bp']) / 1e3:.0f} kb unique "
                     f"({r['dataset']})", loc='left', fontsize=7)
        if i < n - 1: ax.set_xlabel('')
    fig.legend(handles=leg, loc='upper center', ncol=6, fontsize=6, bbox_to_anchor=(0.5, 0.999))
    fig.subplots_adjust(left=0.02, right=0.99, top=1 - 0.45 / (1.55 * n + 0.5), bottom=0.35 / (1.55 * n + 0.5), hspace=0.75)
    for ext in ('png', 'pdf', 'svg'): fig.savefig(os.path.join(FIG, f'fig_showcase_readthrough.{ext}'), dpi=200)
    plt.close(fig)
    print(f'showcase figures: {n}', file=sys.stderr)


# --------------------------------------------------------------------------------------------------------------- (b)
TIER_CLASS = [('definitive, >= 2 datasets', '#1b6e3c'), ('definitive', '#4daf6a'), ('confirmed, not definitive', '#a6d8b0'),
              ('supported (v1)', '#f0e3a0'), ('assessed, not confirmed', '#c9c9c9'), ('not assessable', '#ececec')]


def _tier_class(c):
    d, v = c['definitive_verdict'], c['v1_verdict']
    if d == 'definitive_multi_dataset': return TIER_CLASS[0][0]
    if d == 'definitive': return TIER_CLASS[1][0]
    if v == 'confirmed': return TIER_CLASS[2][0]
    if v == 'supported': return TIER_CLASS[3][0]
    if str(c['v1_class']).startswith('not_assessable'): return TIER_CLASS[5][0]
    return TIER_CLASS[4][0]


def summary(a):
    import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from scipy.stats import beta
    style(); C, Rr, U, fai = _load()
    C['tclass'] = [_tier_class(c) for c in C.to_dict('records')]
    fig = plt.figure(figsize=(7.2, 6.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], width_ratios=[1.15, 1], hspace=0.75, wspace=0.55)
    # a: claims by cohort and tier
    ax = fig.add_subplot(gs[0, 0])
    coh = [('HG002_Q100', 'HG002 Q100'), ('HG002_trio', 'HG003/HG004\n(transmitting)'), ('1KG_ONT', '1KG-ONT\n(42 people)')]
    for i, (k, lab) in enumerate(coh):
        g = C[C.cohort == k]
        if k == 'HG002_trio': g = g[g.v1_class == 'transmitting_parent']
        x0 = 0
        for t, col in TIER_CLASS:
            n = int((g.tclass == t).sum())
            if n:
                ax.barh(i, n, left=x0, color=col, edgecolor='white', lw=0.5, height=0.62)
                if n >= (3 if k != '1KG_ONT' else 14): ax.text(x0 + n / 2, i, str(n), ha='center', va='center', fontsize=6,
                                                               color='white' if col in ('#1b6e3c', '#4daf6a') else 'black')
                x0 += n
        ax.text(x0 * 1.02 + 2, i, f'{len(g)} claims', va='center', fontsize=6)
    ax.set_yticks(range(len(coh))); ax.set_yticklabels([l for _, l in coh]); ax.invert_yaxis(); ax.set_xscale('symlog', linthresh=30)
    ax.set_xticks([0, 10, 20, 30, 100, 300, 800]); ax.set_xticklabels(['0', '10', '20', '30', '100', '300', '800'])
    ax.set_xlabel('claims (linear to 30, log above)')
    ax.set_title('Definitive claims: a subset of v1-confirmed')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c_, label=t) for t, c_ in TIER_CLASS], loc='upper left', bbox_to_anchor=(-0.25, -0.3), ncol=3, fontsize=6)
    # b: controls
    ax = fig.add_subplot(gs[0, 1])
    K = pd.read_csv(os.path.join(REPO, 'evidence', 'controls.tsv'), sep='\t')
    ctl = [('HG002 other haplotype,\neligible reads', K[K.cohort == 'HG002_Q100'].n_eligible.sum(), K[K.cohort == 'HG002_Q100'].n_conclusive.sum()),
           ('1KG-ONT native,\nread x model tests', K[K.cohort == '1KG_ONT'].n_eligible.sum(), K[K.cohort == '1KG_ONT'].n_conclusive.sum())]
    hgel = int(Rr[(Rr.claim_id.str.startswith('control|HG002')) & Rr.v1_verdict.isin(['contradicts'])].shape[0])
    for i, (lab, n, k) in enumerate(ctl):
        ub = 1.0 if n == 0 else float(beta.ppf(0.95, k + 1, n - k))
        ax.plot([0, ub * 100], [i, i], color='#555555', lw=1.2); ax.plot([k / n * 100 if n else 0], [i], 'o', color='black', ms=4)
        ax.text(ub * 100 + 0.05, i, f'{k}/{n:,} conclusive\n95% upper bound {ub * 100:.2f}%', va='center', fontsize=6)
    ax.set_yticks(range(len(ctl))); ax.set_yticklabels([c_[0] for c_ in ctl]); ax.invert_yaxis(); ax.set_xlim(-0.05, 4.2)
    ax.set_ylim(len(ctl) - 0.4, -0.6)
    ax.set_xlabel('false conclusive rate (%)'); ax.set_title('No conclusive read in negative controls')
    # c: HG002 dataset independence of definitive reads
    ax = fig.add_subplot(gs[1, 0])
    H = C[(C.cohort == 'HG002_Q100') & C.definitive_verdict.str.startswith('definitive')].copy()
    hr = Rr[(Rr.cohort == 'HG002_Q100') & (Rr.tier == 'definitive') & (Rr.counted == True)]  # noqa: E712
    ds_order = ['ont_ul_dorado', 'ont_lc24_ulk', 'ont_duplex', 'ont_q28', 'ont_lc24_apk', 'ont_q28_lc24', 'ont_r941_recall',
                'hifi_revio_3cells', 'hifi_revio_pbmay24', 'hifi_revio_30hr', 'hifi_vega_24hr', 'hifi_sprq_30hr', 'hifi_sprq_24hr', 'hifi_pb_sprq', 'hifi_dc']
    ds_order = [d for d in ds_order if d in set(hr.dataset)]
    H['lab'] = [f"{c['arm']} {c['haplotype'][:3]}" for c in H.to_dict('records')]
    H = H.sort_values('n_definitive', ascending=False)
    M = np.array([[int(((hr.claim_id == cid) & (hr.dataset == d)).sum()) for d in ds_order] for cid in H.claim_id])
    ax.imshow(np.where(M > 0, np.log10(M + 1), np.nan), cmap='Greens', vmin=0, vmax=2.2, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j]: ax.text(j, i, str(M[i, j]), ha='center', va='center', fontsize=5.5, color='white' if M[i, j] >= 20 else 'black')
    ax.set_xticks(range(len(ds_order))); ax.set_xticklabels([d.replace('_', ' ') for d in ds_order], rotation=55, ha='right', fontsize=6)
    ax.set_yticks(range(len(H))); ax.set_yticklabels(H.lab, fontsize=6)
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False); ax.tick_params(length=0)
    ax.set_title(f"HG002 definitive reads per dataset ({int((H.n_definitive_datasets >= 2).sum())} of {len(H)} claims: >= 2 datasets)", fontsize=7)
    # d: trio
    ax = fig.add_subplot(gs[1, 1])
    TR = C[C.cohort == 'HG002_trio'].copy(); TR['hg'] = TR.claim_id.str.split('|').str[1]
    hc = C[(C.cohort == 'HG002_Q100') & (C.v1_verdict.isin(['confirmed', 'supported']))].sort_values('n_definitive', ascending=False)
    vals = {'definitive_multi_dataset': 3, 'definitive': 3, 'confirmed_not_definitive': 2}
    def sc(r):
        return vals.get(r['definitive_verdict'], 1 if r['v1_verdict'] == 'supported' else 0)
    M2 = []
    for cid in hc.claim_id:
        t = TR[(TR.hg == cid) & (TR.v1_class == 'transmitting_parent')]; o = TR[(TR.hg == cid) & (TR.v1_class == 'other_parent')]
        M2.append([sc(t.iloc[0]) if len(t) else -1, sc(o.iloc[0]) if len(o) else -1])
    M2 = np.array(M2)
    from matplotlib.colors import ListedColormap
    cm = ListedColormap(['#ececec', '#f0e3a0', '#a6d8b0', '#1b6e3c'])
    ax.imshow(M2, cmap=cm, vmin=-0.5, vmax=3.5, aspect='auto')
    txt = {0: 'not confirmed', 1: 'supported', 2: 'confirmed', 3: 'definitive'}
    for i in range(M2.shape[0]):
        for j in range(2): ax.text(j, i, txt[M2[i, j]], ha='center', va='center', fontsize=5.5, color='white' if M2[i, j] == 3 else 'black')
    ax.set_xticks([0, 1]); ax.set_xticklabels(['transmitting parent', 'other parent'], fontsize=6)
    ax.set_yticks(range(len(hc))); ax.set_yticklabels([f"{c['arm']} {c['haplotype'][:3]}" for c in hc.to_dict('records')], fontsize=6)
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False); ax.tick_params(length=0)
    nt = int((M2[:, 0] >= 2).sum()); ntd = int((M2[:, 0] == 3).sum())
    ax.set_title(f'Trio: transmitting parent confirms {nt} of {len(hc)}', fontsize=7)
    for i, l in enumerate('abcd'):
        fig.axes[i].text(-0.02, 1.08, l, transform=fig.axes[i].transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')
    os.makedirs(FIG, exist_ok=True)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.14)
    for ext in ('png', 'pdf', 'svg'): fig.savefig(os.path.join(FIG, f'fig_evidence_summary.{ext}'), dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------------------- (c)
MAX_ROWS = 12


def claim_panel(c, rows, U, fai, fn):
    import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    n = min(len(rows), MAX_ROWS); kb = 1e-3
    fig, ax = plt.subplots(figsize=(7.2, 0.62 + 0.24 * n))
    xmax = 0
    for i, r in enumerate(rows[:n]):
        u = U[f"{c['claim_id']}|{r['dataset']}|{r['read_id']}"]; s0 = u['s0']; A = _units_split(u['units'])
        y = n - 1 - i; an = u['anchor_read']
        end = (an[1] if an else u['read_len']) - s0; xmax = max(xmax, end)
        ax.add_patch(Rectangle(((r['switch_ci_read_lo'] - s0) * kb, y - 0.42), (r['switch_ci_read_hi'] - r['switch_ci_read_lo']) * kb, 0.84,
                               color=CI_BAND, alpha=0.6, lw=0, zorder=0))
        ax.plot([0, (u['read_len'] - s0) * kb], [y - 0.28, y - 0.28], color='#cfcfcf', lw=0.8, zorder=1)
        if an: ax.plot([(an[0] - s0) * kb, (an[1] - s0) * kb], [y - 0.28, y - 0.28], color=CA, lw=2.2, solid_capstyle='butt', zorder=2)
        for lab, col in ((2, CO), (0, CP), (1, CD)):
            m = A[:, 1] == lab
            if m.any(): ax.vlines((A[m, 0] - s0) * kb, y - 0.15, y + 0.3, color=col, lw=0.5, rasterized=True, zorder=3)
        ax.text(1.0, y, f"{r['dataset']}", transform=ax.get_yaxis_transform(), fontsize=5.5, va='center', ha='left', color='#444444')
    sw = str(c['asm_projection_switch_ci_sub'])
    if '-' in sw and sw.split('-')[0] not in ('', 'None', 'nan'):
        a0, a1 = (float(x) for x in sw.split('-')[:2])
        ax.axvspan(a0 * kb, a1 * kb, color='black', alpha=0.08, lw=0)
        ax.text(a1 * kb, n - 0.45, ' assembly switch', fontsize=5.5, va='bottom', ha='left', color='#333333')
    ax.set_xlim(0, xmax * kb * 1.01); ax.set_ylim(-0.6, n - 0.4 + 0.35); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('position from the telomere boundary (kb)')
    more = f'; {len(rows) - n} more in reads.tsv' if len(rows) > n else ''
    ax.set_title(f"{_label(c)}: {len(rows)} definitive reads ({_end_label(c['donor_D'])} alleles orange, {_end_label(c['proximal_P'])} blue; "
                 f"switch interval yellow, unique anchor green{more})", fontsize=7)
    fig.subplots_adjust(left=0.02, right=0.86, top=1 - 0.32 / (0.62 + 0.24 * n), bottom=0.36 / (0.62 + 0.24 * n) + 0.05)
    fig.savefig(fn + '.png', dpi=150); fig.savefig(fn + '.svg'); plt.close(fig)


def claims(a):
    style(); C, Rr, U, fai = _load()
    od = os.path.join(FIG, 'claims'); os.makedirs(od, exist_ok=True)
    D = C[C.definitive_verdict.str.startswith('definitive')]
    for c in D.to_dict('records'):
        rr = Rr[(Rr.claim_id == c['claim_id']) & (Rr.tier == 'definitive') & (Rr.counted == True)].copy()  # noqa: E712
        rr = rr.assign(_u=rr.seg_D_units.fillna(0) + rr.seg_P_units.fillna(0)).sort_values('_u', ascending=False)
        claim_panel(c, rr.to_dict('records'), U, fai, os.path.join(od, c['slug']))
    print(f'{len(D)} claim panels', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['extract', 'showcase', 'summary', 'claims', 'all'])
    a = ap.parse_args()
    for s in (['extract', 'showcase', 'summary', 'claims'] if a.stage == 'all' else [a.stage]): globals()[s](a)


if __name__ == '__main__':
    main()
