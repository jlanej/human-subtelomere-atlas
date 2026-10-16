#!/usr/bin/env python3
"""Calibration summary tables and figures (evidence/calibration/).

Inputs: the training and held-out runs of run_calibration.py (reads.tsv, claims.tsv, units.json.gz) and site_error_rates.py.
Outputs (evidence/calibration/): reads_<half>.tsv.gz, claims_calibration.tsv, read_metrics.tsv, sensitivity_by_length.tsv,
fig_allele_tracks.{png,pdf}, fig_sensitivity_specificity.{png,pdf}.
Definitions (CRITERIA.md section 5):
  eligible carrier read = terminal telomere and read length >= s0_read + switch_sub + anchor_min_bp (a read that could cross the
  switch and carry a full anchor); sensitivity = conclusive and switch-concordant / eligible.
  false conclusive rate = conclusive reads / reads among strictly native controls (all reads; and terminal reads that reach the
  claimed switch + anchor length); 95% upper bound = exact (Clopper-Pearson) one-sided.
"""
import argparse, os, sys, json, gzip
import numpy as np, pandas as pd
from scipy.stats import beta
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DON, PRO, OTH = '#D55E00', '#0072B2', '#9a9a9a'
PLAT = {'ONT_R10_UL': ('ONT R10 ultra-long', '#000000'), 'ONT_R10_duplex': ('ONT R10 duplex', '#56B4E9'),
        'HiFi_Revio': ('PacBio HiFi Revio', '#009E73'), 'ONT_R9.4.1': ('ONT R9.4.1 (1KG)', '#CC79A7')}
NEG = ('hg002_otherhap_native', 'ont1kg_native_strict')


def cp_upper(k, n, a=0.05):
    return 1.0 if n == 0 else (1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k)))


def cp_ci(k, n, a=0.05):
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1)); hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def style():
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'DejaVu Sans',
                         'axes.titlelocation': 'left', 'savefig.dpi': 300, 'pdf.fonttype': 42})


def load(run):
    R = pd.read_csv(os.path.join(run, 'reads.tsv'), sep='\t', low_memory=False)
    C = pd.read_csv(os.path.join(run, 'claims.tsv'), sep='\t')
    U = json.load(gzip.open(os.path.join(run, 'units.json.gz'), 'rt'))
    return R, C, U


def metrics(R, th):
    R = R.copy()
    R['eligible'] = R.terminal & R.assembly_switch_sub.notna() & (R.read_len >= R.s0_read + R.assembly_switch_sub.fillna(0) + th['anchor_min_bp'])
    out = []
    for half in ['train', 'heldout', 'all']:
        X = R if half == 'all' else R[R.half == half]
        for plat in sorted(X.platform.unique()) + ['all']:
            Y = X if plat == 'all' else X[X.platform == plat]
            car = Y[Y.role == 'carrier']; el = car[car.eligible]
            k = int(((el.verdict == 'conclusive') & (el.switch_concordant == True)).sum())  # noqa: E712
            neg = Y[Y.role.isin(NEG)]; negel = neg[neg.eligible]
            kn = int((neg.verdict == 'conclusive').sum()); kne = int((negel.verdict == 'conclusive').sum())
            lo, hi = cp_ci(k, len(el))
            out.append(dict(half=half, platform=plat, carrier_reads=len(car), carrier_eligible=len(el), conclusive_concordant_eligible=k,
                            sensitivity=round(k / len(el), 4) if len(el) else None, sens_ci_lo=round(lo, 4), sens_ci_hi=round(hi, 4),
                            carrier_conclusive_any=int((car.verdict == 'conclusive').sum()),
                            carrier_contradicts=int((car.verdict == 'contradicts').sum()),
                            neg_reads=len(neg), neg_conclusive=kn, neg_false_conclusive_rate=round(kn / len(neg), 5) if len(neg) else None,
                            neg_fcr_upper95=round(cp_upper(kn, len(neg)), 5) if len(neg) else None,
                            neg_eligible=len(negel), neg_eligible_conclusive=kne,
                            neg_eligible_fcr_upper95=round(cp_upper(kne, len(negel)), 5) if len(negel) else None,
                            neg_consistent=int((neg.verdict == 'consistent').sum()), neg_contradicts=int((neg.verdict == 'contradicts').sum())))
    return R, pd.DataFrame(out)


def by_length(R):
    bins = [0, 25e3, 50e3, 100e3, 200e3, 1e9]; lab = ['<25 kb', '25-50 kb', '50-100 kb', '100-200 kb', '>200 kb']
    car = R[(R.role == 'carrier') & R.terminal & R.assembly_switch_sub.notna() & (R.read_len >= R.s0_read + R.assembly_switch_sub)].copy()
    car['len_bin'] = pd.cut(car.read_len, bins, labels=lab)
    g = car.groupby(['platform', 'len_bin'], observed=False)
    t = g.apply(lambda x: pd.Series(dict(n=len(x), conclusive=int(((x.verdict == 'conclusive') & (x.switch_concordant == True)).sum()),  # noqa: E712
                                         switch_pass=int(x.switch_pass.sum()), anchor_pass=int(x.anchor_pass.sum())))).reset_index()
    t['sensitivity'] = (t.conclusive / t.n).round(4)
    return t


def fig_tracks(sel, out):
    """sel: list of (title, [(label, row, units)]) for 3 exemplars"""
    style()
    nrow = sum(len(r) for _, r in sel)
    fig, axes = plt.subplots(len(sel), 1, figsize=(7.2, 0.42 * nrow + 0.9 * len(sel)),
                             gridspec_kw=dict(height_ratios=[len(r) for _, r in sel], hspace=0.75))
    for ax, (title, rows) in zip(np.atleast_1d(axes), sel):
        xmax = 0
        for i, (lab, r, u) in enumerate(rows):
            y = len(rows) - 1 - i; s0 = r['s0_read']
            for kind, mk, dy in (('s', '|', 0.0), ('b', 's', 0.0)):
                for al, col in (('D', DON), ('P', PRO), ('o', OTH)):
                    xs = [(p - s0) / 1e3 for p, a, k, *rest in u if k == kind and a == al]
                    if xs:
                        ax.scatter(xs, [y + dy] * len(xs), marker=mk, s=(18 if kind == 's' else 5), c=col, lw=0.6 if kind == 's' else 0,
                                   zorder=3)
            if r.get('switch_pass') and pd.notna(r.get('switch_ci_sub_lo')):
                ax.add_patch(plt.Rectangle((r['switch_ci_sub_lo'] / 1e3, y - 0.38), max(0.15, (r['switch_ci_sub_hi'] - r['switch_ci_sub_lo']) / 1e3),
                                           0.76, color='#f0e442', alpha=0.55, lw=0, zorder=1))
            xmax = max(xmax, max([(p - s0) / 1e3 for p, *_ in u] or [0]))
            ax.text(-0.01, y, lab, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6.5)
        sw = rows[0][1].get('assembly_switch_sub')
        if pd.notna(sw): ax.axvline(sw / 1e3, color='k', ls='--', lw=0.7, zorder=2)
        aj = rows[0][1].get('_asm')
        if aj: ax.axvspan(aj[0] / 1e3, max(aj[1], aj[0] + 150) / 1e3, ymin=0.97, ymax=1.0, color='k', lw=0, zorder=2)
        ax.set_yticks([]); ax.spines['left'].set_visible(False)
        ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_xlim(-1, min(xmax, rows[0][1]['assembly_switch_sub'] / 1e3 * 2.2 + 10) if pd.notna(sw) else xmax)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel('position on read from the telomere boundary (kb)')
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker='|', ls='', color=DON, ms=7, label='donor (D) allele, SNV site'),
         Line2D([], [], marker='|', ls='', color=PRO, ms=7, label='proximal (P) allele, SNV site'),
         Line2D([], [], marker='s', ls='', color=DON, ms=3, label='D-private 200 bp bin'),
         Line2D([], [], marker='s', ls='', color=PRO, ms=3, label='P-private 200 bp bin'),
         Line2D([], [], marker='|', ls='', color=OTH, ms=7, label='neither allele'),
         plt.Rectangle((0, 0), 1, 1, color='#f0e442', alpha=0.55, label='read switch, 95% interval'),
         Line2D([], [], color='k', ls='--', lw=0.7, label='catalogue switch'),
         plt.Rectangle((0, 0), 1, 1, color='k', label='assembly window, same unit model (bar at top)')]
    fig.legend(handles=h, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=6.5)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.95, bottom=0.12)
    for ext in ('png', 'pdf'): fig.savefig(f'{out}.{ext}', bbox_inches='tight')
    plt.close(fig)


def fig_sens_spec(L, Mt, C, out):
    style()
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1.2, 1, 1], wspace=0.85))
    labs = ['<25 kb', '25-50 kb', '50-100 kb', '100-200 kb', '>200 kb']
    for p, (name, col) in PLAT.items():
        t = L[L.platform == p]
        if not len(t): continue
        t = t.set_index('len_bin').reindex(labs)
        x = np.arange(len(labs)); ok = t.n.fillna(0) >= 5
        a1.plot(x[ok], t.sensitivity[ok], '-o', color=col, ms=3, lw=1, label=name)
    a1.set_xticks(range(len(labs))); a1.set_xticklabels(labs, rotation=30, ha='right')
    a1.set_ylim(-0.03, 1.03); a1.set_ylabel('fraction conclusive and concordant')
    a1.set_title('Read sensitivity rises with length'); a1.legend(frameon=False, loc='upper left', fontsize=6)
    a1.set_xlabel('read length\n(terminal carrier reads reaching the switch)')
    m = Mt[(Mt.half == 'all') & (Mt.platform != 'all') & (Mt.neg_reads > 0)]
    y = np.arange(len(m))
    for i, r in enumerate(m.itertuples()):
        col = PLAT.get(r.platform, (r.platform, 'k'))[1]
        a2.plot([0, r.neg_fcr_upper95 * 100], [i, i], color=col, lw=1.2)
        a2.plot([r.neg_false_conclusive_rate * 100], [i], 'o', color=col, ms=4)
        a2.text(r.neg_fcr_upper95 * 100 + 0.3, i, f'{r.neg_conclusive}/{r.neg_reads}', va='center', fontsize=6.5)
    a2.set_yticks(y); a2.set_yticklabels([PLAT.get(p, (p,))[0] for p in m.platform]); a2.set_ylim(-0.6, len(m) - 0.4)
    a2.set_xlim(0, max(10, (m.neg_fcr_upper95.max() * 100) * 1.6)); a2.set_xlabel('false conclusive reads (%)\nline = exact 95% upper bound')
    a2.set_title('Strictly native controls: 0 conclusive')
    cc = C.copy(); cc['label'] = cc.claim_id.str.replace('HG002#', '').str.replace('_PATERNAL', ' pat').str.replace('_MATERNAL', ' mat') \
        .str.replace(r'^\d#chr', '', regex=True).str.replace('#L', ' p').str.replace('#R', ' q')
    cc = cc.sort_values(['half', 'n_conclusive_concordant'])
    col = {'confirmed': '#000000', 'supported': '#777777', 'not_confirmed': '#ffffff', 'contradicted': '#D55E00'}
    for i, r in enumerate(cc.itertuples()):
        a3.barh(i, max(r.n_conclusive_concordant, 0.0), color='#444444' if r.half == 'train' else '#56B4E9', height=0.7)
        a3.text(max(r.n_conclusive_concordant, 0) + 0.8, i, r.verdict.replace('_', ' '), va='center', fontsize=5.5)
    a3.set_yticks(range(len(cc))); a3.set_yticklabels(cc.label, fontsize=5.5); a3.set_ylim(-0.6, len(cc) - 0.4)
    a3.set_xlabel('conclusive, switch-\nconcordant reads per claim'); a3.set_xlim(0, cc.n_conclusive_concordant.max() * 1.45)
    a3.set_title('HG002 claims (grey train, blue held-out)')
    for ext in ('png', 'pdf'): fig.savefig(f'{out}.{ext}', bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--train', required=True); ap.add_argument('--heldout', required=True)
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--out', default=os.path.join(REPO, 'evidence', 'calibration'))
    ap.add_argument('--exemplars', default='HG002#2#chr9_MATERNAL#L,HG002#2#chr17_MATERNAL#L,HG002#1#chr8_PATERNAL#L')
    a = ap.parse_args()
    th = json.load(open(a.criteria))['thresholds']; os.makedirs(a.out, exist_ok=True)
    Rt, Ct, Ut = load(a.train); Rh, Ch, Uh = load(a.heldout)
    keep = [c for c in Rt.columns if c not in ('genome_alignments', 'proj_hits')]
    for R, h in ((Rt, 'train'), (Rh, 'heldout')):
        R[keep].to_csv(os.path.join(a.out, f'reads_{h}.tsv.gz'), sep='\t', index=False, compression='gzip')
    R = pd.concat([Rt, Rh], ignore_index=True); C = pd.concat([Ct, Ch], ignore_index=True)
    C.to_csv(os.path.join(a.out, 'claims_calibration.tsv'), sep='\t', index=False)
    R, Mt = metrics(R, th); Mt.to_csv(os.path.join(a.out, 'read_metrics.tsv'), sep='\t', index=False)
    L = by_length(R); L.to_csv(os.path.join(a.out, 'sensitivity_by_length.tsv'), sep='\t', index=False)
    print(Mt.to_string())
    # exemplar tracks: up to 4 conclusive carrier reads (longest, distinct datasets first) + 1 strictly native control read
    sel = []
    U = {**{('train', k): v for k, v in Ut.items()}, **{('heldout', k): v for k, v in Uh.items()}}
    for cid in a.exemplars.split(','):
        X = R[R.claim_id == cid]; half = X.half.iloc[0]
        car = X[(X.role == 'carrier') & (X.verdict == 'conclusive') & (X.switch_concordant == True)]  # noqa: E712
        pick = pd.concat([car[car.platform == p].sort_values('read_len', ascending=False).head(2 if p == 'ONT_R10_UL' else 1)
                          for p in ['ONT_R10_UL', 'ONT_R10_duplex', 'HiFi_Revio']]).head(4)
        ctl = X[X.role.isin(NEG) & X.terminal].assign(_c=lambda d: d.verdict == 'contradicts').sort_values(['_c', 'n_units'], ascending=False).head(1)
        rows = []
        for r in list(pick.itertuples()) + list(ctl.itertuples()):
            d = r._asdict(); u = U[(half, f"{int(d['job'])}|{d['read']}")]
            tag = 'native control ' + (d['sample'] if d['role'].startswith('ont1kg') else 'HG002 other haplotype') if d['role'] in NEG else 'HG002 carrier'
            rows.append((f"{tag}\n{PLAT[d['platform']][0]}, {d['read_len'] / 1e3:.0f} kb: {d['verdict']}", d, u))
        c = X.iloc[0]
        hap = 'paternal' if 'PATERNAL' in cid else 'maternal'
        arm = cid.split('#')[2].split('_')[0][3:] + ('p' if cid.endswith('#L') else 'q')
        asm = C[C.claim_id == cid].iloc[0]
        rows[0][1]['_asm'] = (asm.asm_switch_ci_sub_lo, asm.asm_switch_ci_sub_hi) if asm.asm_switch_pass else None
        sel.append((f"HG002 {hap} {arm} ({half}): telomere-adjacent block {c.D[3:]}-like, then {c.P[3:]}-like", rows))
    fig_tracks(sel, os.path.join(a.out, 'fig_allele_tracks'))
    fig_sens_spec(L, Mt, C, os.path.join(a.out, 'fig_sensitivity_specificity'))


if __name__ == '__main__':
    main()
