"""Step 7: headline summary (data/validation/tan2026_summary.json) and figure (figures/validation/tan2026_concordance.png/.pdf).
Inputs: work/tan2026/summary_crosswalk.json, summary_cohorts.json (steps 5-6), data/validation/tan2026_crosswalk.tsv,
        data/validation/tan2026_cohort_forms.tsv.
Figure: (a) agreement of catalogue calls on the same haplotype sequence, by assembler pairing (primary pairs: verified
        same-haplotype sequence, telomere present in both; open circles: hifiasm vs Verkko pairs whose telomere boundary call is
        identical; grey crosses: hifiasm vs Verkko against catalogue revision 1, when its summary is available locally); (b) base-level switch positions (both sides base-level, independent-assembler pairs) after projection through
        the pairwise alignment and as measured from each record's own telomere boundary; (c) structural-allele frequencies in the pooled independent cohorts versus HPRC Release 2;
        (d) frequencies of four common alleles per cohort with 95% Wilson intervals.
Usage : python tan2026_07_report.py
"""
import os, sys, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
from tan2026_06_cohorts import KEY_ALLELES

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 8,
                     'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
CMP = {'hifiasm_vs_verkko': ('hifiasm vs Verkko (independent assembler)', '#1f5fa8'),
       'verkko_vs_hifiasm': ('Verkko vs hifiasm (independent assembler)', '#7aa6d6'),
       'hifiasm_vs_hifiasm': ('hifiasm vs hifiasm (same assembly)', '#8c8c8c')}
TIER = {'A': '#1f5fa8', 'B': '#d9822b'}

REV1 = f'{T.TW}/rev1/tan2026_summary.json'   # summary of the same cross-check against catalogue revision 1 (kept locally)
R1 = None

def before_after(SX, SC):
    """headline numbers against catalogue revision 1 (5 kb-bin mosaic, single best CHM13 end) and revision 2 (track-based, paralog
    families). Revision-1 S-state agreement compared state + best end; revision 2 compares state + family."""
    if not os.path.exists(REV1): return dict(note='revision-1 summary not available locally (work/tan2026/rev1/)')
    o = json.load(open(REV1)); a, b = o['same_haplotype'], SX; out = {}
    pairs = [('primary pairs', 'n_pairs', 'n_pairs'), ('S5 state', 'S5_state_and_end_agree', 'S5_state_and_family_agree'),
             ('S10 state', 'S10_state_and_end_agree', 'S10_state_and_family_agree'), ('S20 state', 'S20_state_and_end_agree', 'S20_state_and_family_agree'),
             ('form', 'form_agree', 'form_agree'), ('core structure', 'core_structure_agree', 'core_structure_agree'),
             ('structural allele', 'sa_agree', 'sa_agree'), ('allele compatible', 'sa_compatible', 'sa_compatible'),
             ('Tan unassigned', 'frac_tan_unassigned', 'frac_tan_unassigned'), ('base-level switches compared', 'n_switches_base_both', 'n_switches_base_both'),
             ('switches identical after projection', 'switch_base_frac_proj_exact', 'switch_base_frac_proj_exact'),
             ('switches identical from own boundary', 'switch_base_frac_raw_exact', 'switch_base_frac_raw_exact')]
    out['same_haplotype_primary'] = {nm: dict(revision1=a['primary'].get(k1), revision2=b['primary'].get(k2)) for nm, k1, k2 in pairs}
    out['allele_agreement_by_boundary_delta'] = {k: dict(revision1=a['by_boundary_delta'].get(k, {}).get('sa_agree'), revision2=v.get('sa_agree'),
                                                     n_revision2=v.get('n_pairs')) for k, v in b['by_boundary_delta'].items()}
    out['allele_agreement_by_tier'] = {k: dict(revision1=a['by_tier'].get(k, {}).get('sa_agree'), revision2=v.get('sa_agree'), n_revision2=v.get('n_pairs'))
                                       for k, v in b['by_tier'].items()}
    oc = o['independent_cohorts']
    out['independent_cohorts'] = {k: dict(revision1=oc.get(k), revision2=SC.get(k)) for k in
                                  ('n_primary', 'frac_assigned_by_cohort', 'replication_tierA_common', 'tierAB_observed_in_independent',
                                   'freq_correlation_non_native', 'n_sa_tests', 'n_sa_tests_q_lt_0.05', 'n_form_tests', 'n_form_tests_q_lt_0.05', 'novel')}
    out['note'] = ('Revision-1 S-state agreement = same state and same best CHM13 end; revision 2 = same state and same paralog family. '
                   'Allele IDs differ between revisions (id_history.tsv); agreements are within-revision.')
    return out

def catalogue_id():
    """catalogue release used for the comparison: version (PROVENANCE.json) and SHA-256 of the files read"""
    import catalogue_common as C
    cat = os.path.join(T.REPO, 'data', 'catalogue'); out = {}
    try: out['version'] = json.load(open(os.path.join(cat, 'PROVENANCE.json'))).get('release', '')
    except Exception: out['version'] = ''
    out['sha256'] = {f: C.sha256_file(os.path.join(cat, f)) for f in ('ends.tsv', 'alleles_structural.tsv', 'allele_frequencies_by_superpop.tsv', 'SHA256SUMS')
                     if os.path.exists(os.path.join(cat, f))}
    return out

def main():
    SX = json.load(open(f'{T.TW}/summary_crosswalk.json')); SC = json.load(open(f'{T.TW}/summary_cohorts.json'))
    global R1
    R1 = json.load(open(REV1)) if os.path.exists(REV1) else None
    X = pd.read_csv(f'{T.VAL}/tan2026_crosswalk.tsv', sep='\t', low_memory=False)
    R = pd.read_csv(f'{T.VAL}/tan2026_cohort_forms.tsv', sep='\t')
    out = dict(source=dict(study='Tan et al. 2026, bioRxiv', doi='10.64898/2026.09.03.748401', data_doi='10.5281/zenodo.21807156',
                           data_license='CC BY 4.0', archive='subtelomere_500kb_2kbtelo.zip'),
               catalogue=catalogue_id(), method='catalogue revision 2: track-based mosaic, paralog-family labels (subtel_atlas.tracks)',
               revision1_vs_revision2=before_after(SX, SC), same_haplotype=SX, independent_cohorts=SC)
    json.dump(out, open(f'{T.VAL}/tan2026_summary.json', 'w'), indent=1, sort_keys=False)

    fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.78, wspace=0.42)
    # (a) agreement by call type
    ax = fig.add_subplot(gs[0, 0])
    items = [('S5_state_and_family_agree', 'S5 state'), ('S10_state_and_family_agree', 'S10 state'), ('S20_state_and_family_agree', 'S20 state'),
             ('form_agree', 'form'), ('core_structure_agree_common_length', 'core structure'), ('sa_agree', 'structural allele'),
             ('sa_compatible', 'allele compatible')]
    y = np.arange(len(items))[::-1]; bc = SX['by_comparison']; off = {k: d for k, d in zip(CMP, (-0.12, 0.08, 0.28))}
    for k, (lab, col) in CMP.items():
        if k not in bc: continue
        v = [bc[k].get(i) for i, _ in items]
        ax.scatter([np.nan if a is None else a for a in v], y + off[k], s=18, color=col, zorder=3, label=f"{lab}, n={bc[k]['n_pairs']:,}")
    be = SX['by_comparison_and_boundary'].get('hifiasm_vs_verkko|exact')
    if be:
        v = [be.get(i) for i, _ in items]
        ax.scatter([np.nan if a is None else a for a in v], y - 0.3, s=18, facecolor='white', edgecolor=CMP['hifiasm_vs_verkko'][1], zorder=3,
                   label=f"hifiasm vs Verkko, identical boundary call, n={be['n_pairs']:,}")
    if R1 and 'hifiasm_vs_verkko' in R1['same_haplotype']['by_comparison']:
        r1 = R1['same_haplotype']['by_comparison']['hifiasm_vs_verkko']; k1 = {'S5_state_and_family_agree': 'S5_state_and_end_agree',
              'S10_state_and_family_agree': 'S10_state_and_end_agree', 'S20_state_and_family_agree': 'S20_state_and_end_agree'}
        v = [r1.get(k1.get(i, i)) for i, _ in items]
        ax.scatter([np.nan if a is None else a for a in v], y - 0.46, s=14, marker='x', color='#b0b0b0', zorder=2, lw=0.8,
                   label='hifiasm vs Verkko, catalogue revision 1 (5 kb bins, single CHM13 end)')
    base = SX['self_assignment_baseline']['frac_rule_returns_own_allele']
    ax.plot([base, base], [y[5] - 0.4, y[5] + 0.4], color='black', lw=1, label=f'allele-rule ceiling (Release 2 self-assignment, {base:.3f})')
    ax.set_yticks(y); ax.set_yticklabels([b for _, b in items]); ax.set_xlabel('fraction of same-haplotype pairs agreeing')
    lo = min([bc[k].get(i) or 1 for k in bc for i, _ in items]); ax.set_xlim(max(0, lo - 0.05), 1.005)
    ax.set_title('Catalogue calls reproduce on the same haplotype', loc='left')
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(-0.05, -0.2), fontsize=6, ncol=1)
    ax.text(-0.32, 1.06, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold')
    # (b) switch deltas (base-level on both sides, independent-assembler primary pairs; from tan2026_05_compare.py)
    ax = fig.add_subplot(gs[0, 1]); SW = pd.read_csv(f'{T.TW}/switch_deltas_primary.tsv', sep='\t')
    SW = SW[SW.base.astype(bool) & SW.comparison.isin(['hifiasm_vs_verkko', 'verkko_vs_hifiasm'])]
    for col, lab, c in (('d_proj', 'after alignment projection', CMP['hifiasm_vs_verkko'][1]), ('d_raw', 'from each telomere boundary', '#8c8c8c')):
        d = SW[col].dropna().abs().values
        xs = np.sort(d) + 1; ax.step(xs, np.arange(1, len(xs) + 1) / len(xs), where='post', color=c, label=lab)
    d = SW.d_proj.dropna().abs().values
    ax.set_xscale('log'); ax.set_xlim(0.8, 3e4); ax.set_xticks([1, 11, 101, 1001, 10001]); ax.set_xticklabels(['0', '10', '100', '1k', '10k'])
    ax.text(0.97, 0.30, f'{(d == 0).mean():.0%} at the identical base\nn={len(d):,} switches in\n{SW.end_id.nunique():,} independent-assembler pairs',
            transform=ax.transAxes, ha='right', fontsize=7)
    ax.legend(frameon=False, loc='lower right', fontsize=6)
    ax.set_xlabel('|switch position difference| (bp)'); ax.set_ylabel('cumulative fraction'); ax.set_ylim(0, 1.02)
    ax.set_title('Switches fall at the same base', loc='left')
    ax.text(-0.3, 1.06, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold')
    # (c) frequency scatter
    ax = fig.add_subplot(gs[1, 0])
    A = R[(R.kind == 'structural_allele') & (R.group == 'ALL_INDEPENDENT') & (R.hprc_reference == 'ALL')]
    for t, col in TIER.items():
        a = A[A.tier == t]; ax.scatter(a.hprc_freq, a.freq, s=8, color=col, alpha=0.7, lw=0, label=f'tier {t} (n={len(a)})')
    ax.plot([0, 1], [0, 1], color='#999999', lw=0.7, zorder=0)
    for aid, xy in zip(KEY_ALLELES[:3], ((0.12, 0.74), (0.60, 0.33), (0.05, 0.47))):
        a = A[A.item == aid]
        if len(a): ax.annotate(f"{aid} ({a.form.iloc[0]})", (a.hprc_freq.iloc[0], a.freq.iloc[0]), xytext=xy, textcoords='data', fontsize=6,
                               ha='left', arrowprops=dict(arrowstyle='-', lw=0.5, color='#444444', shrinkA=1, shrinkB=2))
    fc = SC['freq_correlation_non_native']
    ax.text(0.03, 0.95, f"non-native alleles: Pearson r = {fc['pearson_r']:.2f} (n={fc['n_alleles']})", transform=ax.transAxes, fontsize=7, va='top')
    ax.set_xlabel('HPRC Release 2 frequency (arm-specific)'); ax.set_ylabel('frequency in HGSVC3, CPC, APC, JASA')
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02); ax.legend(frameon=False, loc='lower right')
    ax.set_title('Allele frequencies replicate in independent cohorts', loc='left')
    ax.text(-0.32, 1.06, 'c', transform=ax.transAxes, fontsize=10, fontweight='bold')
    # (d) key alleles per cohort
    ax = fig.add_subplot(gs[1, 1]); keys = list(KEY_ALLELES[:4])
    grp = [('HGSVC3', '#6a3d9a'), ('CPC', '#33a02c'), ('APC', '#b15928'), ('JASA', '#e7298a')]
    for i, aid in enumerate(keys):
        h = R[(R.item == aid) & (R.group == 'HGSVC3') & (R.hprc_reference == 'ALL')]
        if len(h): ax.plot([h.hprc_freq.iloc[0]] * 2, [i - 0.38, i + 0.38], color='black', lw=1.2, zorder=1)
        for j, (g, col) in enumerate(grp):
            a = R[(R.item == aid) & (R.group == g) & (R.hprc_reference == 'ALL')]
            if not len(a): continue
            yy = i - 0.27 + 0.18 * j
            ax.errorbar(a.freq.iloc[0], yy, xerr=[[a.freq.iloc[0] - a.ci95_lo.iloc[0]], [a.ci95_hi.iloc[0] - a.freq.iloc[0]]], fmt='o', ms=3,
                        color=col, lw=0.8, label=g if i == 0 else None)
    ax.set_yticks(range(len(keys))); lab = {r.item: f"{r.item}\n{r.form}" for r in R[R.kind == 'structural_allele'].drop_duplicates('item').itertuples()}
    ax.set_yticklabels([lab.get(k, k) for k in keys]); ax.invert_yaxis(); ax.set_xlim(0, 1.0)
    ax.set_xlabel('frequency (95% CI); black bar = HPRC R2')
    ax.legend(frameon=False, loc='lower right', fontsize=6)
    ax.set_title('Common alleles in every cohort', loc='left')
    ax.text(-0.42, 1.06, 'd', transform=ax.transAxes, fontsize=10, fontweight='bold')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{T.FIG}/tan2026_concordance.{ext}', dpi=300, bbox_inches='tight')
    print('written', f'{T.VAL}/tan2026_summary.json', f'{T.FIG}/tan2026_concordance.png')

if __name__ == '__main__':
    main()
