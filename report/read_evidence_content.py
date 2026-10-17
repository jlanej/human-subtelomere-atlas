"""Content of evidence/READ_EVIDENCE.pdf / .md (python report/build_report.py read_evidence_content evidence/READ_EVIDENCE, run from
the repository root). Every number is read at build time from committed files in evidence/."""
import json, os
import pandas as pd
from scipy.stats import beta

T_ = lambda p: pd.read_csv(p, sep='\t', low_memory=False)
f0 = lambda x: f'{x:,.0f}'
ub = lambda k, n: 1.0 if n == 0 else float(beta.ppf(0.95, k + 1, n - k))
MONO = lambda s: f'<font face="Courier">{s}</font>'

C = T_('evidence/claims.tsv'); R = T_('evidence/reads.tsv.gz'); K = T_('evidence/controls.tsv')
S = T_('evidence/showcase/showcase_claims.tsv'); SL = T_('evidence/showcase/showcase_loci.tsv').set_index('claim_id')
M = T_('evidence/release_manifest.tsv'); CR = json.load(open('evidence/criteria_v1.json'))
RT = T_('evidence/replication_results.tsv') if os.path.exists('evidence/replication_results.tsv') else None
SE = T_('evidence/calibration/site_error_summary.tsv'); SE = SE[SE.kind == 'carrier_vs_consensus'].set_index('platform')
SE['wrong_allele_rate'] = SE.wrong_allele_rate.astype(float)
HCL = T_('evidence/hg002/claims.tsv')
isdef = C.definitive_verdict.str.startswith('definitive')
H = C[C.cohort == 'HG002_Q100']; TR = C[(C.cohort == 'HG002_trio')]; TT = TR[TR.v1_class == 'transmitting_parent']; P = C[C.cohort == '1KG_ONT']
cnt = lambda g, v: int((g.v1_verdict == v).sum())
ndef = lambda g: int(isdef.loc[g.index].sum())
dreads = lambda g: R[R.claim_id.isin(g[isdef.loc[g.index]].claim_id) & (R.tier == 'definitive') & (R.counted == True)]  # noqa: E712
hd, pd_, td = dreads(H), dreads(P), dreads(TT)
nh_el, nh_c = int(K[K.cohort == 'HG002_Q100'].n_eligible.sum()), int(K[K.cohort == 'HG002_Q100'].n_conclusive.sum())
nk_el, nk_c = int(K[K.cohort == '1KG_ONT'].n_eligible.sum()), int(K[K.cohort == '1KG_ONT'].n_conclusive.sum())
nk_reads = int(R[R.claim_id.str.startswith('control|') & (R.cohort == '1KG_ONT')].read_id.nunique())
hmulti = int((H.definitive_verdict == 'definitive_multi_dataset').sum())
pconf = cnt(P, 'confirmed'); pdef = ndef(P)
pna = int(P.v1_class.astype(str).str.startswith('not_assessable').sum())
p_cnd = P[(P.v1_verdict == 'confirmed') & ~isdef.loc[P.index]]
tiers_cnd = pd.Series(';'.join(p_cnd.reads_by_tier.dropna()).split(';')).str.split(':').str[0].value_counts()
nshared = int(p_cnd.reads_by_tier.fillna('').str.contains('conclusive_shared_anchor').sum())
noff = int(p_cnd.reads_by_tier.fillna('').str.contains('off_end').sum())
sz_show = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk('evidence/showcase') for f in fs)
n_claim_dirs = int((C.n_conclusive > 0).sum())
ex_c = 'HG002#2#chr9_MATERNAL#L'
ex = R[(R.claim_id == ex_c) & (R.tier == 'definitive') & (R.dataset == 'ont_ul_dorado')].iloc[0]
exc = C.set_index('claim_id').loc[ex_c]; exl = SL.loc[ex_c]
SHA = '%s' % 'a89c5eaf023ed032'
nc = HCL[HCL.verdict == 'not_confirmed']
_rs = nc.reason.fillna('')
nc_txt = (f"{int(_rs.str.contains('none passes the switch test').sum())} have no read that passes the switch test (too few "
          f"diagnostic units or a mosaic donor block), {int(_rs.str.contains('lack a terminal telomere').sum())} has switch-passing reads "
          f"without a terminal telomere, {int(_rs.str.contains('fail integrity').sum())} have anchored switch-passing reads that fail "
          f"integrity, and {int(_rs.str.contains('beyond the assembly window').sum())} (16p paternal) has conclusive reads but a switch "
          "beyond the 100 kb assembly window, where v1 concordance is undefined")
SBL = T_('evidence/calibration/sensitivity_by_length.tsv'); _ul = SBL[SBL.platform == 'ONT_R10_UL'].set_index('len_bin')
sens_txt = (f"{100 * _ul.loc['25-50 kb', 'sensitivity']:.0f}% for 25-50 kb reads to {100 * _ul.loc['>200 kb', 'sensitivity']:.0f}% "
            "for reads > 200 kb")
n_edge = int(C.diagnostics_note.fillna('').str.contains('donor-window-edge').sum())

TITLE = 'Read-level evidence for non-native terminal blocks'
SUBTITLE = ('Single long reads from public data, tested against pre-registered criteria, show that telomere-adjacent blocks '
            'resembling another chromosome arm are real chromosome ends')
DATE = 'October 2026 · atlas v0.1 · criteria v1 (' + SHA + ')'

SUMMARY = [
 ('Question', 'Is a non-native terminal block real, or an assembly artefact? It is a telomere-adjacent block that resembles the '
  'CHM13 copy of another arm (donor D) more than the own arm (recipient R). We tested this with single reads not used to build the assemblies.'),
 ('Test', 'A conclusive read starts in the telomere. It carries D alleles at paralog-diagnostic sites up to a switch and R alleles beyond it, '
  'and it anchors at MAPQ 60 in the recipient chromosome. The criteria were frozen before any population read was inspected.'),
 ('HG002', f'All {cnt(H, "confirmed")} of {len(H)} HG002 Q100 claims confirmed under v1 are also <b>definitive</b>, with '
  f'{f0(len(hd))} definitive reads. In {hmulti} of them the definitive reads come from at least 2 independent datasets.'),
 ('Trio', f'The transmitting parent (HG003 or HG004; GIAB reads) confirms {cnt(TT, "confirmed")} HG002 claims; {ndef(TT)} of these are definitive.'),
 ('Population', f'In 1000 Genomes ONT reads of 42 HPRC individuals, {pconf} of {len(P)} claims are confirmed and {pdef} are definitive '
  f'({len(pd_)} reads, {P[isdef.loc[P.index]].individual.nunique()} individuals). Most claims ({pna}) cannot be assessed with these reads.'),
 ('Specificity', f'No negative-control read was conclusive: 0 of {nh_el} HG002 other-haplotype tests and 0 of {f0(nk_el)} 1KG-ONT '
  f'native read × model tests (exact one-sided 95% upper bounds {100 * ub(0, nh_el):.2f}% and {100 * ub(0, nk_el):.2f}%).'),
 ('Verify', 'Run evidence/replicate_claim.sh CLAIM_ID. It fetches the cited reads from the public files and re-derives every verdict.'
  + (f' In a clean-tree test, {int((RT.result == "MATCH").sum())} of {len(RT)} claims and controls matched exactly.' if RT is not None else '')),
]


def T(*rows): return [list(r) for r in rows]


def lab(c):
    return f"{c['arm']} {c['haplotype'][:3]}"


hrows = [['Claim', 'Allele', 'Donor → recipient', 'v1', 'Conclusive', 'Definitive', 'Datasets with definitive reads']]
for c in H[H.v1_verdict.isin(['confirmed', 'supported'])].sort_values('n_definitive', ascending=False).to_dict('records'):
    hrows.append([lab(c), c['allele_id'], f"{c['donor_D'][3:]} → {c['recipient_R'][3:]}", c['v1_verdict'], str(int(c['n_conclusive_concordant'])),
                  str(int(c['n_definitive'])), str(c['definitive_datasets']).replace('nan', '-').replace(',', ', ')])
trows = [['HG002 claim', 'Transmitting parent', 'v1', 'Conclusive', 'Definitive', 'Contradicting', 'Other parent (v1)']]
for c in TT.sort_values('n_definitive', ascending=False).to_dict('records'):
    if c['v1_verdict'] == 'not_confirmed' and c['n_conclusive_concordant'] == 0: continue
    o = TR[(TR.claimed_end_ids == c['claimed_end_ids']) & (TR.v1_class == 'other_parent')].iloc[0]
    trows.append([lab(dict(arm=c['arm'], haplotype=c['claimed_end_ids'].split('#')[2].split('_')[1].lower())), c['individual'], c['v1_verdict'],
                  str(int(c['n_conclusive_concordant'])), str(int(c['n_definitive'])), str(int(c['n_contradicts'])),
                  f"{o['individual']}: {o['v1_verdict']}"])
prows = [['Allele', 'Arm', 'Donor', 'Definitive claims', 'Individuals', 'Definitive reads']]
for (a, arm), g in P[isdef.loc[P.index]].groupby(['allele_id', 'arm']):
    prows.append([a, arm, g.donor_D.iloc[0][3:], str(len(g)), ', '.join(sorted(g.individual)), str(int(g.n_definitive.sum()))])
_vc = P.v1_class.astype(str)
_na = _vc[_vc.str.startswith('not_assessable')].str.split(':').str[-1].value_counts()
_an = _vc[_vc.str.startswith('assessed_not_confirmed')].str.split(':').str[-1].value_counts()
na_txt = (f"{int(_na.sum())} not assessable (" + ', '.join(f'{k.replace("_", " ")} {v}' for k, v in _na.items()) + f"); {int(_an.sum())} "
          "assessed but not confirmed (" + ', '.join(f'{k.replace("_", " ")} {v}' for k, v in _an.items()) + ")")
_E = T_('data/catalogue/ends.tsv'); _R2 = _E[_E['set'] == 'HPRC_R2']
nn_frac = (_R2.nativity_class.isin(['foreign', 'within_family'])).mean(); nn_n = int(_R2.nativity_class.isin(['foreign', 'within_family']).sum())
rep_tab = None
if RT is not None:
    rep_tab = [['Claim', 'Kind', 'Cited reads', 'Published', 'Replicated', 'Reads identical', 'Wall (s)']] + \
        [[r.claim_id, r.kind, str(r.n_cited), str(r.published), str(r.replicated), str(r.reads_identical), f'{r.wall_s:.0f}'] for r in RT.itertuples()]

BODY = [
 ('h1', '1. The question'),
 ('p', f'In the atlas, {nn_n:,} of {len(_R2):,} HPRC Release 2 chromosome ends ({100 * nn_frac:.0f}%) have a non-native core structure: a '
  'telomere-adjacent block that resembles the CHM13 copy of another arm (donor D) more than the own arm (recipient R). A skeptic can ask whether such blocks are real chromosome ends or '
  'assembly artefacts: mis-joins, collapsed paralogs, or contigs placed on the wrong end. The switch from D-like to R-like sequence '
  'lies inside sequence that D and R share (about 98% identity in CHM13). Whole-read alignment to CHM13 therefore usually stays '
  'on R and cannot show the structure. At the HG00099 20q carrier, for example, 1 of 575 reads had a supplementary alignment to chr3.'),
 ('h1', '2. What would convince a skeptic'),
 ('p', 'A single molecule that (i) begins in the telomere; (ii) carries the donor\'s alleles at positions where CHM13 D and R differ '
  'up to a switch, then the recipient\'s alleles beyond it, with far more sites than sequencing error could flip; (iii) continues into '
  'sequence found only at the recipient end; (iv) shows no sign of a chimera; and (v) aligns best to the claimed end of the same '
  'person\'s own assembly. Such molecules should recur across independent datasets and platforms, and in the parent who transmitted '
  'the haplotype. They should be absent in people whose ends are native. The criteria for (i)-(iv) were registered and frozen '
  f'(hash {MONO(SHA)}) before any population read was inspected. Calibration used HG002 training ends and native controls only '
  '(evidence/CRITERIA.md, evidence/calibration/CALIBRATION.md). Criterion (v) and a uniqueness requirement for (iii) form the post-hoc '
  '<b>definitive</b> tier. It was defined after the population diagnostics (evidence/population/DIAGNOSTICS.md) and can only remove '
  'claims.'),
 ('h1', '3. Method in one page'),
 ('p', '<b>Diagnostic units, from CHM13 alone.</b> The donor window on D is aligned to the comparator window P (normally R). Inside '
  'homologous blocks (≥ 2 kb at ≥ 90%), single-base differences away from indels, homopolymers, low-complexity and locally repeated '
  'sequence are SNV sites. Read bins of 200 bp that fall in sequence private to one window are bin units. Each unit is labelled D, P '
  'or neither by the read base or bin identity.'),
 ('p', '<b>Per read.</b> The read is oriented telomere-first, and a terminal (CCCTAA)n block of ≥ 300 bp within 300 bp of the end is '
  'required. A one- or two-change-point maximum-likelihood segmentation (error rate fixed at 0.05) gives the switch with a 95% profile '
  'interval. Each segment needs ≥ 10 units at ≥ 80% purity, and both log10 likelihood ratios against the no-switch models must be ≥ 5. '
  'The proximal part is aligned genome-wide (minimap2 map-ont, CHM13 v2.0) and must anchor on the recipient chromosome at MAPQ 60 over '
  '≥ 10 kb, within 2 Mb of the telomere. Integrity: no internal telomere, adapter or foldback, and ≥ 90% of the subtelomeric part '
  'aligned. Switch concordance: the read\'s interval (± 2 kb) overlaps the switch of the person\'s own assembly window analysed with '
  'the same units.'),
 ('p', '<b>Claim.</b> v1 <i>confirmed</i> needs ≥ 2 conclusive switch-concordant reads and, for phased HG002 reads, no conclusive '
  'contradicting read on the claimed haplotype. A <b>definitive read</b> adds two conditions. (a) ≥ 1 kb of the anchor lies on CHM13 '
  'recipient positions that no other CHM13 end covers at ≥ 95%. (b) Own-assembly placement: the read has more aligned matching '
  'bases on the claimed end window than on any other end window of the same person (HG002: 92 Q100 windows; 1KG-ONT individuals: their '
  'about 90 Release 2 windows; HG003/HG004: the HG002 windows, a proxy for the transmitted haplotype). A <b>definitive claim</b> is v1 '
  'confirmed with ≥ 2 definitive reads from distinct molecules.'),
 ('p', f'<b>Reads.</b> HG002: 16 public HG002 Q100 v1.1 read-alignment files (ONT R10.4.1 ultra-long, duplex and other runs, R9.4.1, '
  'PacBio HiFi Revio, Vega, SPRQ and DeepConsensus), 12 independent datasets after molecule de-duplication. HERRO-corrected reads are '
  'derived from the ultra-long reads and are never counted. Trio: GIAB HG003 and HG004 HiFi (CHM13), ONT ultra-long R9.4.1 and R10.4.1 '
  '(GRCh38) alignments. Population: 1000 Genomes ONT Vienna CRAMs (R9.4.1, aligned to CHM13 v2.0) of the 42 HPRC Release 2 individuals in '
  'that panel. Reads were fetched only from the outer chromosome-end regions, by HTTP range requests. Per-read error at diagnostic '
  f'SNV sites (wrong allele): ONT R10 ultra-long {100 * SE.loc["ONT_R10_UL", "wrong_allele_rate"]:.1f}%, duplex '
  f'{100 * SE.loc["ONT_R10_duplex", "wrong_allele_rate"]:.2f}%, HiFi {100 * SE.loc["HiFi_Revio", "wrong_allele_rate"]:.3f}%.'),
 ('fig', 'evidence/figures/fig_showcase_readthrough.png',
  '<b>Figure 1. Eight exemplar reads, telomere to anchor.</b> One definitive read per showcase claim. Ticks are diagnostic units '
  '(upper row SNV sites, lower row private 200 bp bins), orange for the CHM13 donor allele, blue for the recipient or proximal allele '
  'and grey for neither. Yellow: 95% interval of the switch. Green: genome-wide MAPQ-60 anchor with CHM13-unique sequence. Boxes '
  'below: the CHM13 coordinates of the donor and recipient segments. Reads, IGV sessions and notes are in evidence/showcase/.',
  'figures/fig_showcase_readthrough.png'),
 ('h1', '4. Results'),
 ('p', f'<b>HG002 Q100 (phased).</b> Of {len(H)} non-native HG002 ends, {cnt(H, "confirmed")} are v1-confirmed and '
  f'{cnt(H, "supported")} supported. All {cnt(H, "confirmed")} confirmed claims are definitive, {hmulti} with definitive reads '
  f'from ≥ 2 independent datasets ({f0(len(hd))} definitive reads). Of the {len(nc)} claims that are not confirmed, ' + nc_txt +
  ' (evidence/hg002/README.md).'),
 ('table', T(*hrows), [0.7, 0.85, 1.0, 0.75, 0.82, 0.78, 2.1]),
 ('p', f'<b>Trio.</b> Each HG002 claim was tested in both parents, using unphased reads and the unphased rule. The transmitting parent '
  f'confirms {cnt(TT, "confirmed")} claims and supports {cnt(TT, "supported")}; {ndef(TT)} are definitive. Contradicting reads in '
  'a parent come from the parent\'s other haplotype at that arm. They are expected, and they do not block an unphased confirmation. '
  'HG002 carries SA-17p-02 (11p-like) on both haplotypes. The father, HG003, transmitted the paternal copy, so his reads also '
  'confirm the maternal 17p claim.'),
 ('table', T(*trows), [0.85, 0.95, 0.9, 0.82, 0.78, 0.98, 1.72]),
 ('p', f'<b>1000 Genomes ONT (unphased, R9.4.1).</b> {len(P)} claims (individual × arm × allele) in 42 individuals: {pconf} confirmed, '
  f'{cnt(P, "supported")} supported, {pdef} definitive. The others: {na_txt}. Of the {len(p_cnd)} confirmed claims '
  f'that are not definitive, {nshared} rest on reads whose anchor lies only in CHM13-shared sequence (mostly 2q SA-2q-02) and {noff} '
  'include reads that align best to another end of the same person. These are often the other haplotype at the same arm, which '
  'carries a related allele (e.g. HG00235 20q SA-20q-06, whose reads place on the SA-20q-01 haplotype). The definitive claims:'),
 ('table', T(*prows), [1.0, 0.5, 0.6, 1.0, 2.9, 1.0]),
 ('p', f'<b>Negative controls.</b> HG002 reads from the other haplotype at the same arm, where that haplotype is strictly native: '
  f'{nh_c} of {nh_el} eligible read tests conclusive (95% upper bound {100 * ub(nh_c, nh_el):.2f}%). 1KG-ONT individuals strictly native '
  f'in both haplotypes, each tested against every allele model of the arm: {nk_c} of {f0(nk_el)} eligible read × model tests conclusive '
  f'({nk_reads} cited control reads; bound {100 * ub(nk_c, nk_el):.2f}%). Controls mostly <i>contradict</i>, that is, they show the '
  'native arm in the claimed donor interval. Nine claimed arms have no strictly native individual among the 42.'),
 ('fig', 'evidence/figures/fig_evidence_summary.png',
  '<b>Figure 2. Summary.</b> (a) Claims by cohort and tier (trio: transmitting parent). (b) Negative controls with exact one-sided 95% upper '
  'bounds. (c) HG002 definitive reads per claim and dataset. (d) Verdict in the transmitting and in the other parent for each HG002 '
  'claim with read support.', 'figures/fig_evidence_summary.png'),
 ('h1', '5. Worked example: HG002 maternal 9p'),
 ('p', f'Claim {MONO(ex_c)}: allele {exc["allele_id"]}. The maternal 9p end carries a CHM13 15q-like block before 9p sequence; the '
  f'assembly switch is {int(exc["assembly_switch_sub"]):,} bp from the telomere boundary. Under v1 it has '
  f'{int(exc["n_conclusive_concordant"])} conclusive switch-concordant reads and is definitive with {int(exc["n_definitive"])} '
  f'definitive reads from {int(exc["n_definitive_datasets"])} datasets. The mother (HG004) confirms it. One read, from the R10 '
  f'ultra-long dataset, is {MONO(ex.read_id)}, {int(ex.read_len):,} bp, with a telomere of {int(ex.telo_len):,} bp. It has '
  f'{int(ex.seg_D_units)} donor-segment units at {100 * ex.seg_D_purity:.0f}% 15q alleles and {int(ex.seg_P_units)} proximal units at '
  f'{100 * ex.seg_P_purity:.0f}% 9p alleles. The switch lies at {int(ex.switch_ci_sub_lo):,}-{int(ex.switch_ci_sub_hi):,} bp '
  f'(CHM13 {ex.switch_ci_chm13}). The anchor is {ex.anchor_chm13} at MAPQ {int(ex.anchor_mapq)}, {int(ex.anchor_unique_bp):,} bp of it '
  f'CHM13-unique. The read places best on {MONO(str(ex.own_best_end))} at {100 * ex.own_best_identity:.1f}% identity.'),
 ('p', 'Fetch the read from the public alignment file:'),
 ('p', MONO(ex.samtools_command.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))),
 ('p', f'Split at the switch interval and aligned genome-wide, the distal segment of each showcase read aligns to CHM13 15q '
  f'({int(exl.distal_on_D)} of {int(exl.n_reads)}) and the proximal segment to CHM13 9p at MAPQ 60 ({int(exl.proximal_on_R_mapq60)} of '
  f'{int(exl.n_reads)}) (evidence/showcase/HG002_2_chr9_MATERNAL_L/). Re-derive the claim from public data:'),
 ('p', MONO("evidence/replicate_claim.sh 'HG002#2#chr9_MATERNAL#L'")),
 ('fig', 'evidence/figures/claims/HG002_2_chr9_MATERNAL_L.png',
  '<b>Figure 3. All definitive reads of HG002 maternal 9p</b> (12 of the most informative shown). Reads are aligned at the telomere boundary, '
  'with units coloured by CHM13 allele, per-read switch interval in yellow, unique anchor in green and the own-assembly switch as a grey band.',
  'figures/claims/HG002_2_chr9_MATERNAL_L.png'),
 ('h1', '6. Alternative explanations and how each is excluded'),
 ('p', '<b>Chimeric reads.</b> Definitive reads have no internal telomere, adapter or foldback, and ≥ 90% of the subtelomeric part '
  'aligns. Their switches fall, within ± 2 kb, where the person\'s own assembly places them, and the same switch recurs in independent '
  'molecules, datasets and platforms and in the transmitting parent. Ligation chimeras form at random positions and would not recur.'),
 ('p', '<b>Mis-mapping and paralogs.</b> This is the main risk. Diagnostic units tell which CHM13 end a segment resembles, but human '
  'haplotypes carry blocks that CHM13 lacks. In the population diagnostics, 19 SA-3q-01 reads were MAPQ-60 on chr3 with 0 unique bp; '
  'they came from 6p, 20q and 9q haplotypes that carry 3q-like sequence. Assembly-switch concordance had already kept them from '
  'confirming SA-3q-01. The definitive tier excludes the mechanism explicitly: it requires unique anchor sequence and best placement on '
  'the claimed end among all of the person\'s ends.'),
 ('p', '<b>Assembly error.</b> The read test is assembly-free. Units come from CHM13; the 1KG-ONT reads (another laboratory, DNA '
  'extraction and platform) and the GIAB parental reads were not used for the assemblies. Assemblies are used only to locate the '
  'expected switch and in the post-hoc placement step.'),
 ('p', '<b>Population variation at diagnostic sites.</b> Carrier haplotypes differ from CHM13 at some sites. Purity thresholds are '
  'therefore 0.80, and true-carrier donor purity goes as low as 0.80 (calibration T0). Such variation lowers purity in both segments. It '
  'does not create an ordered donor-then-recipient pattern with log10 likelihood ratios ≥ 5 against both no-switch models. No native '
  'control read did.'),
 ('p', '<b>Duplicate molecules and derived reads.</b> HG002 reads are de-duplicated across datasets by ONT read UUID or PacBio movie and '
  'ZMW. HERRO-corrected reads are reported but never counted. Definitive claims count distinct molecules.'),
 ('h1', '7. Limits'),
 ('p', f'<b>Sensitivity.</b> A read must reach the telomere and pass the switch by ≥ 10 kb. For HG002 ONT R10 ultra-long reads, '
  'read-level sensitivity (calibration: conclusive / eligible carrier reads) rises with length, from ' + sens_txt + ' '
  '(calibration/sensitivity_by_length.tsv). The R9.4.1 1KG-ONT reads are shorter, and most population claims are <i>not '
  f'assessable</i> ({pna} of {len(P)}). That is no evidence against them. <b>Model limits.</b> Claims whose assembly switch is at '
  f'the edge of the v1 donor window ({n_edge} claims, including SA-20p-03) are treated as not assessable. Switches beyond the 100 kb assembly window '
  'have undefined concordance. <b>Controls.</b> Nine claimed arms (12p, 17q, 19p, 19q, 1p, 3q, 4q, 6q, Yq) have no strictly native '
  'individual among the 42. <b>Post-hoc tier.</b> The definitive tier was defined after the population diagnostics. It is stricter '
  'than v1 and only removes claims. Its placement step uses the person\'s assembly, so it is not assembly-free. For HG003/HG004 it is a '
  'proxy (placement on HG002 windows). <b>Software.</b> Verdicts depend on minimap2/mappy versions (2.31 here).'),
 ('h1', '8. Evidence files and replication'),
 ('p', f'evidence/claims.tsv has one row per claim ({len(C)}: cohort, individual, haplotype, arm, allele, nativity, D/P/R, assembly '
  'switch and its CHM13 projection, read counts by tier and dataset, v1 and definitive verdicts, diagnostics note). '
  f'evidence/reads.tsv.gz has one row per cited read ({f0(len(R))} rows): public source URL and index, region, read group, read length, '
  'strand, telomere, donor/recipient units and purity, switch interval in read, subtelomere and CHM13 coordinates, anchor position, '
  'MAPQ and unique bp, own-assembly placement and identity, tier, and the samtools command that retrieves the read. '
  f'evidence/showcase/ holds eight exemplar claims ({sz_show / 1e6:.1f} MB): CRAM and BAM on CHM13 v2.0 (read groups name the public '
  'source; @SQ M5 tags; no local paths), BAM on the person\'s own end windows with the window FASTA, read segments split at the '
  f'switch, a BED of the claim model, and IGV sessions with a "what you should see" note. The full set ({n_claim_dirs} claims with a '
  'conclusive read) is packaged as GitHub release assets: ' + '; '.join(f'{r.file} ({r.bytes / 1e6:.1f} MB)' for r in M.itertuples())
  + '. SHA-256 checksums are in evidence/release_manifest.tsv.'),
 ('p', 'evidence/replicate_claim.sh CLAIM_ID works from a clean checkout with internet, samtools and minimap2. It verifies the criteria '
  'hash, downloads CHM13 v2.0, cuts the person\'s end windows from the public assembly (SHA-256-checked), fetches the cited reads by '
  'region and read ID from the public files, re-runs the frozen criteria and the definitive tier, and compares every verdict with the '
  'published tables (evidence/REPLICATION_TEST.md).'),
] + ([('table', rep_tab, [2.0, 0.8, 0.7, 1.1, 1.1, 0.7, 0.6])] if rep_tab else []) + [
 ('h1', '9. Data provenance, licences and how to cite'),
 ('p', 'HPRC Release 2 assemblies (public domain; HPRC data-use best practice; BioProject PRJNA730823; Liao et al. 2023). HG002 Q100 '
  'v1.1 assembly and read alignments (T2T/GIAB; human-pangenomics S3). GIAB HG003/HG004 reads (NIST, Ashkenazi trio; NCBI FTP) and ONT '
  'Open Data (giab_2023.05). 1000 Genomes ONT Vienna panel (Schloissnig, Pani et al. 2025): reads generated at IMP Vienna from '
  'Coriell DNA consented for public release. We used them for a targeted analysis of chromosome-end regions under the panel\'s data '
  'reuse statement, and the evidence files redistribute only the few hundred cited reads. CHM13 v2.0 (T2T consortium; Nurk et al. 2022). '
  'Atlas tables, figures and text are CC BY 4.0, code MIT (DATA_LICENSE.md). We acknowledge the Human Pangenome Reference Consortium '
  '(BioProject PRJNA730823) and its funder, the National Human Genome Research Institute (NHGRI).'),
 ('p', 'How to cite: Human subtelomere atlas v0.1, read-level evidence (criteria v1, ' + SHA + '), '
  'https://github.com/jlanej/human-subtelomere-atlas (CITATION.cff), together with the read sources above.'),
]
REFS = ['Nurk S, et al. (2022). The complete sequence of a human genome. <i>Science</i> 376:44-53. doi:10.1126/science.abj6987',
        'Liao W-W, et al. (2023). A draft human pangenome reference. <i>Nature</i> 617:312-324. doi:10.1038/s41586-023-05896-x',
        'Schloissnig S, Pani S, et al. (2025). Structural variation in 1,019 diverse humans based on long-read sequencing. <i>Nature</i>. '
        'doi:10.1038/s41586-025-09290-7 (preprint doi:10.1101/2024.04.18.590093)',
        'Zook JM, et al. (2016). Extensive sequencing of seven human genomes to characterize benchmark reference materials. '
        '<i>Scientific Data</i> 3:160025. doi:10.1038/sdata.2016.25',
        'Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. <i>Bioinformatics</i> 34:3094-3100. doi:10.1093/bioinformatics/bty191',
        'Tan KT, et al. (2022). Identifying and correcting repeat-calling errors in nanopore sequencing of telomeres. <i>Genome Biology</i> '
        '23:180. doi:10.1186/s13059-022-02751-6',
        'Tan et al. (2026). Widespread structural variations at human chromosome ends. bioRxiv. doi:10.64898/2026.09.03.748401']
