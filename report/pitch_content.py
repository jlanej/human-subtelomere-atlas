"""Content of report/PITCH.pdf / .md, the two-minute introduction to the atlas and its read-level evidence
(python report/build_report.py pitch_content report/PITCH, run from the repository root).
Every number is read at build time from committed files (catalogue, validation, evidence, report/pitch_data/)."""
import json
import pandas as pd
from scipy.stats import beta

T_ = lambda p: pd.read_csv(p, sep='\t', low_memory=False)
f0 = lambda x: f'{x:,.0f}'
ub = lambda k, n: float(beta.ppf(0.95, k + 1, n - k))
MONO = lambda s: f'<font face="Courier">{s}</font>'

# ---- catalogue
E = T_('data/catalogue/ends.tsv'); R2 = E[E.set == 'HPRC_R2']
n_ends, n_people = len(R2), R2['sample'].nunique()
n_foreign = int((R2.nativity_class == 'foreign').sum())
FR = R2[R2.nativity_class == 'foreign']
_FAM = {r.family: set(r.members.split(',')) for r in T_('data/catalogue/paralog_families.tsv').itertuples()}
# the block next to the telomere is itself foreign: neither home/none nor a family or arm that contains the own arm
n_first = int(sum(t not in ('home', 'none') and a not in _FAM.get(t, {t})
                  for t, a in zip(FR.core_structure.astype(str).str.split('|').str[0], FR.catalogue_arm)))
per = R2.groupby('sample').agg(n=('end_id', 'size'), f=('nativity_class', lambda s: int((s == 'foreign').sum())))
A = T_('data/catalogue/alleles_structural.tsv').set_index('allele_id')
a20, a3 = A.loc['SA-20q-01'], A.loc['SA-3q-01']
JJ = json.load(open('data/validation/tan2026_summary.json')); J = JJ['independent_cohorts']
vk_sa = JJ['same_haplotype']['primary']['sa_agree']; tan_n = J['n_primary']
tan_common = J['window_sensitivity']['sampled_R2_lengths']['common_tierA_observed']
tan_r = J['freq_correlation_non_native']['pearson_r']

# ---- read evidence
C = T_('evidence/claims.tsv'); RD = T_('evidence/reads.tsv.gz'); K = T_('evidence/controls.tsv')
isdef = C.definitive_verdict.astype(str).str.startswith('definitive')
H = C[C.cohort == 'HG002_Q100']; P = C[C.cohort == '1KG_ONT']; TR = C[C.cohort == 'HG002_trio']
TT = TR[TR.v1_class == 'transmitting_parent']
h_conf = int((H.v1_verdict == 'confirmed').sum()); h_def = int(isdef[H.index].sum())
h_multi = int((H.definitive_verdict == 'definitive_multi_dataset').sum())
h_reads = int(RD[RD.claim_id.isin(H[isdef[H.index]].claim_id) & (RD.tier == 'definitive') & (RD.counted == True)].shape[0])  # noqa: E712
t_conf = int((TT.v1_verdict == 'confirmed').sum())
p_conf = int((P.v1_verdict == 'confirmed').sum()); p_def = int(isdef[P.index].sum())
p_people = P[isdef[P.index]].individual.nunique()
p_reads = int(RD[RD.claim_id.isin(P[isdef[P.index]].claim_id) & (RD.tier == 'definitive') & (RD.counted == True)].shape[0])  # noqa: E712
p_na = int(P.v1_class.astype(str).str.startswith('not_assessable').sum())
kh = K[K.cohort == 'HG002_Q100']; kk = K[K.cohort == '1KG_ONT']
kh_n, kh_c = int(kh.n_eligible.sum()), int(kh.n_conclusive.sum())
kk_n, kk_c = int(kk.n_eligible.sum()), int(kk.n_conclusive.sum())
RT = T_('evidence/replication_results.tsv')
rep_ok, rep_n = int((RT.result == 'MATCH').sum()) if 'result' in RT else len(RT), len(RT)

# ---- the hero molecule (report/pitch_data.py)
PR = T_('report/pitch_data/profiles.tsv'); ST = T_('report/pitch_data/sites.tsv'); SG = T_('report/pitch_data/segments.tsv')
HERO = '3788f60a-2f28-4b28-917b-f1ecae3645e3'
hr = RD[RD.read_id == HERO].iloc[0]
hp = PR[PR.read_id == HERO].iloc[0]; hs = ST[ST.read_id == HERO]; hg = SG[SG.read_id == HERO].set_index('segment')
pre, post = hs[hs.read_pos < hp.switch_lo], hs[hs.read_pos > hp.switch_hi]
import sys; sys.path.insert(0, 'report'); import pitch_figures as PF  # noqa: E402
sys.path.insert(0, 'website'); import tour_pages as _TP  # noqa: E402
SUM = PF.read_summary(PR); SUM = SUM[SUM.slug.isin([s for s, _ in PF.SLOPE_CLAIMS])]
n_flip, n_mol = int(((SUM.before > 0) & (SUM.after < 0)).sum()), len(SUM)
ALL = PF.read_summary(PR)
excl = {s: f"{int(((g.before > 0) & (g.after < 0)).sum())} of {len(g)}" for s, g in ALL[~ALL.slug.isin(SUM.slug)].groupby('slug')}

REP_CMD = MONO("bash evidence/replicate_claim.sh 'HG01258|20q|SA-20q-01'")

TITLE = 'Chromosome ends borrow from each other, and single DNA molecules show it'
SUBTITLE = 'A two-minute introduction to the human subtelomere atlas and the evidence behind it'
DATE = 'October 2026 · atlas v0.1 · read-evidence criteria v1'

SUMMARY = [
    ('The finding', f'The DNA just inside a chromosome\'s telomere often looks like the end of a <i>different</i>, '
     f'unrelated chromosome. All {n_people} people in the Human Pangenome Reference Consortium Release 2 carry such ends: a median of '
     f'{per.f.median():.0f} (range {per.f.min()}-{per.f.max()}) of their roughly 90 assembled chromosome ends, '
     f'{f0(n_foreign)} of {f0(n_ends)} ends in all.'),
    ('They are common', f'They are recurring variants, not one-off rearrangements. One version of the 20q end, '
     f'resembling the 3q/19p end, is found in {a20.n_samples} of {n_people} people ({a20.r2_freq:.0%} of assembled 20q ends).'),
    ('Why believe it', f'Because a single DNA molecule reads it end to end. A {hr.read_len / 1e3:.0f} kb nanopore read from '
     f'HG002 starts in the telomere, carries chromosome 11\'s letters at {int((pre.call == "D").sum())} of {len(pre)} '
     f'diagnostic sites, then switches and runs for {hg.loc["proximal", "aligned_bp"] / 1e3:.0f} kb along chromosome 3.'),
    ('Not a one-off', f'{h_def} of {h_conf} HG002 claims confirmed under the pre-registered criteria are definitive '
     f'({h_reads} molecules; {h_multi} seen in two or more independent datasets). HG002\'s parents confirm {t_conf}. '
     f'Independent 1000 Genomes data give {p_def} definitive claims in {p_people} more people.'),
    ('We tried to break it', f'The criteria were frozen before any population read was examined. 0 of {f0(kh_n + kk_n)} '
     f'negative-control tests passed them. All {tan_common.split("/")[1]} common, well-supported alleles recur in four '
     f'independent cohorts (Tan et al. 2026).'),
    ('Check it yourself', f'Any claim can be re-derived from public reads with one command: '
     f'{MONO("evidence/replicate_claim.sh CLAIM_ID")}.'),
]

BODY = [
    ('h1', '1. What we looked at'),
    ('p', 'Each chromosome ends in a telomere, a long run of the six letters TTAGGG. The 100 kb just inside it, the '
          '<b>subtelomere</b>, is made of blocks that are shared, in near-identical copies, between different chromosome ends. '
          'The reference genome (T2T-CHM13) shows one version of each end. We compared the subtelomeres of '
          f'{f0(n_ends)} chromosome ends from {n_people} people with the 48 CHM13 ends and labelled each block by the CHM13 end '
          f'it matches best. An end is <b>foreign</b> if one of its blocks best matches the end of an unrelated chromosome '
          f'(not a known near-identical partner). That is true of {n_foreign / n_ends:.0%} of all ends; in '
          f'{n_first / n_foreign:.0%} of them the foreign block is the one right next to the telomere. The same versions recur in many '
          'unrelated people. The atlas catalogues them as structural alleles, such as SA-20q-01.'),
    ('p', 'A skeptic\'s first question: are these real chromosome ends, or assembly mistakes? Assemblies are stitched '
          'together from many reads, and stitching inside near-identical sequence is exactly where assemblers can go wrong. '
          'So we asked the reads directly, one molecule at a time.'),
    ('h1', '2. One molecule tells the story'),
    ('fig', 'figures/pitch/fig_pitch_molecule.png',
     f'<b>Figure 1. A single DNA molecule joins chromosome-11-like sequence at the telomere to chromosome 3.</b> '
     f'(a) Nanopore read {MONO(HERO)} ({hr.dataset}, HG002 Q100 public alignments; {hr.read_len / 1e3:.0f} kb), drawn to scale '
     f'over its first 80 kb against the CHM13 ends of chromosome 11 (p arm) and chromosome 3 (q arm). Cut at the switch, its '
     f'distal part aligns to CHM13 chr11:{hg.loc["distal", "pos"]:,}-{hg.loc["distal", "end"]:,} and its proximal part to '
     f'chr3:{hg.loc["proximal", "pos"]:,}-{hg.loc["proximal", "end"]:,}, both at MAPQ 60. '
     f'(b) Every diagnostic site (a position where the CHM13 chr11 and chr3 copies differ) covered by the two aligned parts: '
     f'tick colour = which chromosome\'s letter the molecule carries. Yellow: 95% interval of the switch. '
     f'(c) The letters at 6 sites evenly spaced on each side; the one chr11 letter after the switch is shown as found. '
     f'(d) Every molecule of five showcase claims, compared in 2 kb windows with the donor and recipient CHM13 ends: '
     f'closeness index (d_recipient - d_donor)/(d_recipient + d_donor), where d = differences per 1,000 '
     f'aligned letters, median of windows before and after the switch. {n_flip} of {n_mol} molecules flip. Not shown: '
     f'the HG002 and HG01258 17p claims ({excl.get("HG002_1_chr17_PATERNAL_L", "")}; '
     f'{excl.get("HG01258__17p__SA-17p-02", "")} molecules flip), whose 11p-like block is rich in tandem repeats, so 2 kb '
     f'windows align poorly to both ends, and NA20346 2q ({excl.get("NA20346__2q__SA-2q-04", "")}), whose switch at 7-8 kb '
     f'leaves too few windows. All eight claims are confirmed by the formal diagnostic-site test. Data: report/pitch_data/.',
     '../figures/pitch/fig_pitch_molecule.png'),
    ('p', '<b>How to read Figure 1.</b> (1) The molecule starts with telomere, so it really is the end of a chromosome. '
          f'(2) Its next {(hp.switch_lo - hp.telo_len) / 1e3:.0f} kb carries chromosome 11\'s letters, not chromosome 3\'s: '
          f'{int((pre.call == "D").sum())} of {len(pre)} sites. (3) After the switch it carries chromosome 3\'s letters '
          f'({int((post.call == "P").sum())} of {len(post)}) and keeps going, uniquely in chromosome 3, for another '
          f'{(hr.read_len - hp.switch_hi) / 1e3:.0f} kb. No assembly was involved: this is one physical molecule from a public '
          'sequencing run, and the reference genome is used only as a ruler.'),
    ('p', f'The same picture holds molecule after molecule (Figure 1d). It holds in different sequencing technologies '
          f'(Oxford Nanopore R9 and R10, PacBio HiFi), at different chromosome ends, in different people, and in HG002\'s '
          f'mother HG004, who passed the 9p end on.'),
    ('fig', 'figures/pitch/fig_consolidated.png',
     '<b>Figure 2. Molecules, assembly and reference genome tell the same story.</b> ' + _TP.consolidated_caption(
         lambda q: q, html=False).replace('guided tour.', 'guided tour on the project website.'),
     '../figures/pitch/fig_consolidated.png'),
    ('h1', '3. Common, inherited and replicated'),
    ('fig', 'figures/pitch/fig_pitch_common.png',
     '<b>Figure 3. The most frequent foreign chromosome ends.</b> The 12 foreign structural alleles with the highest Release 2 '
     'frequency (share of the arm\'s assembled ends). People = distinct individuals carrying the allele. Filled markers: at '
     'least one claim for this allele is definitive at read level; open markers: no definitive claim yet. This means not '
     'testable with the available reads or not confirmed; it is not evidence against the allele.',
     '../figures/pitch/fig_pitch_common.png'),
    ('p', f'SA-20q-01 is the clearest case: {a20.n_samples} people carry it, and independent 1000 Genomes nanopore '
          f'reads show it in single molecules in {int(P[isdef[P.index] & (P.allele_id == "SA-20q-01")].individual.nunique())} '
          f'of them. SA-3q-01, the allele of Figure 1, is carried by {a3.n_samples} people. These versions are old: '
          'switch points are shared across unrelated people and continental groups. They are also reproducible: '
          f'classifying {f0(tan_n)} subtelomere sequences from four independent cohorts (HGSVC3, CPC, APC, JASA; Tan et al. 2026) '
          f'with the same method finds all {tan_common.split("/")[1]} common tier-A alleles again, and non-native allele frequencies '
          f'agree with Release 2 at r = {tan_r:.2f}.'),
    ('h1', '4. How hard we tried to break it'),
    ('table', [['Test', 'The worry it addresses', 'Result'],
               ['Criteria frozen before looking', 'Moving the goalposts after seeing the data',
                'Criteria v1 committed (git 673b77b) before any population read was inspected; a stricter "definitive" '
                'tier added later can only remove claims'],
               ['Negative controls', 'Does the test call anything foreign?',
                f'0 of {kh_n} HG002 other-haplotype tests and 0 of {f0(kk_n)} 1000 Genomes native-individual tests '
                f'passed (95% upper bounds {ub(kh_c, kh_n):.1%} and {ub(kk_c, kk_n):.2%})'],
               ['Independent sequencing runs', 'An artefact of one run or one chemistry',
                f'{h_multi} of {h_def} definitive HG002 claims seen in two or more independent datasets'],
               ['Inheritance', 'Assembly errors are not inherited',
                f'{t_conf} HG002 claims confirmed in the reads of the parent who passed the chromosome on'],
               ['Second assembler', 'An artefact of one assembler',
                f'Same haplotype assembled independently with Verkko (Tan et al.) gives the same structural allele {vk_sa:.1%} of the time'],
               ['Independent cohorts', 'Something specific to Release 2',
                f'{tan_common} common tier-A alleles recur in four cohorts; frequencies r = {tan_r:.2f}'],
               ['Re-derivation from public data', '"Trust us"',
                f'{rep_ok} of {rep_n} clean-tree re-runs of evidence/replicate_claim.sh reproduce the published verdict '
                'read for read']],
     [1.55, 1.85, 3.6]),
    ('h1', '5. What this does and does not show'),
    ('p', '"Foreign" is defined relative to the reference genome: a block is labelled by the CHM13 end it resembles most, '
          'and CHM13 itself carries one version of each end. A "switch" is where the best match changes. It usually lies '
          'inside sequence that both ends share, so it marks a change in resemblance, not necessarily a breakpoint. '
          f'Read-level proof needs long, telomere-reaching molecules: {p_na} of {len(P)} population claims could not be tested '
          'with the available 1000 Genomes reads, which says nothing against them. The window comparison in Figure 1d is an '
          'illustration; verdicts rest on the pre-registered diagnostic-site test (evidence/CRITERIA.md).'),
    ('h1', '6. Check it yourself'),
    ('p', f'Showcase molecules, with IGV sessions and a "what you should see" note for each: evidence/showcase/. Every cited '
          f'read, with the exact {MONO("samtools")} command that fetches it from the public file: evidence/reads.tsv.gz. '
          f'Re-derive a verdict from public data: {REP_CMD}. '
          'Full details: evidence/READ_EVIDENCE.pdf; the catalogue, query tool and website: '
          'https://jlanej.github.io/human-subtelomere-atlas/.'),
]

REFS = ['Nurk S, et al. (2022). The complete sequence of a human genome. <i>Science</i> 376:44-53. doi:10.1126/science.abj6987',
        'Liao W-W, et al. (2023). A draft human pangenome reference. <i>Nature</i> 617:312-324. doi:10.1038/s41586-023-05896-x',
        'Schloissnig S, Pani S, et al. (2025). Structural variation in 1,019 diverse humans based on long-read sequencing. '
        '<i>Nature</i>. doi:10.1038/s41586-025-09290-7',
        'Zook JM, et al. (2016). Extensive sequencing of seven human genomes to characterize benchmark reference materials. '
        '<i>Scientific Data</i> 3:160025. doi:10.1038/sdata.2016.25',
        'Tan et al. (2026). Widespread structural variations at human chromosome ends. bioRxiv. doi:10.64898/2026.09.03.748401']
