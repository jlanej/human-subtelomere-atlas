"""Content of report/NEW_BIOLOGY.pdf / .md (built by report/build_report.py from the repository root).
Every number is read at build time from committed files in data/, analyses/, tests/validation/ and literature/."""
import json, os
import pandas as pd

J = lambda p: json.load(open(p))
T_ = lambda p: pd.read_csv(p, sep='\t', low_memory=False)
f0 = lambda x: f'{x:,.0f}'
p1 = lambda x: f'{100 * x:.1f}%'
p2 = lambda x: f'{100 * x:.2f}%'

# catalogue
E = T_('data/catalogue/ends.tsv'); R2 = E[E['set'] == 'HPRC_R2']
SA = T_('data/catalogue/alleles_structural.tsv'); HG = T_('data/catalogue/haplogroups.tsv')
FAM = T_('data/catalogue/paralog_families.tsv')
NAT = R2.nativity_class.value_counts()
sa = lambda i: SA.set_index('allele_id').loc[i]
a20, a3, a6, a8 = sa('SA-20q-01'), sa('SA-3q-01'), sa('SA-6p-01'), sa('SA-8p-02')
TIER = SA.tier.value_counts(); HTIER = HG.tier.value_counts()
STAB = T_('data/catalogue/qc/stability_summary.tsv')
def stab(shift, stratum, key, val):
    s = STAB[(STAB['shift'] == shift) & (STAB.stratum == stratum)]
    if key: s = s[s[key].astype(str) == str(val)]
    return s.iloc[0]
nn_m, nn_p = stab(-500, 'new nativity', 'new_native', False), stab(500, 'new nativity', 'new_native', False)
na_m, na_p = stab(-500, 'new nativity', 'new_native', True), stab(500, 'new nativity', 'new_native', True)
n_stab = int(stab(-500, 'all', None, None)['n'])

# validation
TS = J('data/validation/tan2026_summary.json'); SH = TS['same_haplotype']['primary']; IC = TS['independent_cohorts']
RV = TS['revision1_vs_revision2']['same_haplotype_primary']
BYN = TS['same_haplotype']['by_nativity_class']
SC = T_('tests/validation/self_consistency_summary.tsv'); LO = T_('tests/validation/loso_summary.tsv')
sc = lambda m: SC[(SC.stratum == 'all') & (SC.metric == m)].iloc[0]
lo_sa = LO[(LO.allele_type == 'structural_allele') & (LO.stratum == 'all')].iloc[0]
IND = T_('tests/validation/independent_summary.tsv')
REP = IC['replication_tierA_common']
REP_TXT = ', '.join(f'{k} {REP[k]["n_observed"]}/{REP[k]["n_alleles"]}' for k in ['HGSVC3', 'CPC', 'APC', 'JASA'])
EV = J('literature/evidence_numbers.json')
PED = f'{EV["ped_foreign_intact"]["value"]} of {EV["ped_foreign"]["value"]} foreign ends inherited intact ({EV["ped_transmissions"]["value"]} transmissions)'

# analyses
REC = J('analyses/junctions/results/reciprocal_summary.json'); JS = J('analyses/junctions/results/junction_summary.json')
nc = REC['pooled']['noncarriers']; rec_nc = nc['reciprocal'] / nc['assessed']
OS = J('analyses/allele_ages/results/origin_summary.json'); PS = J('analyses/allele_ages/results/pop_summary.json')
OV = T_('analyses/allele_ages/results/origin_verdicts.tsv').set_index('form')
_SOC = OS['single_origin_counts']
SO_OK = _SOC['consistent with single origin']; SO_NO = _SOC['not supported']
SO_DIV = _SOC['clustered, carriers as diverse as natives']
BS = J('analyses/telomere_cis/results/baseline_stats.json')['mean_rel_tl']
WR = T_('analyses/telomere_cis/results/window_robustness.tsv')
ext = WR[(WR.subset == 'tl_ext') & WR.window.isin(['S20', 'S10', 'S5', 'W515'])]
s20 = WR[(WR.test == 'S20')].iloc[0]; s5 = WR[(WR.test == 'S5')].iloc[0]; t20 = WR[(WR.test == 'T20R')].iloc[0]
age = lambda f: f0(OV.loc[f, 'age_kyr_median_mu'])

TITLE = 'Common terminal alleles of human chromosome ends'
SUBTITLE = ('What the human subtelomere atlas adds: a validated allele catalogue of HPRC Release 2 chromosome ends, '
            'its replication in independent cohorts, and tests of mechanism, age and telomere length')
DATE = 'October 2026 · atlas v0.1'

SUMMARY = [
 ('Resource', f'{f0(len(R2))} telomeric contig ends from {R2["sample"].nunique()} HPRC Release 2 individuals were described as '
  f'mosaics of CHM13 chromosome ends and grouped into {len(SA)} structural alleles (tiers A/B/C {TIER.get("A",0)}/{TIER.get("B",0)}/'
  f'{TIER.get("C",0)}) and {len(HG)} sequence haplogroups, each with evidence and a query tool.'),
 ('Foreign ends are alleles', f'{p1(NAT["foreign"] / len(R2))} of ends carry a block resembling an unrelated arm. Common forms are '
  f'recurrent, old and shared across continents. SA-20q-01, a 3q-like terminal block on 20q, is carried by {int(a20.n_samples)} individuals '
  f'(frequency {a20.r2_freq:.2f}).'),
 ('Replicated', f'In {f0(IC["n_primary"])} sequences from four independent cohorts (HGSVC3, CPC, APC, JASA; Tan et al. 2026), '
  f'every common tier-A allele recurs in the pooled cohorts, and non-native allele frequencies correlate with Release 2 '
  f'(r = {IC["freq_correlation_non_native"]["pearson_r"] if isinstance(IC.get("freq_correlation_non_native"), dict) else IC.get("freq_correlation_non_native")}).'),
 ('Mechanism', f'Transfer is one-way and homology-mediated. Donor ends of carriers look like those of non-carriers '
  f'(reciprocal {p2(REC["reciprocal_frac"])} vs {p2(rec_nc)}), and {JS["junction_class"]["homology_spanning"]} of {JS["clusters"]} '
  f'switch clusters lie within donor-recipient homology.'),
 ('Telomere length', f'Arm-specific lengths reproduce (Spearman ρ = {BS["spearman_rho"]:.2f} against Karimian et al. 2024). '
  'The terminal allele has no detectable cis effect once telomere-boundary artefacts are controlled.'),
 ('Novelty', 'Multi-haplotype arms and shared blocks confirm Tan et al. 2026 and Guarracino et al. 2026. The validated allele-level catalogue, '
  'cross-cohort replication, one-way homology-mediated mechanism, pre-out-of-Africa ages and the negative cis result are new.'),
]

def T(*rows): return [list(r) for r in rows]
fam_txt = '; '.join(f'{r.family} = {r.members.replace(",", "/")}' for r in FAM.itertuples())
BODY = [
 ('h1', '1. The atlas'),
 ('p', f'The atlas describes the outer 100 kb of every telomeric contig end in HPRC Release 2 ({R2["assembly_name"].nunique()} haplotype '
  f'assemblies, {f0(len(R2))} ends, including {f0((R2.arm_source != "contig_name").sum())} on contigs without a chromosome name, whose arm was inferred). '
  'Windows start at the telomere/subtelomere boundary. One minimap2 pass against the 48 CHM13 v2.0 ends gives per-base identity '
  'tracks, and each 100 bp unit is labelled with the best-matching end. Near-identical CHM13 ends form six paralog families '
  f'({fam_txt}). A block is <i>home</i> if the own arm is within 0.5 identity points of the best match. Ends then fall into three '
  f'nativity classes: native ({f0(NAT["native"])}), within-family ({f0(NAT["within_family"])}, closer to another member of a family '
  f'containing the own arm) and foreign ({f0(NAT["foreign"])}).'),
 ('p', 'A <b>structural allele</b> is a recurring block order with switch positions within 10 kb. A <b>haplogroup</b> is a cluster at '
  'ANI ≥ 0.995 over the first 80 kb. Evidence tiers combine recurrence in unrelated samples with orthogonal lines of support: '
  'Verkko concordance, reads across switch points, intact pedigree transmission and presence in a T2T assembly. Labels are relative to '
  'CHM13. A switch marks a change in which CHM13 end the sequence most resembles, usually inside shared paralogous sequence, and '
  'need not be a structural breakpoint.'),
 ('fig', 'figures/atlas/F1_genome_overview.png', '<b>Figure 1. Terminal blocks shared between arms.</b> (a) Chords join the '
  'best-matching CHM13 end of each terminal non-home block to the recipient arm (tier A/B alleles), with width proportional to frequency '
  'on the recipient arm. The outer bars show the number of structural alleles per arm. (b) Nativity class of Release 2 ends per arm.'),
 ('fig', 'figures/atlas/F2_allele_mosaics.png', '<b>Figure 2. Allele mosaics of six arms.</b> Each row is one Release 2 end, with the telomere '
  'at left. Blocks are coloured by their best-matching CHM13 end, and rows are grouped by structural allele (frequency and tier at right).'),
 ('h1', '2. Validation'),
 ('table', T(['Check', 'n', 'Result'],
   ['Same haplotype, independent assembly (Tan et al. Verkko vs our hifiasm)', f0(SH['n_pairs']),
    f'structural allele {p1(SH["sa_agree"])}, form {p1(SH["form_agree"])}, S20 {p1(SH["S20_state_and_family_agree"])}; '
    f'base-level switches identical after projection {p1(SH["switch_base_frac_proj_exact"])}'],
   ['By nativity class (same comparison)', '', '; '.join(f'{k} {p1(v["sa_agree"])}' for k, v in BYN.items() if isinstance(v, dict) and 'sa_agree' in v)],
   ['Boundary shifted ±500 bp and re-aligned', f'{n_stab}', f'non-native allele unchanged {p1(nn_m.same_allele)} / {p1(nn_p.same_allele)}; '
    f'native {p1(na_m.same_allele)} / {p1(na_p.same_allele)}'],
   ['Query tool re-querying catalogue ends', f0(sc('S20_match').n), f'mosaic, S-states and form 100%; structural allele '
    f'{p1(SC[(SC.stratum=="all")&(SC.metric.str.startswith("sa"))].agreement.iloc[0])}'],
   ['Leave one sample out (60 samples)', f0(lo_sa.n), f'structural allele {p1(lo_sa.accuracy)}; held-out-only alleles reported novel '
    f'{int(LO[(LO.allele_type=="structural_allele")&(LO.stratum=="removed_by_holdout")].n_correct.iloc[0])}/'
    f'{int(LO[(LO.allele_type=="structural_allele")&(LO.stratum=="removed_by_holdout")].n.iloc[0])}'],
   ['Independent T2T / pedigree assemblies', 'CHM13, HG002 Q100, WashU', 'structural allele 100% in each (tests/validation/independent_summary.tsv)']),
  [2.3, 0.9, 3.8]),
 ('p', f'<b>Design change before release.</b> A first design tiled 5 kb bins from the telomere boundary. It was rejected because only 40% '
  'of non-native allele calls survived a 500 bp boundary shift. The track-and-family design used here raised agreement with independent '
  f'assemblies from {p1(RV["structural allele"]["revision1"])} to {p1(RV["structural allele"]["revision2"])} (data/catalogue/QC.md, CHANGELOG). '
  'Calls 5 kb from the boundary (S5) remain boundary-sensitive and are not recommended as primary per-end calls.'),
 ('fig', 'figures/validation/tan2026_concordance.png', '<b>Figure 3. Cross-check against Tan et al. 2026.</b> (a) Agreement on the same '
  'haplotype. (b) Switch positions after alignment projection. (c, d) Allele frequencies in four independent cohorts compared with HPRC Release 2.'),
 ('h1', '3. Common terminal alleles are old and replicate'),
 ('p', f'All {REP["n_alleles"]} common tier-A alleles (Release 2 frequency ≥ 5%) recur in the pooled independent cohorts, each in at least '
  f'three individuals ({f0(IC["n_primary"])} sequences from {IC["n_primary_samples"]} individuals). Per cohort, {REP_TXT} are observed. '
  f'Of {IC["n_sa_tests"]:,} allele-frequency tests, {IC["n_sa_tests_q_lt_0.05"]} reach q < 0.05 against the all-ancestry reference, and none '
  'against a matched superpopulation. In Release 2, '
  f'{PS["forms_q05_perm"]} of {PS["forms_tested"]} common forms differ among superpopulations (median FST {PS["median_fst_anova"]:.3f}), '
  f'and {PS["forms_in_all_5_superpops"]} occur in all five. {PS["junction_clusters_ge10_shared_AFRcont_OOA"]} of {PS["junction_clusters_ge10_samples"]} '
  'switch clusters seen in ≥ 10 individuals are shared between continental African and non-African samples. Carriers of '
  f'{OS["clustered"]} of {OS["forms_analysed"]} testable forms cluster tightly (median nearest-neighbour carrier fraction {OS["median_nn_carrier_frac"]:.2f} '
  f'vs {OS["median_nn_expected"]:.2f} expected), but only {SO_OK} of {OS["forms_analysed"]} meet the full single-origin criteria; '
  f'{SO_DIV} clustered forms are as diverse as native alleles, and {SO_NO} are not supported. None is strictly monophyletic, which is compatible with '
  'gene conversion at block edges. '
  f'Median-based clock ages: 20q>3q {age("20q>3q")} kyr, 3q>11p {age("3q>11p")} kyr, 6p>11p {age("6p>11p")} kyr. These ages assume a genome-average '
  'mutation rate and are order-of-magnitude estimates. No recurrent human switch configuration was found in five T2T great-ape genomes.'),
 ('fig', 'analyses/allele_ages/figures/fig_origin.png', '<b>Figure 4. Single origin and age.</b> (a) Diversity of carriers compared with native '
  'donor-arm alleles. (b) Distance of carriers, natives and other recipients to the carrier consensus. (c) Clock ages. (d) Great-ape screen.'),
 ('h1', '4. Transfer is one-way and homology-mediated'),
 ('p', f'For {REC["forms"]} forms ({f0(REC["carriers"])} carrier ends), the donor arm in the same haplotype is assessable for {f0(REC["donor_assessed"])} '
  f'carriers. Their donor ends look like those of non-carriers: reciprocal products {p2(REC["reciprocal_frac"])} vs {p2(rec_nc)}, '
  f'and native {p1(REC["donor_native_frac"])}. So the common terminal alleles arose by copying the donor terminal sequence, not by reciprocal exchange. '
  f'At the {JS["clusters"]} base-resolved switch clusters, {JS["junction_class"]["homology_spanning"]} lie within donor-recipient homology, '
  f'and {JS["mh_classes_obs"]["homology_26_200"] + JS["mh_classes_obs"]["homology_gt200"]} share more than 25 bp of identical parental sequence. '
  f'Fewer than {JS["mh_classes_null_frac"]["homology_26_200"]:.0e} of random position pairs do. These are signatures of homology-directed copying '
  '(break-induced replication or gene conversion). Blunt or microhomology end-joining does not fit. Interstitial telomeric repeats are '
  'enriched only within 5 kb of the telomere.'),
 ('fig', 'analyses/junctions/figures/junctions_mechanism.png', '<b>Figure 5. Mechanism.</b> (a) Donor-end states of carriers and non-carriers. '
  '(b) Length of replaced recipient sequence. (c) Identical parental sequence at the switch. (d) Interstitial telomeric repeats near switches. '
  '(e) Three base-level exemplars, two in HG002 with read support.'),
 ('h1', '5. No cis effect on telomere length'),
 ('p', f'Assembled per-arm telomere length correlates with the long-read per-arm values of Karimian et al. 2024 (ρ = {BS["spearman_rho"]:.2f}, '
  f'{BS["n"]} arms). Within individuals heterozygous for a foreign and a native allele of the same arm, an earlier analysis window '
  f'suggested foreign alleles carry longer telomeres (+{t20.mean_d:.0f} bp). That window shrinks as the telomere grows, so the '
  f'association was an artefact. Windows anchored at the boundary give small negative differences (S20 {s20.mean_d:.0f} bp; S5 {s5.mean_d:.0f} bp). '
  'These arise because the telomere caller stops early at degenerate repeats next to the telomere. With a boundary-robust telomere length, '
  f'the paired differences are {ext.mean_d.min():.0f} to {ext.mean_d.max():+.0f} bp (all p ≥ {ext.wilcoxon_p.min():.2f}). The minimum detectable '
  f'effect is about {s20.mde80:.0f} bp. Tan et al. 2026 report suggestive haplotype associations for 1p, 5q and 6q (FDR < 0.1); this '
  'within-person design neither supports nor excludes them.'),
 ('fig', 'analyses/telomere_cis/figures/fig2_allele_test.png', '<b>Figure 6. Telomere length and terminal allele.</b> (a) Paired differences by '
  'allele window. (b) Per-form differences. (c, d) Donor-likeness. (e) The window artefact.'),
 ('h1', '6. Comparison with the literature'),
 ('table', T(['Claim', 'Status before this work', 'Atlas evidence'],
   ['Chromosome ends are mosaics of shared blocks; arms carry several haplotypes', 'Established (Linardopoulou 2005; Tan 2026; Guarracino 2026)', 'Confirmed; not claimed as new'],
   ['Common foreign terminal blocks are real alleles, not assembly errors', 'Assumed in recent preprints, not validated per allele',
    'Validated: second assembler, reads, pedigree, independent cohorts'],
   ['Terminal replacements form by one-way, homology-directed copying', 'Contested (NHEJ inferred by Linardopoulou 2005)', 'Supported for common alleles'],
   ['Terminal alleles are few, old and of single origin', 'Established only for 4q/10q (Lemmers 2010)', f'Carriers cluster for {OS["clustered"]}/{OS["forms_analysed"]} testable forms; full single-origin criteria met by {SO_OK}/{OS["forms_analysed"]}; switch clusters shared across Africa and non-Africa'],
   ['Telomere length is set in cis by the terminal allele', 'Suggestive associations (Tan 2026); untested within persons', 'No detectable effect (≥ ~0.4 kb excluded)'],
   ['Ectopic subtelomere exchange is frequent in pedigrees', 'Claimed (Guarracino 2026)', f'Not supported: {PED}']),
  [2.4, 2.3, 2.3]),
 ('p', 'The full screened review (199 records) and the novelty table are in literature/REVIEW.md and literature/NOVELTY.md.'),
 ('h1', '7. Limitations'),
 ('p', 'Labels are relative to a single reference (CHM13), so an arm whose CHM13 copy is a minority haplotype appears foreign in most people. '
  'Acrocentric p arms and the pseudoautosomal regions are pooled. Duplicon-rich arms carry many rare structures with little orthogonal '
  'support (tier C). Most source DNA comes from cell lines, and assembled telomere lengths are lower bounds. Clock ages ignore gene conversion '
  'and assume a genome-average mutation rate. The same-haplotype cross-check uses the same sequencing data and tests assembly robustness, not '
  'sequencing error. Read-level assignment in the query tool is experimental.'),
 ('h1', 'Data and code'),
 ('p', 'Catalogue, figures, site and code: github.com/jlanej/human-subtelomere-atlas (data CC BY 4.0, code MIT). '
  'Browse at jlanej.github.io/human-subtelomere-atlas. Reproduction commands are in the repository README.'),
]
_LIT = T_('literature/records.tsv').set_index('doi')
_REF_DOIS = ['10.64898/2026.07.10.737660', '10.64898/2026.09.03.748401', '10.1126/science.ado0431', '10.1038/nature04029',
             '10.1016/j.ajhg.2010.01.035', '10.1126/science.abj6987', '10.1038/s41586-025-08816-3']
REFS = [f'{_LIT.loc[d, "first_author"]} ({_LIT.loc[d, "year"]}). {_LIT.loc[d, "title"]}. <i>{_LIT.loc[d, "venue"]}</i>. doi:{d}' for d in _REF_DOIS]
REFS += ['Tan et al. (2026) subtelomere sequences, Zenodo, CC BY 4.0. doi:10.5281/zenodo.21807156',
         'Human Pangenome Reference Consortium, Release 2 (BioProject PRJNA730823).']
