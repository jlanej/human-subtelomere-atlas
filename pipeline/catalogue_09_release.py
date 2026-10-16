"""Step 3: attach validation evidence and flags, assign tiers, compute frequencies and write the v0.1 release tables.
Inputs : work/windows/ends_windows.tsv, work/catalogue/{arms,mosaic_ends,mosaic_blocks,core_blocks,sa_members,sa_alleles,
         haplogroup_members,haplogroups_work}.tsv, data/metadata/samples.tsv (sibling track; read only), and earlier evidence
         tables in $ATLAS_INPUT_DIR/results: subtel_foreign_tips.tsv (prior 20 kb calls, vk_* Verkko concordance, lay_* hifiasm
         layout support), subtel_hg002_ends.tsv (HG002 ONT/HiFi reads crossing switches), subtel_pedigree_transmission.tsv;
         tipdup/tip_placement.tsv (junk_hit = outer junk sequence matching elsewhere in CHM13).
Per-end evidence:
  vk_status   concordant (Verkko contig of the same haplotype matches >= 99.9% and, where tested, places the switch identically) /
              discordant / outer_only (switch beyond the compared window) / not_tested
  lay_status  supported (hifiasm read layout matches and >= 2 reads span the outer 1 kb or the switch lies beyond the one-read
              zone) / one_read_zone (switch inside the zone covered by a single read) / mismatch / match_only / not_tested
  reads_status (HG002 Q100 ends) supported (>= 3 ONT UL reads cross the switch and conflicting reads <= 10% of crossing) / weak / not_tested
  ped_status  (pedigree_aux ends) intact / not_intact transmission (subtel_pedigree_transmission.tsv), child end
Flags: artefact_candidate (4 ends named in the task brief), junk_gt_1kb, junk_chimeric (junk_hit), short_telomere (< 1 kb),
       partial_window, arm_inferred (and arm_low_confidence / arm_conflict), complex_structure (>= 4 core blocks),
       dup_arm_in_haplotype (another end of the same haplotype has the same home arm), vk_discordant, one_read_zone, layout_mismatch
Tiers (structural alleles and haplogroups; evidence lines: R = recurrence in >= 3 unrelated samples; orthogonal O = Verkko
       concordance | read support (HG002 reads or hifiasm layout) | intact pedigree transmission of a member structure/haplogroup |
       presence in a T2T assembly (HG002 Q100; CHM13 counts only for non-native structural alleles and for haplogroups)):
       A = R and >= 1 O, or >= 2 O;  B = R only, or 2 samples and >= 1 O;  C = otherwise (singleton = 1 sample).
Fixed-window states: S5/S10/S20 (5/10/20 kb from the telomere boundary) are the primary per-end states; form = home>foreign arm
       of the S20 call. Concordance is reported against the boundary-anchored states of the telomere_cis track
       (analyses/telomere_cis/results/fixed_window_states.tsv, read only; S20 primary) and, secondarily, the earlier t20 calls.
Frequencies: HPRC Release 2 ends only (CHM13, HG002 Q100 and pedigree ends excluded), denominators = Release 2 ends of the arm
       with an allele call of that type; by superpopulation from samples.tsv.
Outputs: data/catalogue/{ends.tsv, mosaic_blocks.tsv, alleles_structural.tsv, haplogroups.tsv, allele_frequencies_by_superpop.tsv,
         representatives.fa.gz, catalogue.json}
"""
import os, sys, json, gzip, collections, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
VERSION = 'v0.1'
ARTEFACTS = [('HG02647', 'pat', 'CM086565.1', 'R'), ('HG04228', 'mat', 'CM088508.1', 'L'),
             ('HG03486', 'pat', 'CM091971.1', 'R'), ('HG03710', 'pat', 'CM086867.1', 'L')]
os.makedirs(C.CAT, exist_ok=True)
rd = lambda f, **k: pd.read_csv(f, sep='\t', low_memory=False, **k)
W = rd(f'{C.WIN}/ends_windows.tsv'); A = rd(f'{C.CW}/arms.tsv'); ME = rd(f'{C.CW}/mosaic_ends.tsv')
SM = rd(f'{C.CW}/sa_members.tsv'); SA = rd(f'{C.CW}/sa_alleles.tsv').rename(columns={'n_full': 'n_defining'})
HM = rd(f'{C.CW}/haplogroup_members.tsv'); HW = rd(f'{C.CW}/haplogroups_work.tsv').rename(columns={'n_full': 'n_clustered'})
S = rd(f'{C.REPO}/data/metadata/samples.tsv')
F = rd(C.inp('results/subtel_foreign_tips.tsv')); H2 = rd(C.inp('results/subtel_hg002_ends.tsv'))
TR = rd(C.inp('results/subtel_pedigree_transmission.tsv')); PL = rd(C.inp('tipdup/tip_placement.tsv'))

E = W.merge(A.drop(columns=['contig_arm']), on='end_id').merge(ME.drop(columns=['catalogue_arm']), on='end_id', how='left')
E['q'] = np.where(E.set == 'HPRC_R2', E.assembly + '|' + E.contig + '|' + E.side, '')
BS = f'{C.CW}/boundary_stable.tsv'
BSC = ['boundary_stable', 'boundary_stable_S20', 'boundary_stable_allele', 'S20_state_minus500', 'S20_state_plus500', 'allele_minus500', 'allele_plus500']
if os.path.exists(BS):
    E = E.merge(rd(BS)[['end_id'] + BSC], on='end_id', how='left')
else:
    for c in BSC: E[c] = ''
E = E.merge(S[['sample', 'population', 'superpopulation', 'sex']], on='sample', how='left')
E = E.merge(SM[['end_id', 'structure', 'switches', 'n_core_blocks', 'sub_len', 'allele_id', 'assignment', 'candidates']].rename(
    columns={'structure': 'core_structure', 'switches': 'core_switches_sub', 'allele_id': 'structural_allele_id',
             'assignment': 'sa_assignment', 'candidates': 'sa_candidates'}), on='end_id', how='left')
E = E.merge(HM[['end_id', 'haplogroup_id', 'assignment', 'ani_to_medoid', 'nn_ani', 'medoid_aln_ident', 'medoid_aln_frac']].rename(
    columns={'assignment': 'hg_assignment', 'nn_ani': 'hg_nn_ani'}), on='end_id', how='left')
# prior calls and evidence (Release 2)
pc = ['q', 't20_state', 't20_best_end', 'form', 'switch_final', 'switch_res', 'junction', 'junction_cluster', 'junction_cluster_samples',
      'vk_verkko_hit', 'vk_ident', 'vk_switch_covered', 'vk_junction_status', 'lay_raw_match', 'lay_reads_span_1kb', 'lay_one_read_zone',
      'lay_switch_beyond_zone']
E = E.merge(F[pc].rename(columns={c: f'prior_{c}' for c in pc[1:9]}), on='q', how='left')
def vk(r):
    if r.vk_verkko_hit != r.vk_verkko_hit: return 'not_tested'
    if r.vk_verkko_hit is False or r.vk_verkko_hit == 'False' or (r.vk_ident == r.vk_ident and r.vk_ident < 0.999): return 'discordant'
    if r.vk_junction_status == 'outer_differs': return 'discordant'
    if r.vk_junction_status == 'window_limit': return 'outer_only'
    return 'concordant'
def lay(r):
    if r.lay_raw_match != r.lay_raw_match: return 'not_tested'
    if str(r.lay_raw_match) == 'False': return 'mismatch'
    if str(r.lay_switch_beyond_zone) == 'False': return 'one_read_zone'
    if str(r.lay_switch_beyond_zone) == 'True' or (r.lay_reads_span_1kb == r.lay_reads_span_1kb and r.lay_reads_span_1kb >= 2): return 'supported'
    return 'match_only'
E['vk_status'] = [vk(r) for r in E.itertuples()]; E['lay_status'] = [lay(r) for r in E.itertuples()]
# HG002 Q100 read evidence
H2['end_id'] = ['HG002#{}#{}#{}'.format(1 if c.endswith('PATERNAL') else 2, c, 'L' if e == 'p' else 'R') for c, e in zip(H2.contig, H2.end)]
H2['reads_status'] = np.where(H2.ont_ul_dorado_crossing.isna(), 'not_tested',
                              np.where((H2.ont_ul_dorado_crossing >= 3) & (H2.ont_ul_dorado_conflict.fillna(0) <= 0.1 * H2.ont_ul_dorado_crossing), 'supported', 'weak'))
E = E.merge(H2[['end_id', 'reads_status', 'ont_ul_dorado_crossing', 'ont_ul_dorado_conflict']].rename(
    columns={'ont_ul_dorado_crossing': 'hg002_ont_crossing', 'ont_ul_dorado_conflict': 'hg002_ont_conflict'}), on='end_id', how='left')
E['reads_status'] = E.reads_status.fillna('not_tested')
# pedigree transmissions (auxiliary ends)
def ped_id(x):
    h, arm = x.split('|'); s, hh = h.split('_'); return f'{s}#{hh}#chr{arm[:-1]}#{"L" if arm.endswith("p") else "R"}'
TR['child_id'] = TR.child_end.map(ped_id); TR['parent_id'] = TR.best_parent_end.map(ped_id)
sa_of = dict(zip(E.end_id, E.structural_allele_id)); hg_of = dict(zip(E.end_id, E.haplogroup_id))
TR['child_sa'] = TR.child_id.map(sa_of); TR['parent_sa'] = TR.parent_id.map(sa_of)
TR['child_hg'] = TR.child_id.map(hg_of); TR['parent_hg'] = TR.parent_id.map(hg_of)
pst = {r.child_id: ('intact' if r.intact else 'not_intact') for r in TR.itertuples()}
E['ped_status'] = E.end_id.map(pst).fillna(pd.Series(np.where(E.set == 'pedigree_aux', 'parent_or_untested', 'not_applicable'), index=E.index))
# flags
junk_hit = set((PL.asm + '|' + PL.contig + '|' + PL.end)[PL.junk_hit.notna()])
art = {f'{s}#{1 if h in ("pat", "hap1") else 2}#{c}#{e}' for s, h, c, e in ARTEFACTS}
art_found = set(E.end_id) & art
# duplicate arm in a haplotype: same specific arm twice (pooled-group ends with an inferred arm are not compared)
armkey = np.where((E.arm_source == 'contig_name') | (E.catalogue_arm == E.home_arm), E.home_arm.astype(str), E.end_id)
hapkey = E['sample'].astype(str) + '#' + E.set + '#' + E.haplotype.astype(str) + '#' + pd.Series(armkey, index=E.index)
dup = hapkey.map(hapkey.value_counts()) > 1
def flags(r, d):
    f = []
    if r.end_id in art: f.append('artefact_candidate')
    if r.junk_len > 1000: f.append('junk_gt_1kb')
    if r.q in junk_hit: f.append('junk_chimeric')
    if r.telo_len < 1000: f.append('short_telomere')
    if r.telomere_to_window_end: f.append('telomere_to_window_end')
    if r.window_extended: f.append('window_extended')
    if r.window_len < C.WINDOW or (r.sub_len == r.sub_len and r.sub_len < 80000): f.append('partial_window')
    if r.arm_source != 'contig_name': f.append('arm_inferred')
    if r.arm_confidence == 'low': f.append('arm_low_confidence')
    if r.arm_confidence == 'conflict': f.append('arm_conflict')
    if r.n_core_blocks == r.n_core_blocks and r.n_core_blocks >= 4: f.append('complex_structure')
    if d: f.append('dup_arm_in_haplotype')
    if r.vk_status == 'discordant': f.append('vk_discordant')
    if r.lay_status == 'one_read_zone': f.append('one_read_zone')
    if r.lay_status == 'mismatch': f.append('layout_mismatch')
    return ','.join(f)
E['flags'] = [flags(r, d) for r, d in zip(E.itertuples(), dup)]
# t20 concordance with the earlier calls
def t20cmp(r):
    if r.set != 'HPRC_R2' or r.arm_source != 'contig_name' or r.prior_t20_state != r.prior_t20_state: return ''
    if r.t20_state != r.prior_t20_state: return 'state_differs'
    if r.t20_state == 'foreign' and r.t20_best_end != r.prior_t20_best_end: return 'foreign_end_differs'
    return 'concordant'
E['t20_vs_prior'] = [t20cmp(r) for r in E.itertuples()]
E['form'] = np.where(E.S20_state == 'foreign', E.home_arm.astype(str) + '>' + E.S20_family.astype(str), '')   # over paralog families
from subtel_atlas import tracks as TK
FAMS = TK.load_families(f'{C.CAT}/paralog_families.tsv')
E['nativity_class'] = [TK.nativity_class(s_, a_, FAMS) for s_, a_ in zip(E.core_structure, E.catalogue_arm)]   # native / within_family / foreign
E['prior_switch_sub'] = E.prior_switch_final - E.s0
# boundary-anchored states of the telomere_cis track (if present)
CIS = f'{C.REPO}/analyses/telomere_cis/results/fixed_window_states.tsv'
if os.path.exists(CIS):
    X = rd(CIS)[['q', 'S5_state', 'S10_state', 'S20_state', 'S20_best_end']].rename(columns=lambda c: c if c == 'q' else f'cis_{c}')
    E = E.merge(X, on='q', how='left')
    def scmp(r):
        if not isinstance(r.cis_S20_state, str) or r.catalogue_arm != r.home_arm or r.arm_source != 'contig_name': return ''
        if r.S20_state != r.cis_S20_state: return 'state_differs'
        if r.S20_state == 'foreign' and str(r.S20_best_end)[3:] != str(r.cis_S20_best_end): return 'foreign_end_differs'
        return 'concordant'
    E['S20_vs_cis'] = [scmp(r) for r in E.itertuples()]
else:
    for c in ('cis_S5_state', 'cis_S10_state', 'cis_S20_state', 'cis_S20_best_end', 'S20_vs_cis'): E[c] = ''

# ---- allele-level evidence and tiers
def aggregate(T, key, members, native_col=None):
    out = []
    for aid, g in members.groupby(key):
        core = g[g.set != 'pedigree_aux']; r2 = core[core.set == 'HPRC_R2']
        samples = set(core['sample'])          # HG002 Q100 and HG002 R2 are the same sample
        ev = []
        if (core.vk_status == 'concordant').any(): ev.append('verkko')
        if (core.lay_status == 'supported').any() or (core.reads_status == 'supported').any(): ev.append('reads')
        if aid in T: ev.append('pedigree')
        native = native_col.get(aid, False) if native_col else False
        if (core.set == 'HG002_Q100').any() or ((core.set == 'CHM13').any() and not native): ev.append('t2t_assembly')
        rec = len(samples) >= 3; nO = len(ev)
        tier = 'A' if (rec and nO >= 1) or nO >= 2 else ('B' if rec or (len(samples) == 2 and nO >= 1) else 'C')
        fl = core['flags'].fillna('')
        out.append({key: aid, 'n_ends': len(core), 'n_full': int((~core['flags'].fillna('').str.contains('partial_window')).sum()),
                    'n_r2_ends': len(r2), 'n_samples': len(samples), 'n_chm13': int((core.set == 'CHM13').sum()),
                    'n_hg002_q100': int((core.set == 'HG002_Q100').sum()), 'n_pedigree_aux': int((g.set == 'pedigree_aux').sum()),
                    'n_vk_concordant': int((core.vk_status == 'concordant').sum()), 'n_vk_discordant': int((core.vk_status == 'discordant').sum()),
                    'n_layout_supported': int((core.lay_status == 'supported').sum()), 'n_hg002_reads_supported': int((core.reads_status == 'supported').sum()),
                    'n_pedigree_intact': T.get(aid, 0), 'evidence': ','.join(['recurrent'] * rec + ev), 'tier': tier,
                    'singleton': len(samples) == 1,
                    'n_artefact_candidates': int(fl.str.contains('artefact_candidate').sum()),
                    'n_flagged_junk': int(fl.str.contains('junk_').sum()), 'n_arm_inferred': int(fl.str.contains('arm_inferred').sum())})
    return pd.DataFrame(out)
ped_sa = collections.Counter(r.child_sa for r in TR.itertuples() if r.intact and isinstance(r.child_sa, str) and r.child_sa == r.parent_sa)
ped_hg = collections.Counter(r.child_hg for r in TR.itertuples() if r.intact and isinstance(r.child_hg, str) and r.child_hg == r.parent_hg)
mem_sa = E[E.structural_allele_id.notna() & (E.structural_allele_id != '')]
mem_hg = E[E.haplogroup_id.notna() & (E.haplogroup_id != '')]
SAx = SA.merge(aggregate(ped_sa, 'structural_allele_id', mem_sa, dict(zip(SA.allele_id, SA.is_native))).rename(
    columns={'structural_allele_id': 'allele_id'}), on='allele_id', how='left')
HGx = HW.merge(aggregate(ped_hg, 'haplogroup_id', mem_hg), on='haplogroup_id', how='left')
# arm totals and overall frequency (Release 2)
r2 = E[E.set == 'HPRC_R2']
den_sa = r2[r2.structural_allele_id.fillna('') != ''].groupby('catalogue_arm').size(); den_hg = r2[r2.haplogroup_id.fillna('') != ''].groupby('catalogue_arm').size()
SAx['r2_freq'] = (SAx.n_r2_ends / SAx.arm.map(den_sa)).round(4); HGx['r2_freq'] = (HGx.n_r2_ends / HGx.arm.map(den_hg)).round(4)
# structural allele <-> haplogroup cross-tab
ct = mem_sa[mem_sa.set != 'pedigree_aux'].groupby(['structural_allele_id', 'haplogroup_id']).size()
def top(series_idx, g):
    s = g.sort_values(ascending=False); return ';'.join(f'{k}:{v}' for k, v in s.head(5).items())
SAx['haplogroups_top'] = SAx.allele_id.map({a: top(None, g.droplevel(0)) for a, g in ct.groupby(level=0)})
ct2 = mem_hg[mem_hg.set != 'pedigree_aux'].groupby(['haplogroup_id', 'structural_allele_id']).size()
HGx['structural_alleles_top'] = HGx.haplogroup_id.map({a: top(None, g.droplevel(0)) for a, g in ct2.groupby(level=0)})
HGx['structural_purity'] = HGx.haplogroup_id.map({a: round(g.max() / g.sum(), 4) for a, g in ct2.groupby(level=0)})
# prior junction clusters represented in each structural allele
jc = mem_sa[mem_sa.prior_junction_cluster.notna()].groupby('structural_allele_id').prior_junction_cluster.agg(lambda x: ';'.join(f'{k}:{v}' for k, v in x.value_counts().head(3).items()))
SAx['prior_junction_clusters'] = SAx.allele_id.map(jc)
SAx['form'] = [f"{a}>{s.split('|')[0]}" if s.split('|')[0] not in ('home', 'none') else '' for a, s in zip(SAx.arm, SAx.structure)]
SAx['nativity_class'] = [TK.nativity_class(s_, a_, FAMS) for s_, a_ in zip(SAx.structure, SAx.arm)]
REP_MIN_SAMPLES = 3
SAx['representative_in_fasta'] = SAx.n_samples >= REP_MIN_SAMPLES; HGx['representative_in_fasta'] = HGx.n_samples >= REP_MIN_SAMPLES
# ---- frequencies by superpopulation
fr = []
for kind, col, den in (('structural_allele', 'structural_allele_id', den_sa), ('haplogroup', 'haplogroup_id', den_hg)):
    x = r2[r2[col].fillna('') != ''].copy(); x['superpopulation'] = x.superpopulation.fillna('NA')
    for sp, xs in list(x.groupby('superpopulation')) + [('ALL', x)]:
        d = xs.groupby('catalogue_arm').size(); c = xs.groupby(['catalogue_arm', col]).agg(n=('end_id', 'size'), ns=('sample', 'nunique'))
        for (arm, aid), v in c.iterrows():
            fr.append(dict(allele_type=kind, allele_id=aid, arm=arm, superpopulation=sp, n_haplotype_ends=int(v.n), n_samples=int(v.ns),
                           n_arm_ends_assessed=int(d[arm]), frequency=round(v.n / d[arm], 4)))
FR = pd.DataFrame(fr).sort_values(['allele_type', 'arm', 'allele_id', 'superpopulation'])
# ---- tables
endcols = ['end_id', 'set', 'sample', 'haplotype', 'pansn_hap', 'population', 'superpopulation', 'sex', 'assembly', 'assembly_name',
           'assembly_accession', 'assembly_url', 'assembler', 'hifi_input', 'contig', 'contig_len', 'side', 'win_start', 'win_end', 'strand', 'window_len',
           'window_sha256', 'window_extended', 'junk_len', 'telo_len', 'telo_source', 'telomere_to_window_end', 'tiptable_junk_len', 'tiptable_telo_len', 'home_arm', 'catalogue_arm', 'arm_source', 'arm_confidence', 'probe_arm', 'mosaic_arm',
           'profile', 'structure', 'switches_sub', 'switch_res', 'core_structure', 'core_switches_sub', 'n_core_blocks', 'nativity_class', 'sub_len',
           'S5_state', 'S5_family', 'S5_best_end', 'S10_state', 'S10_family', 'S10_best_end', 'S20_state', 'S20_family', 'S20_best_end',
           'S20_best_ident', 'form', 'boundary_stable', 'boundary_stable_S20', 'boundary_stable_allele', 'S20_state_minus500', 'S20_state_plus500',
           'allele_minus500', 'allele_plus500',
           't20_state', 't20_best_end', 't20_best_ident', 'structural_allele_id', 'sa_assignment', 'sa_candidates',
           'haplogroup_id', 'hg_assignment', 'ani_to_medoid', 'hg_nn_ani', 'medoid_aln_ident', 'medoid_aln_frac',
           'vk_status', 'lay_status', 'reads_status', 'hg002_ont_crossing', 'hg002_ont_conflict', 'ped_status', 'flags',
           'cis_S5_state', 'cis_S10_state', 'cis_S20_state', 'cis_S20_best_end', 'S20_vs_cis',
           'prior_t20_state', 'prior_t20_best_end', 'prior_form', 'prior_switch_final', 'prior_switch_sub', 'prior_junction_cluster',
           'prior_junction_cluster_samples', 't20_vs_prior', 'q']
E = E.sort_values(['set', 'end_id'])
C.write_tsv_atomic(E[endcols].rename(columns={'q': 'tip_key_r2'}), f'{C.CAT}/ends.tsv')
MB = rd(f'{C.CW}/mosaic_blocks.tsv'); CB = rd(f'{C.CW}/core_blocks.tsv')
MB = MB.rename(columns={'n_grid': 'n_grid_points'}).assign(level='raw', source_arm=MB.source.map(lambda s: C.arm_of(s) if s != 'none' else 'none'))
MB['family'] = np.where(MB.label.isin(['home', 'none']), '', MB.label)
CB = CB.rename(columns={'core_block': 'block'}).assign(level='core', source_arm=CB.source.map(lambda s: C.arm_of(s) if s != 'none' else 'none'))
bc = ['end_id', 'level', 'block', 'label', 'family', 'source', 'source_arm', 'best_end', 'is_home', 'sub_start', 'sub_end', 'tip_start', 'tip_end',
      'n_grid_points', 'ident', 'mean_tied', 'chm13_dist_start', 'chm13_dist_end', 'next_switch_res']
CB['label'] = np.where(CB.is_home, 'home', CB.source); CB['family'] = np.where(CB.label.isin(['home', 'none']), '', CB.label)
for c in ['label', 'family', 'best_end', 'mean_tied', 'n_grid_points']:
    if c not in CB: CB[c] = np.nan
MBo = pd.concat([MB, CB], ignore_index=True)[bc].sort_values(['end_id', 'level', 'block'], ascending=[True, False, True])
for c in ['chm13_dist_start', 'chm13_dist_end', 'tip_start', 'tip_end', 'n_grid_points']:
    MBo[c] = MBo[c].astype('Int64')
C.write_tsv_atomic(MBo, f'{C.CAT}/mosaic_blocks.tsv')
sac = ['allele_id', 'arm', 'structure', 'form', 'is_native', 'nativity_class', 'n_core_blocks', 'switches_median', 'switches_min', 'switches_max', 'tier', 'evidence',
       'singleton', 'n_ends', 'n_defining', 'n_full', 'n_r2_ends', 'r2_freq', 'n_samples', 'n_chm13', 'n_hg002_q100', 'n_pedigree_aux', 'n_pedigree_intact',
       'n_vk_concordant', 'n_vk_discordant', 'n_layout_supported', 'n_hg002_reads_supported', 'n_artefact_candidates', 'n_flagged_junk',
       'n_arm_inferred', 'representative_end_id', 'representative_in_fasta', 'haplogroups_top', 'prior_junction_clusters']
C.write_tsv_atomic(SAx[sac].sort_values(['arm', 'allele_id']), f'{C.CAT}/alleles_structural.tsv')
hgc = ['haplogroup_id', 'arm', 'tier', 'evidence', 'singleton', 'n_ends', 'n_clustered', 'n_full', 'n_partial', 'n_r2_ends', 'r2_freq', 'n_samples', 'n_chm13',
       'n_hg002_q100', 'n_pedigree_aux', 'n_pedigree_intact', 'n_vk_concordant', 'n_layout_supported', 'n_hg002_reads_supported',
       'medoid_end_id', 'representative_in_fasta', 'min_within_ani', 'mean_within_ani', 'medoid_nearest_other_member_ani', 'min_medoid_aln_ident', 'median_medoid_aln_ident',
       'min_medoid_aln_frac', 'structural_alleles_top', 'structural_purity', 'n_artefact_candidates', 'n_arm_inferred', 'threshold_ani', 'compared_len']
C.write_tsv_atomic(HGx[hgc].sort_values(['arm', 'haplogroup_id']), f'{C.CAT}/haplogroups.tsv')
C.write_tsv_atomic(FR, f'{C.CAT}/allele_frequencies_by_superpop.tsv')
# ---- representative sequences: medoid windows (telomere-anchored, tip frame; the telomere prefix s0 is kept and reported) of every
# structural allele and haplogroup observed in >= REP_MIN_SAMPLES samples, one record per distinct medoid end (header lists all IDs it
# represents). Other windows are retrievable with pipeline/fetch_end_window.py (coordinates + sha256 in ends.tsv).
byend = collections.defaultdict(list)
for r in SAx[SAx.representative_in_fasta].itertuples(): byend[r.representative_end_id].append(r.allele_id)
for r in HGx[HGx.representative_in_fasta].itertuples(): byend[r.medoid_end_id].append(r.haplogroup_id)
seqs = {n: s for n, s in C.readfa(f'{C.WIN}/ends.fa') if n in byend}
s0 = dict(zip(E.end_id, E.s0))
order = sorted(byend, key=lambda e: sorted(byend[e])[0])
with gzip.GzipFile(f'{C.CAT}/representatives.fa.gz.tmp', 'wb', mtime=0, compresslevel=9) as g:
    for e in order:
        ids = sorted(byend[e], key=lambda i: (not i.startswith('SA'), i)); sq = seqs[e]
        g.write(f'>{ids[0]} ids={",".join(ids)} end_id={e} len={len(sq)} subtelomere_start={int(s0[e])}\n'.encode())
        for i in range(0, len(sq), 100): g.write((sq[i:i + 100] + '\n').encode())
os.replace(f'{C.CAT}/representatives.fa.gz.tmp', f'{C.CAT}/representatives.fa.gz')
# ---- compact JSON for the website
cbk = {e: g for e, g in CB.groupby('end_id')}
frj = collections.defaultdict(dict)
for r in FR.itertuples(): frj[r.allele_id][r.superpopulation] = r.frequency
J = dict(version=VERSION, window_bp=C.WINDOW, haplogroup_threshold_ani=float(HW.threshold_ani.iloc[0]), arms={})
for arm in sorted(set(SAx.arm) | set(HGx.arm)):
    sa = SAx[SAx.arm == arm].sort_values('allele_id'); hg = HGx[HGx.arm == arm].sort_values('haplogroup_id')
    J['arms'][arm] = dict(
        n_r2_ends=int((r2.catalogue_arm == arm).sum()),
        structural_alleles=[dict(id=r.allele_id, structure=r.structure, native=bool(r.is_native), nativity_class=r.nativity_class, tier=r.tier, evidence=r.evidence,
                                 n_ends=int(r.n_ends), n_samples=int(r.n_samples), freq=frj[r.allele_id], switches=r.switches_median if isinstance(r.switches_median, str) else '',
                                 representative=r.representative_end_id,
                                 blocks=[dict(source=b.source_arm if b.source_arm != arm else 'home', start=int(b.sub_start), end=int(b.sub_end))
                                         for b in cbk[r.representative_end_id].sort_values('block').itertuples()],
                                 members=sorted(mem_sa[(mem_sa.structural_allele_id == r.allele_id) & (mem_sa.set != 'pedigree_aux')].end_id))
                            for r in sa.itertuples()],
        haplogroups=[dict(id=r.haplogroup_id, tier=r.tier, n_ends=int(r.n_ends), n_samples=int(r.n_samples), freq=frj[r.haplogroup_id],
                          medoid=r.medoid_end_id, structural_alleles=r.structural_alleles_top,
                          members=sorted(mem_hg[(mem_hg.haplogroup_id == r.haplogroup_id) & (mem_hg.set != 'pedigree_aux')].end_id))
                     for r in hg.itertuples()])
C.write_tsv_atomic(SAx[sac].sort_values(['arm', 'allele_id']), f'{C.CAT}/alleles_structural.tsv')
C.write_tsv_atomic(HGx[hgc].sort_values(['arm', 'haplogroup_id']), f'{C.CAT}/haplogroups.tsv')
if os.path.exists(f'{C.CW}/id_history.tsv'):
    import shutil; shutil.copyfile(f'{C.CW}/id_history.tsv', f'{C.CAT}/id_history.tsv.tmp'); os.replace(f'{C.CAT}/id_history.tsv.tmp', f'{C.CAT}/id_history.tsv')
tmp = f'{C.CAT}/catalogue.json.tmp'; json.dump(J, open(tmp, 'w'), separators=(',', ':'), sort_keys=True); os.replace(tmp, f'{C.CAT}/catalogue.json')
TR.to_csv(f'{C.CW}/pedigree_transmission_alleles.tsv', sep='\t', index=False)
print('ends', len(E), 'SA', len(SAx), SAx.tier.value_counts().to_dict(), 'HG', len(HGx), HGx.tier.value_counts().to_dict(),
      'artefact candidates found', len(art_found))
