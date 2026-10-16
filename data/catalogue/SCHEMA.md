# Allele catalogue v0.1 - schema

Coordinates: **tip frame** = position 0 is the outermost base of the contig end, reading inward (telomere -> centromere); R ends
are reverse-complemented. **Subtelomere coordinate** (`sub_*`, `switches_sub`) = tip-frame position minus `s0`
(`s0 = junk_len + telo_len`, the masked telomere plus any outer non-telomeric bases). Contig window coordinates (`win_start`,
`win_end`) are 1-based inclusive on the deposited contig. CHM13 distances (`chm13_dist_*`) are bases from the telomeric end of
the named CHM13 v2.0 chromosome end. Arms are written `1p`, `20q`; CHM13 ends `chr20q`. Missing values are empty cells.

Blocks, switches and structural calls are reported in coordinates from the telomere/subtelomere boundary (`s0`), so they do not depend
on telomere length. A **switch** (junction) is the position where the sequence changes which CHM13 end it most resembles; inside
shared, paralogous subtelomeric sequence this is a change in paralog similarity, not necessarily a structural breakpoint.
'Foreign' is a distance-dependent property: see `profile` (1 kb steps from the boundary) and the S5/S10/S20 states.

Mosaic method (revision 2, atlas/subtel_atlas/tracks.py): each subtelomere is aligned once to the 48 CHM13 ends; per-base match tracks per CHM13 end give a sliding 5 kb identity on a 100 bp grid anchored at the contig tip, so labels do not depend on where the boundary is called. A grid point is `home` if a home end is within 0.5 percentage points of the best identity, otherwise it carries the paralog family of the best end (paralog_families.tsv). Remaining boundary sensitivity is reported in QC.md section 9 (structural-allele stability of re-aligned windows at boundary shifts of -1 kb to +1 kb, by stratum). Before the redesign, the external same-haplotype check found allele agreement of 98.8% with identical boundaries and 73.4% for > 1 kb differences (data/validation/README.md).

In plain language (details and counts: QC.md section 10):

* **Home rule.** At each 100 bp point the label is `home` when the CHM13 copy of the end's own arm (or of any arm of its pooled group)
  is within 0.5 percentage points of the best-matching CHM13 end; home wins ties.
* **Paralog families.** Some CHM13 ends share long, near-identical duplicated blocks, so the single best end can flip between them
  with small sequence or boundary changes. Two CHM13 ends form a family when, across the catalogue, they are tied with each other at
  >= 75% of the points where the less often tied of the two is tied (`paralog_families.tsv`). A non-home point is labelled with the
  family of its best end (`PFnn`), or with that end's arm if it is in no family. A family is a property of the CHM13 ends over the
  whole catalogue, not a statement that its members are tied inside every block; the single best end is kept in `best_end`
  (mosaic_blocks.tsv) and `S*_best_end`, and `mean_tied` gives the local number of tied ends.
* **none.** No CHM13 end reaches 80% identity over the surrounding 5 kb: sequence without a close CHM13 counterpart. It is neither
  home nor foreign evidence. In a core structure a `none` block is >= 10 kb.
* **Nativity class** (`nativity_class`): `native` = only home (and none) blocks; `within_family` = every non-home block carries a family
  that contains the own arm (sequence closer to another member's CHM13 copy of a shared family block than to the own-arm copy);
  `foreign` = at least one block labelled with a family or arm that does not contain the own arm.
* **Worked example.** End HG02129#2#CM085868.1#R (20q, allele SA-20q-01): its outer 31,040 bp best match CHM13 3q. The previous
  build called the form `20q>3q`; it is now `20q>PF04` (core structure `PF04|home`, class foreign) because CHM13 3q and 19p form
  family PF04. Within this block the end is specifically 3q-like (median identity 99.92% to 3q, 98.34% to 19p, 98.36% to its own
  arm 20q), and `S20_best_end` keeps `chr3q`. A 1p end with core structure `home|PF02|none|home` (PF02 = 1p, 5q, 6q) is
  `within_family`: its middle block is closer to the 5q or 6q copy of the shared block than to CHM13 1p.

Catalogue arm groups (`catalogue_arm`): ends whose distal sequence is shared between chromosomes are pooled -
`acroP` = 13p, 14p, 15p, 21p, 22p; `PAR1` = Xp, Yp; `PAR2` = Xq, Yq; every other arm is its own group. All CHM13 ends of a
group count as home sequence. Allele and haplogroup IDs are specific to this release (v0.1); later releases will publish an ID
mapping.

Sets: `HPRC_R2` (HPRC Release 2 assemblies), `CHM13` (T2T CHM13 v2.0), `HG002_Q100` (HG002 Q100 v1.1),
`pedigree_aux` (WashU PAN010/PAN011/PAN027/PAN028 v1.1 assemblies; used only as validation evidence, never counted in alleles'
member counts or frequencies).

## ends.tsv (one row per chromosome end)

| column | definition |
|---|---|
| end_id | stable end identifier `<sample>#<hap>#<contig>#<L/R>` (Release 2: PanSN contig name + side, e.g. `HG00097#1#CM094060.1#L`; CHM13 `CHM13#0#chr1#L`; HG002 Q100 `HG002#1#chr1_PATERNAL#L` (1 = paternal, 2 = maternal); pedigree `PAN027#pat#chr1#L`) |
| set | HPRC_R2, CHM13, HG002_Q100 or pedigree_aux |
| sample | sample ID |
| haplotype | haplotype label of the assembly (hap1/hap2/mat/pat; CHM13 haploid; HG002 pat/mat; pedigree h1/h2/mat/pat). This is the assembly's own label: it denotes parent of origin only for trio-phased assemblies (mat/pat); hap1/hap2 labels are assembly-specific and are not comparable between assemblers or releases (hifiasm and Verkko labels of the same sample agree for only about half of the ends outside trio-phased samples) |
| pansn_hap | PanSN haplotype number (1/2; 0 = not applicable) |
| population, superpopulation, sex | from data/metadata/samples.tsv (1000 Genomes labels); empty for CHM13 and pedigree ends |
| assembly | assembly label used in the analysis (e.g. `HG00097.hap1.hifiasm`) |
| assembly_name | deposited assembly name (e.g. `HG00097_hap1_hprc_r2_v1.0.1`) or source file name |
| assembly_accession | GenBank assembly accession (Release 2: from the HPRC assembly index; CHM13 GCA_009914755.4) |
| assembly_url | public source of the assembly FASTA (s3:// paths are served at https://human-pangenomics.s3.amazonaws.com/) |
| assembler | hifiasm / verkko (Release 2), T2T, WashU_v1.1 |
| hifi_input | HiFi input type recorded in the Release 2 tip table |
| contig, contig_len | contig name and length (bp) |
| side | L = contig start, R = contig end |
| win_start, win_end | window on the contig, 1-based inclusive |
| strand | + (L ends, as deposited) or - (R ends, reverse-complemented into tip frame) |
| window_len | window length (bp); 100,000 unless the contig is shorter, or extended (see window_extended) |
| window_extended | the window was extended beyond 100 kb to ceil((s0 + 80 kb) / 5 kb) x 5 kb because the telomere boundary lies beyond 20 kb (keeps >= 80 kb of subtelomere) |
| window_sha256 | SHA-256 of the upper-case tip-frame window sequence (no header/newlines) |
| junk_len | bases outside the telomere at the very tip (before the first telomeric repeat block) |
| telo_len | length of the terminal (CCCTAA)n block in tip frame |
| telo_source | window_call: junk_len/telo_len called on this window with subtel_atlas.telomere.find_block (all ends) |
| telomere_to_window_end | the telomere block reaches the end of the window (no subtelomere observed; none in v0.1) |
| tiptable_junk_len, tiptable_telo_len | values in the Release 2 tip table (called on 20 kb windows; capped near 20 kb), kept for traceability |
| home_arm | the end's own arm: from the contig name, or inferred (see arm_source) |
| catalogue_arm | arm or pooled group used for alleles/haplogroups |
| arm_source | contig_name, probe (10 kb probes at 250 kb / 1 Mb mapped to CHM13), mosaic (best CHM13 end of the proximal window), none |
| arm_confidence | named; high (two probes agree); medium (one informative probe); conflict (probes disagree at group level); mosaic_group (mosaic call inside a pooled group); low (other mosaic calls); none |
| probe_arm, mosaic_arm | arm proposed by the probe rule / the mosaic rule (recorded for all ends with data) |
| profile | track label every 1 kb of subtelomere (every 10th grid point): O = home, F = foreign (a paralog family or single CHM13 end), . = none |
| structure | raw block labels telomere -> centromere (`home`, a paralog family `PFnn`, the arm of a CHM13 end in no family, `none`), `|`-separated |
| switches_sub | raw block boundaries (subtelomere coordinates), comma-separated |
| switch_res | per raw boundary: base (base-level change point on the label tracks) or grid (100 bp grid boundary) |
| core_structure | structure after dropping blocks < 10 kb (defines structural alleles); tokens as in structure |
| core_switches_sub | core block boundaries (subtelomere coordinates) |
| n_core_blocks | number of core blocks |
| nativity_class | native / within_family / foreign, from core_structure and the paralog families (see above; QC.md section 10); empty if no core structure |
| sub_len | subtelomere length available in the window (window_len - s0) |
| S5_state, S10_state, S20_state | fixed-window states of the first 5 / 10 / 20 kb of subtelomere = majority track label of the window's grid points: own (home; a home end is the best end at most of those points), tie (home, but another end is best at most points), foreign, none |
| S5_family, S10_family, S20_family | paralog family (or single-end arm) of a foreign window |
| S5_best_end, S10_best_end, S20_best_end, S20_best_ident | most frequent best single CHM13 end over the window's grid points (auxiliary, for continuity) and the mean best identity |
| boundary_stable | boundary_stable_S20 and boundary_stable_allele |
| boundary_stable_S20, boundary_stable_allele | the S20 call (state and family) / the structural allele is unchanged when the boundary is moved by -500 and +500 bp (pipeline/catalogue_08e_boundary_stable.py: identity tracks shifted without re-alignment; agreement with re-aligned windows is reported in QC.md) |
| allele_minus500, allele_plus500 | structural allele with the boundary moved by -500 / +500 bp |
| S20_state_minus500, S20_state_plus500 | S20 state with the boundary moved by -500 / +500 bp |
| form | `<home_arm>><S20_family>` for S20-foreign ends (families from paralog_families.tsv; single ends by arm) |
| t20_state, t20_best_end, t20_best_ident | the earlier window tip[s0:20 kb] (shrinks as the telomere grows), same track rule; kept only for comparison with earlier calls |
| structural_allele_id | SA-<catalogue_arm>-NN, empty if unassigned |
| sa_assignment | full (defines the allele), partial_unique / partial_ambiguous (short window, prefix-compatible), arm_lowconf, aux (pedigree), unassigned_*. With several compatible alleles the one with the nearest median switch positions is chosen, then the one with most defining ends |
| sa_candidates | all compatible alleles for partial/aux assignments |
| haplogroup_id | HG-<catalogue_arm>-NN, empty if unassigned |
| hg_assignment | full (clustered), partial, arm_lowconf, aux, unassigned_* |
| ani_to_medoid | k-mer ANI to the haplogroup medoid (clustered members) |
| hg_nn_ani | ANI to the nearest clustered end of the arm |
| medoid_aln_ident, medoid_aln_frac | alignment (mappy asm20) of the end's compared subtelomere to its haplogroup medoid: gap-compressed identity and aligned fraction of the end |
| vk_status | Verkko concordance: concordant / discordant / outer_only / not_tested |
| lay_status | hifiasm read-layout support: supported / one_read_zone / mismatch / match_only / not_tested |
| reads_status | HG002 Q100 ends: supported (>= 3 ONT UL reads cross the switch, conflicting <= 10%) / weak / not_tested |
| hg002_ont_crossing, hg002_ont_conflict | ONT UL reads crossing the switch / conflicting (HG002 Q100 foreign ends) |
| ped_status | pedigree_aux child ends: intact / not_intact transmission; parent_or_untested; not_applicable |
| flags | comma-separated: artefact_candidate, junk_gt_1kb, junk_chimeric, short_telomere (< 1 kb), telomere_to_window_end, window_extended, partial_window, arm_inferred, arm_low_confidence, arm_conflict, complex_structure (>= 4 core blocks), dup_arm_in_haplotype, vk_discordant, one_read_zone, layout_mismatch |
| cis_S5_state, cis_S10_state, cis_S20_state, cis_S20_best_end | boundary-anchored states from analyses/telomere_cis/results/fixed_window_states.tsv (independent implementation; named Release 2 tips) |
| S20_vs_cis | concordant / state_differs / foreign_end_differs (named ends of non-pooled arms) |
| prior_t20_state, prior_t20_best_end, prior_form, prior_switch_final, prior_switch_sub, prior_junction_cluster, prior_junction_cluster_samples | earlier t20 foreign-end calls for named Release 2 tips (switch_final in tip frame; prior_switch_sub = switch_final - s0) |
| t20_vs_prior | t20 state vs the earlier t20 call: concordant / state_differs / foreign_end_differs (named Release 2 ends) |
| tip_key_r2 | Release 2 tip key `<assembly>|<contig>|<side>` used by the earlier tables |

## mosaic_blocks.tsv (one row per block; `level` raw or core)

| column | definition |
|---|---|
| end_id | end |
| level | raw (track blocks: runs of grid labels, runs < 2 kb merged, base-level boundaries) or core (blocks >= 10 kb) |
| block | block number from the telomere (1-based) |
| label, family | block label (`home`, family/arm token, `none`); family = the label of a foreign block |
| source | home blocks: the most frequent best home CHM13 end; other blocks: the label token (as used for core structures) |
| best_end | most frequent best single CHM13 end over the block's grid points (auxiliary) |
| source_arm | arm of source |
| is_home | block is home sequence |
| sub_start, sub_end | block span in subtelomere coordinates |
| tip_start, tip_end | block span in tip frame (raw) |
| n_grid_points, mean_tied | 100 bp grid points in the block; mean number of CHM13 ends within the tie margin per point (raw) |
| ident | mean best sliding 5 kb identity over the block's grid points (raw) |
| chm13_dist_start, chm13_dist_end | position on best_end (distance from its telomere) of the first and last matched base of the block (raw) |
| next_switch_res | resolution of the switch to the next block (base/grid) |

## alleles_structural.tsv (one row per structural allele)

| column | definition |
|---|---|
| allele_id | SA-<catalogue_arm>-NN, numbered by number of defining ends |
| arm | catalogue arm |
| structure | core structure (`home` = native; foreign tokens are paralog families `PFnn` or single-end arms) |
| form | `<arm>><first block label>` when the outermost core block is foreign |
| is_native | structure is exactly `home` |
| nativity_class | native / within_family / foreign, from the allele structure (native includes structures with `none` blocks besides home) |
| n_core_blocks | number of core blocks |
| switches_median, switches_min, switches_max | per core switch: median, min and max subtelomere coordinate among defining ends |
| tier | A / B / C (see below) |
| evidence | evidence lines present: recurrent (>= 3 samples), verkko, reads, pedigree, t2t_assembly |
| singleton | observed in one sample |
| n_ends | assigned ends (HPRC_R2 + CHM13 + HG002_Q100; full + partial + arm_lowconf) |
| n_defining | defining ends (full-length, reliable arm call) |
| n_full | assigned ends without the partial_window flag |
| n_r2_ends | assigned Release 2 ends |
| r2_freq | n_r2_ends / Release 2 ends of the arm with a structural allele call |
| n_samples | distinct samples among assigned ends (HG002 Q100 and HG002 Release 2 are one sample) |
| n_chm13, n_hg002_q100, n_pedigree_aux | members from these sets (pedigree ends are not counted in n_ends) |
| n_pedigree_intact | intact pedigree transmissions in which parent and child ends both carry the allele |
| n_vk_concordant, n_vk_discordant, n_layout_supported, n_hg002_reads_supported | members with that evidence status |
| n_artefact_candidates, n_flagged_junk, n_arm_inferred | members with those flags |
| representative_end_id | medoid defining end (smallest summed switch distance) |
| representative_in_fasta | the representative window is included in representatives.fa.gz (alleles seen in >= 3 samples) |
| haplogroups_top | up to five haplogroups of the members with counts |
| prior_junction_clusters | earlier junction clusters among members with counts |

Tiers. Evidence lines: R = present in >= 3 unrelated samples (no parent-child pairs exist among Release 2 samples in
samples.tsv); orthogonal O = verkko (a member's Verkko assembly agrees), reads (HG002 Q100 reads cross the member's switch, or the
hifiasm read layout supports the member's end), pedigree (intact transmission with parent and child carrying the allele),
t2t_assembly (a member in HG002 Q100, or in CHM13 for non-native alleles and haplogroups). **A** = R and >= 1 O, or >= 2 O.
**B** = R without O, or 2 samples with >= 1 O. **C** = all others (`singleton` marks one-sample alleles).

## haplogroups.tsv (one row per sequence haplogroup)

| column | definition |
|---|---|
| haplogroup_id | HG-<catalogue_arm>-NN, numbered by size |
| arm | catalogue arm |
| tier, evidence, singleton | as for structural alleles |
| n_clustered | ends clustered into the haplogroup (defining ends) |
| n_ends, n_full, n_partial, n_r2_ends, r2_freq, n_samples, n_chm13, n_hg002_q100, n_pedigree_aux, n_pedigree_intact, n_vk_concordant, n_layout_supported, n_hg002_reads_supported, n_artefact_candidates, n_arm_inferred | as for structural alleles |
| medoid_end_id | member with the highest mean ANI to the other clustered members |
| representative_in_fasta | the medoid window is included in representatives.fa.gz (haplogroups seen in >= 3 samples) |
| min_within_ani, mean_within_ani | minimum / mean pairwise k-mer ANI among clustered members |
| medoid_nearest_other_member_ani | highest ANI between the medoid and an end of another haplogroup of the arm |
| min_medoid_aln_ident, median_medoid_aln_ident, min_medoid_aln_frac | member-to-medoid alignment validation |
| structural_alleles_top | structural alleles of the members with counts |
| structural_purity | fraction of members carrying the most common structural allele |
| threshold_ani | UPGMA cut used (ANI) |
| compared_len | subtelomere length compared (bp from the telomere boundary) |

## allele_frequencies_by_superpop.tsv

| column | definition |
|---|---|
| allele_type | structural_allele or haplogroup |
| allele_id, arm | allele and catalogue arm |
| superpopulation | AFR, AMR, EAS, EUR, SAS, NA (no label) or ALL |
| n_haplotype_ends | Release 2 ends carrying the allele in the superpopulation |
| n_samples | distinct samples among them |
| n_arm_ends_assessed | Release 2 ends of the arm with a call of that type in the superpopulation |
| frequency | n_haplotype_ends / n_arm_ends_assessed |

## representatives.fa.gz

One record per distinct representative end of every structural allele and haplogroup observed in >= 3 samples (file-size
limit): the end's 100 kb tip-frame window (whole contig if shorter). Header `>ID ids=<all allele/haplogroup IDs it represents>
end_id=<end> len=<bp> subtelomere_start=<s0>`; the telomere occupies [0, s0). Any other window (including singletons) can be
retrieved and checksum-verified with `python pipeline/fetch_end_window.py <end_id>`.

## catalogue.json

`{version, window_bp, haplogroup_threshold_ani, arms: {<catalogue_arm>: {n_r2_ends, structural_alleles: [{id, structure, native,
tier, evidence, n_ends, n_samples, freq: {superpopulation: frequency}, switches, representative, blocks: [{source, start, end}],
members: [end_id]}], haplogroups: [{id, tier, n_ends, n_samples, freq, medoid, structural_alleles, members}]}}}`.
Block coordinates are subtelomere coordinates of the representative's core blocks.

## paralog_families.tsv

One row per multi-member paralog family: family (PFnn, numbered by size), n_members, members (arms of the CHM13 ends), n_edges, min/median
co-tie rate of its edges, cotie_points, tied_points, and the build parameters tau, min_points, margin, window_bp. Families are connected
components of the graph joining two CHM13 ends that are co-tied (identity within the margin of the best, sampled every 1 kb over all defining
ends) at >= min_points points and in >= tau of the points where the rarer end is tied. Ends not listed form single-end families named by their arm.

## id_history.tsv

One row per structural-allele (SA) or haplogroup (HG) ID: type, arm, id, status (unchanged / kept_changed / new / retired),
merged (the cluster received >= 2 previous IDs' members, >= 2 ends each), split (the previous ID's members went to >= 2 clusters),
n_members, previous_sources (previous IDs of the members with counts), now_in (where the previous ID's members are now).
IDs are inherited from the previous build by reciprocal best member overlap (pipeline/catalogue_08b_stable_ids.py); a structural-allele ID is
inherited only by an allele of the same nativity, so a native ID never passes to a non-native allele.

## qc/stability_summary.tsv, qc/stability_tau_scan.tsv

Committed copies of the boundary-shift acceptance tables (QC.md section 9). `stability_summary.tsv`: one row per boundary shift
(`shift`, bp, relative to the called boundary) and stratum (`stratum`; the stratum value is in `native` = previous-build nativity,
`new_native` = nativity of the allele at shift 0 in this build, `tier` = previous tier, or `arm_class` = duplicon_rich / other):
`n` ends and the fraction whose core structure, structural allele and S20 state are unchanged (`same_structure`, `same_allele`,
`same_S20`), each with its Wilson 95% interval as `lo-hi` (`*_ci`). `stability_tau_scan.tsv`: `same_structure` and `same_S20` by
paralog-family threshold `tau` (`n_families`, size of the `largest` family; tau 1.01 = no families), `shift` and previous-build
nativity (`native`), with `n` ends.

`public_window_verification.tsv` (end_id, set, stratum, source URL, window_len, status OK/MISMATCH/MISSING, selection) and
`public_rebuild_subset_check.tsv` (rebuild table, release column, n compared, n identical) / `.json` record the verification of
the public rebuild path (pipeline/CATALOGUE.md, 'Rebuild from public data').

## MANIFEST.tsv / SHA256SUMS

Path, size, SHA-256, row count and description of every released file; `SHA256SUMS` is `sha256sum -c` compatible from the repository root.
