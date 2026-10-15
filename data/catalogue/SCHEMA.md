# Allele catalogue v0.1 - schema

Coordinates: **tip frame** = position 0 is the outermost base of the contig end, reading inward (telomere -> centromere); R ends
are reverse-complemented. **Subtelomere coordinate** (`sub_*`, `switches_sub`) = tip-frame position minus `s0`
(`s0 = junk_len + telo_len`, the masked telomere plus any outer non-telomeric bases). Contig window coordinates (`win_start`,
`win_end`) are 1-based inclusive on the deposited contig. CHM13 distances (`chm13_dist_*`) are bases from the telomeric end of
the named CHM13 v2.0 chromosome end. Arms are written `1p`, `20q`; CHM13 ends `chr20q`. Missing values are empty cells.

All bins, blocks, switches and structural calls are anchored at the telomere/subtelomere boundary (`s0`), so they do not depend
on telomere length. A **switch** (junction) is the position where the sequence changes which CHM13 end it most resembles; inside
shared, paralogous subtelomeric sequence this is a change in paralog similarity, not necessarily a structural breakpoint.
'Foreign' is a distance-dependent property: see `profile` (5 kb bins from the boundary) and the S5/S10/S20 states.

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
| haplotype | haplotype label of the assembly (hap1/hap2/mat/pat; CHM13 haploid; HG002 pat/mat; pedigree h1/h2/mat/pat) |
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
| window_len | window length (bp); 100,000 unless the contig is shorter |
| window_sha256 | SHA-256 of the upper-case tip-frame window sequence (no header/newlines) |
| junk_len | bases outside the telomere at the very tip (before the first telomeric repeat block) |
| telo_len | length of the terminal (CCCTAA)n block in tip frame |
| telo_source | tip_table (Release 2 tip table, same caller) or catalogue_telomere (pipeline/catalogue_telomere.py) |
| home_arm | the end's own arm: from the contig name, or inferred (see arm_source) |
| catalogue_arm | arm or pooled group used for alleles/haplogroups |
| arm_source | contig_name, probe (10 kb probes at 250 kb / 1 Mb mapped to CHM13), mosaic (best CHM13 end of the proximal window), none |
| arm_confidence | named; high (two probes agree); medium (one informative probe); conflict (probes disagree at group level); mosaic_group (mosaic call inside a pooled group); low (other mosaic calls); none |
| probe_arm, mosaic_arm | arm proposed by the probe rule / the mosaic rule (recorded for all ends with data) |
| profile | label of each 5 kb bin counted from the telomere boundary: O = home, F = another CHM13 end, . = no alignment (after smoothing) |
| structure | raw block sources telomere -> centromere (`home`, CHM13 arm of a foreign block, `none`), `|`-separated |
| switches_sub | raw block boundaries (subtelomere coordinates), comma-separated |
| switch_res | per raw boundary: base (base-level changepoint) or window (bin boundary, +-5 kb) |
| core_structure | structure after dropping blocks < 10 kb (defines structural alleles) |
| core_switches_sub | core block boundaries (subtelomere coordinates) |
| n_core_blocks | number of core blocks |
| sub_len | subtelomere length available in the window (window_len - s0) |
| S5_state, S10_state, S20_state | boundary-anchored fixed-window states of the first 5 / 10 / 20 kb of subtelomere: own (a home CHM13 end scores best), tie (home within 2% of the best foreign score), foreign, none (no alignment covering >= 50%) |
| S5_best_end, S10_best_end, S20_best_end, S20_best_ident | best CHM13 end of each window (and S20 gap-compressed identity) |
| form | `<home_arm>><foreign arm>` for S20-foreign ends |
| t20_state, t20_best_end, t20_best_ident | the earlier definition, tip[s0:20 kb] (shrinks as the telomere grows); kept only for comparison with earlier calls |
| structural_allele_id | SA-<catalogue_arm>-NN, empty if unassigned |
| sa_assignment | full (defines the allele), partial_unique / partial_ambiguous (short window, prefix-compatible), arm_lowconf, aux (pedigree), unassigned_* |
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
| flags | comma-separated: artefact_candidate, junk_gt_1kb, junk_chimeric, short_telomere (< 1 kb), partial_window, arm_inferred, arm_low_confidence, arm_conflict, complex_structure (>= 4 core blocks), dup_arm_in_haplotype, vk_discordant, one_read_zone, layout_mismatch |
| cis_S5_state, cis_S10_state, cis_S20_state, cis_S20_best_end | boundary-anchored states from analyses/telomere_cis/results/fixed_window_states.tsv (independent implementation; named Release 2 tips) |
| S20_vs_cis | concordant / state_differs / foreign_end_differs (named ends of non-pooled arms) |
| prior_t20_state, prior_t20_best_end, prior_form, prior_switch_final, prior_switch_sub, prior_junction_cluster, prior_junction_cluster_samples | earlier t20 foreign-end calls for named Release 2 tips (switch_final in tip frame; prior_switch_sub = switch_final - s0) |
| t20_vs_prior | t20 state vs the earlier t20 call: concordant / state_differs / foreign_end_differs (named Release 2 ends) |
| tip_key_r2 | Release 2 tip key `<assembly>|<contig>|<side>` used by the earlier tables |

## mosaic_blocks.tsv (one row per block; `level` raw or core)

| column | definition |
|---|---|
| end_id | end |
| level | raw (5 kb-bin blocks with base-level switches) or core (blocks >= 10 kb) |
| block | block number from the telomere (1-based) |
| source | CHM13 end best matching the block (for home blocks the best-scoring home end), or none |
| source_arm | arm of source |
| is_home | block is home sequence |
| sub_start, sub_end | block span in subtelomere coordinates |
| tip_start, tip_end | block span in tip frame (raw) |
| n_bins | 5 kb bins in the block (raw) |
| ident | mean gap-compressed identity of the source alignment over the block's bins (raw) |
| chm13_dist_start, chm13_dist_end | span of the source alignments on the source CHM13 end, distance from its telomere (raw) |
| next_switch_res | resolution of the switch to the next block (base/window) |
| next_sw_outer_chm13_dist, next_sw_inner_chm13_dist | CHM13 coordinates (distance from telomere) of the last base before / first base after a base-level switch |
| next_sw_outer_ident_5kb, next_sw_inner_ident_5kb | identity of the outer/inner source within 5 kb of the switch |

## alleles_structural.tsv (one row per structural allele)

| column | definition |
|---|---|
| allele_id | SA-<catalogue_arm>-NN, numbered by number of defining ends |
| arm | catalogue arm |
| structure | core structure (`home` = native) |
| form | `<arm>><first block source>` when the outermost block is foreign |
| is_native | structure is `home` |
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

## MANIFEST.tsv / SHA256SUMS

Path, size, SHA-256, row count and description of every released file; `SHA256SUMS` is `sha256sum -c` compatible from the repository root.
