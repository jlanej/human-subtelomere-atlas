### Revision 1 - telomere boundaries re-called on the full window

(Frozen text of the revision-1 QC report; numbers refer to the revision-1 build, which still used the 5 kb-bin mosaic.)

**v0.1 (revised before publication): telomere boundaries re-called on the full window.** In the first v0.1 build the telomere/subtelomere boundary of Release 2 ends was copied from the Release 2 tip table, whose telomere calls were made on 20 kb windows and therefore stop near 20 kb. All boundaries are now called on the catalogue window with the shared caller (`subtel_atlas.telomere.find_block`, same parameters). Ends whose boundary lies beyond 20 kb received an extended window (boundary + 80 kb, rounded up to 5 kb; range-read from the same public assemblies) so that they keep >= 80 kb of subtelomere. Bins were recomputed for the changed ends (a seeded control set of 300 unchanged ends reproduced its stored bin rows exactly); steps 05-10 were re-run in full. Allele and haplogroup IDs were carried over by reciprocal best member overlap (`id_history.tsv`).

Folded into the same revision: (i) ends with several compatible structural alleles now go to the allele with the nearest median switches, then the most defining ends (previously most defining ends first; section 8); (ii) new columns `boundary_stable`, `S20_state_minus500`, `S20_state_plus500` (section 8), `window_extended`, `telomere_to_window_end`, `tiptable_junk_len`, `tiptable_telo_len`; (iii) SCHEMA states that `haplotype` is the assembly label (parent of origin only when trio-phased) and documents boundary sensitivity as a known limitation; (iv) the query-tool package checksums and golden test outputs were regenerated from this release (the fixture windows themselves are unchanged).

| item | previous | revised |
|---|---|---|
| ends with a changed telomere boundary |  | 199 |
|   of which previous boundary >= 19 kb |  | 199 |
| ends with previous boundary >= 19 kb | 304 | 304 |
| windows extended beyond 100 kb | 0 | 195 |
| ends with S5 / S10 / S20 state changed |  | 57 / 39 / 14 |
| ends with core structure changed |  | 31 |
| ends with structural allele ID changed |  | 44 |
| ends with haplogroup ID changed |  | 29 |
| structural alleles | 1221 | 1215 |
| structural alleles tier A / B / C | 155/249/817 | 153/250/812 |
| haplogroups | 598 | 593 |
| haplogroups tier A / B / C | 130/151/317 | 131/152/310 |
| Release 2 ends with a structural allele / haplogroup | 20115 / 20111 | 20115 / 20111 |
| changed-boundary ends: S20 / core structure / allele / haplogroup changed |  | 14 / 31 / 25 / 13 of 199 |

Shared-code equivalence of step 07: the refactored haplogroup step, run in a scratch copy of the tree on the exact pre-revision inputs (windows reconstructed and SHA-256-verified), reproduced 6/6 original outputs byte for byte (haplogroup_members.tsv, haplogroups_work.tsv, hg_pair_ani_hist.tsv, hg_nn_ani.tsv, hg_calibration.tsv, hg_threshold_scan.tsv).

ID changes (`id_history.tsv`; merged = cluster receiving >= 2 previous IDs with >= 2 ends each, split = previous ID spread over >= 2 clusters):

| type | kept_changed | new | retired | unchanged |
|---|---|---|---|---|
| HG | 15 | 4 | 9 | 574 |
| SA | 44 | 8 | 14 | 1163 |

Merged: HG 1, SA 1; split: HG 2, SA 1.
