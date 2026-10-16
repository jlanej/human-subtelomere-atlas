# subtel-atlas: query tool for the human subtelomere allele catalogue

`subtel-atlas` places telomeric chromosome ends of a new assembly (or pre-cut telomere-first windows, or telomeric long reads)
into catalogue release **v0.1, revision 2** (`data/catalogue/`). For every telomeric end it reports:

* the catalogue arm;
* the end's mosaic of CHM13 v2.0 end segments, labelled by paralog family;
* the S5/S10/S20 terminal states;
* the nearest structural allele (`SA-<arm>-NN`), with its nativity class, evidence tier and population frequency;
* the nearest sequence haplogroup (`HG-<arm>-NN`), with identities;
* whether the allele and S20 calls survive a ±500 bp shift of the telomere boundary;
* novelty flags.

The tool runs the catalogue's own code. The package (`atlas/subtel_atlas`) contains the telomere caller, the window rule, the
track-based mosaic and paralog-family labels (`tracks.py`), and the structure, allele-matching, arm and haplogroup rules. The
catalogue build (`pipeline/catalogue_*.py`) imports them from here. Re-querying 5,981 catalogue windows reproduces every
mosaic, S5/S10/S20 state, core structure and boundary-stability call in `ends.tsv` (`tests/validation/VALIDATION.md`).

The atlas is one of several concurrent catalogues of human subtelomere variation. It complements Tan et al. 2026 (bioRxiv,
doi:10.64898/2026.09.03.748401; 860 haplotypes, 500 kb ends, paralogy-block haplotypes) and Guarracino et al. 2026 (bioRxiv,
doi:10.64898/2026.07.10.737660). Allele, haplogroup and family IDs here are specific to catalogue v0.1 and do not correspond
to the labels used in those studies.

## Install

```bash
pip install "git+https://github.com/jlanej/human-subtelomere-atlas"     # or, from a checkout: pip install .
conda install -c bioconda minimap2                                      # recommended (exact catalogue computation)
subtel-atlas fetch-refs                                                 # ~40 MB into ~/.cache/subtel-atlas/v0.1
```

Python dependencies: `mappy`, `numpy`, `pandas`; `pysam` is optional, for BAM/CRAM input to `reads`. With the `minimap2`
executable (>= 2.24; validated with 2.31) on `PATH`, the alignments are computed exactly as in the catalogue build. Without it
the `mappy` binding is used, which is an approximation (see Limits).

`fetch-refs` builds the reference directory:

* **CHM13 ends.** The outer 500 kb of the 48 CHM13 v2.0 chromosome ends, by HTTP range reads from the public FASTA
  (`https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa`). `--chm13 <FASTA>` reads a local copy instead.
* **Catalogue files.** Copied from the repository checkout, or downloaded from `--catalogue-url` (the repository's raw URL once
  published): `alleles_structural.tsv`, `haplogroups.tsv`, `allele_frequencies_by_superpop.tsv`, `representatives.fa.gz`,
  `paralog_families.tsv` and `PROVENANCE.json`.
* **Verification.** Every file is checked against SHA-256 values shipped with the package (one per CHM13 end, and the release
  checksums of the catalogue files). `subtel-atlas verify-refs` repeats the check. Each query re-verifies the directory;
  `--no-verify-refs` skips this for custom references.

## Usage

```bash
subtel-atlas query assembly.fa.gz -t 4 > ends.tsv              # every telomeric contig end of an assembly (FASTA/FASTQ, gzip ok)
subtel-atlas query assembly.fa.gz --format json -o ends.json  # + blocks, core blocks, top-3 haplogroups, frequencies, metadata
subtel-atlas query --windows tips.fa > tips.tsv               # pre-cut tip-frame windows (first base = contig tip)
subtel-atlas query assembly.fa --chm13-genome chm13v2.0.fa    # arm probes for unnamed contigs (catalogue probe rule)
subtel-atlas reads telomeric_reads.fq.gz > reads.tsv          # EXPERIMENTAL: telomeric long reads (FASTQ/FASTA, or BAM with pysam)
```

* **Contig mode** (default). Each side of each contig is a candidate end, read in tip frame (R ends reverse-complemented). The
  end is reported if the catalogue telomere caller finds a clean (CCCTAA)n block of at least `--min-telomere` bp starting within
  8 kb of the tip (default 500 bp, the Release 2 telomeric-tip definition).
* **Window length.** The window is the outer 100 kb. If the telomere/subtelomere boundary s0 lies beyond 20 kb, the window is
  extended to ceil((s0 + 80 kb) / 5 kb) × 5 kb so that a full-length end keeps 80 kb of subtelomere; this is the catalogue rule,
  reported in `window_extended`. Coordinates (`win_start`, `win_end`) are 1-based on the contig, as in `ends.tsv`.
* **Windows mode** (`--windows`). Each record is already a tip-frame window; the telomere is called on the window. Header tokens
  `arm=20q`, `side=L|R` and `junk_len=<bp> telo_len=<bp>` (a boundary from another caller) are used when present.
* **Arms.** For chromosome-named contigs (`chr7`, `chr7_MATERNAL`, `HG002#1#chr7`), the L end is the p arm and the R end the q
  arm (standard orientation assumed). `--arm-map` supplies arms explicitly; `--no-name-arms` disables name parsing. For unnamed
  contigs, `--chm13-genome` (the CHM13 v2.0 FASTA or a minimap2 `-x asm5` index) applies the catalogue probe rule: 10 kb probes
  at 250 kb and 1 Mb from the tip, so contigs of at least 520 kb for one probe and about 2 Mb for both. Without probes the arm
  comes from the bin-based mosaic rule and is flagged `arm_inferred` (see Limits).
* `-t` sets threads for minimap2 and worker processes. Output is identical for any thread count, and rows follow input order.

## How an end is described

1. **Boundary.** `s0` = outer junk + telomere, from the shared caller. Every downstream coordinate (`sub_*`, switches) counts
   from s0.
2. **Identity tracks.** The subtelomere (window from s0 on) is aligned once to the 48 CHM13 ends, keeping secondary chains.
   For each CHM13 end, a per-base match track gives the identity over a sliding 5 kb, evaluated on a 100 bp grid anchored at
   the contig tip.
3. **Labels.** At each grid point:
   * `home` if a CHM13 end of the end's own arm (or pooled group) is within 0.5 percentage points of the best identity;
   * `none` if no end reaches 80%;
   * otherwise the **paralog family** of the best end.
4. **Blocks and structure.** Runs of a label form blocks. Runs shorter than 2 kb are merged into a neighbour, and boundaries
   are refined at base level. Blocks shorter than 10 kb are then dropped to give the **core structure** (e.g. `PF04|home`),
   which defines structural alleles.
5. **Fixed windows.** S5/S10/S20 take the majority label over the first 5/10/20 kb: own, tie, foreign (with the family in
   `S*_family`) or none. `S*_best_end` keeps the single best CHM13 end.
6. **Allele.** The end matches an allele when it has the same core structure and every switch lies within 10 kb of the
   allele's median. Several candidates: nearest medians first, then most defining ends.
7. **Haplogroup.** k-mer ANI over the first 80 kb of subtelomere against the medoid of each represented haplogroup of the arm,
   with an alignment check against the best one.
8. **Boundary stability.** Steps 3-6 are repeated with the boundary moved by -500 and +500 bp (identity tracks shifted, no
   re-alignment; as `pipeline/catalogue_08e_boundary_stable.py`).

**Paralog families** (`paralog_families.tsv`). Some CHM13 ends share long, near-identical duplicated blocks, so the single
best end can flip between them with small sequence or boundary changes. Two CHM13 ends are joined into a family when, across
the catalogue, they are tied (within 0.5 points) at >= 75% of the points where the less often tied of the two is tied. Release
v0.1 has six multi-member families:

| family | CHM13 ends |
|---|---|
| PF01 | 4p, 14p, 15p, 21p, 22p |
| PF02 | 1p, 5q, 6q |
| PF03 | 11p, 15q |
| PF04 | 3q, 19p |
| PF05 | Xp, Yp |
| PF06 | Xq, Yq |

Any other CHM13 end is its own family, named by its arm (e.g. `12p`). A family label says which shared block the sequence
belongs to; `S*_best_end` and the block `source` keep the best single end, and `mean_tied` gives the local number of tied ends.

**Nativity class** (`sa_nativity_class`, from the matched allele):

* `native`: only home (and none) blocks.
* `within_family`: every non-home block carries a family that contains the end's own arm. The sequence is closer to another
  member's CHM13 copy of a shared block than to the own-arm copy.
* `foreign`: at least one block is labelled with a family or arm that does not contain the own arm.

### Example (real output)

Four catalogue windows (`tests/fixtures/windows.fa.gz`) queried with the full references (`--windows`, minimap2 backend;
selected columns):

| query_id | arm | s0 | S20 (family; best end) | form | core_structure | core switches | allele | nativity | assignment | tier | R2 freq | boundary_stable | haplogroup | HG ANI | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HG02129#2#CM085868.1#R | 20q | 7177 | foreign (PF04; chr3q) | 20q>PF04 | PF04\|home | 31040 | SA-20q-01 | foreign | full_match | A | 0.4829 | True | HG-20q-04 | 0.99942 | sa_multiple_candidates |
| NA21093#2#JBHIJD010000010.1#L | 6p | 1618 | foreign (PF04; chr3q) | 6p>PF04 | PF03\|PF04\|20q\|12p\|20q | 10938,30427,59740,70219 | SA-6p-01 | foreign | full_match | A | 0.1936 | True | HG-6p-01 | 0.99935 | complex_structure |
| CHM13#0#chr1#L | 1p | 889 | own (chr1p) | | home | | SA-1p-01 | native | full_match | A | 0.1218 | True | HG-1p-02 | 0.99856 | short_telomere |
| HG01071#2#JAHBCE020000040.1#R | 20q | 4704 | foreign (PF04; chr3q) | 20q>PF04 | PF04\|home | 32661 | SA-20q-01 | foreign | partial_ambiguous | A | 0.4829 | True | HG-20q-01 | 0.99815 | hg_close_second, partial_window |

In JSON, the first end has these blocks:
* `PF04` (best end chr3q) over subtelomere 0-31,040 bp, mean identity 0.9981, 1.0 ends tied on average;
* `home` (chr20q) over 31,040-92,823 bp.

It also carries `hg_top` = HG-20q-04 0.99942, HG-20q-01 0.99312, HG-20q-02 0.98824, and the SA-20q-01 Release 2 frequency
by superpopulation: AFR 0.785, AMR 0.383, EAS 0.194, EUR 0.574, SAS 0.352, all 0.483.

Whole-genome contig-mode runs: CHM13 v2.0 (48 telomeric ends) took 42 s (2.2 GB peak memory) with `-t 3`; diploid HG002 Q100 v1.1 (86 ends)
took 78 s (2.4 GB).

## Output columns (TSV; JSON has the same fields plus nested detail)

| column | meaning |
|---|---|
| query_id, input_name, mode, side, contig_len, win_start, win_end, strand | end identity (`<contig>#L/R` in contig mode) and window on the input |
| window_len, window_extended | window length; extended beyond 100 kb by the catalogue rule (boundary beyond 20 kb) |
| telomere_found, junk_len, telo_len, telomere_source, s0, sub_len, full_length | telomere call (s0 = junk + telomere = boundary); subtelomere length in the window; full = window >= 100 kb with >= 80 kb subtelomere |
| home_arm, catalogue_arm, arm_source, arm_confidence | arm and pooled catalogue group (acroP = 13p/14p/15p/21p/22p, PAR1 = Xp/Yp, PAR2 = Xq/Yq); source contig_name / header / arm_map / probe / mosaic; confidence named / high / medium / conflict / mosaic_group / low / none |
| probe_arm, probe_confidence, mosaic_arm, mosaic_arm_score_frac | the probe-rule and mosaic-rule calls (the mosaic rule is always computed) |
| profile | track label every 1 kb of subtelomere: O = home, F = foreign (a family or single CHM13 end), . = none |
| structure, switches_sub, switch_res | raw blocks telomere -> centromere (`home`, `PFnn`, an arm, `none`), their boundaries (bp from s0) and resolution (base / grid) |
| core_structure, core_switches_sub, n_core_blocks | blocks >= 10 kb; this defines structural alleles |
| S5/S10/S20_state, S5/S10/S20_family, S5/S10/S20_best_end, S20_best_ident | first 5/10/20 kb of subtelomere: own / tie / foreign / none (majority grid label); family of a foreign window; most frequent best single CHM13 end; mean best identity |
| form | `<arm>><S20 family>` for S20-foreign ends |
| boundary_stable, boundary_stable_S20, boundary_stable_allele | the S20 call (state + family) and/or the structural allele are unchanged when the boundary moves by -500 and +500 bp |
| S20_state_minus500, S20_state_plus500, allele_minus500, allele_plus500 | the calls at the shifted boundaries |
| structural_allele_id, sa_is_native, sa_nativity_class, sa_assignment, sa_candidates, sa_max_switch_delta | matched allele and its nativity (native / within_family / foreign, when the catalogue provides the class); full_match, partial_unique / partial_ambiguous (short windows: the allele truncated at the observed length matches), unassigned, not_assessed; all compatible alleles; largest switch difference (bp) |
| sa_nearest_same_structure, sa_nearest_delta | closest allele with the same structure even when outside tolerance (novel switch positions) |
| sa_tier, sa_evidence, sa_n_samples, sa_r2_freq, sa_representative_end_id | catalogue evidence tier (A/B/C), evidence lines, samples, Release 2 frequency on the arm, representative end |
| haplogroup_id, hg_status, hg_ani, hg_second_id, hg_second_ani, hg_compared_len | nearest represented haplogroup medoid by k-mer ANI over the first 80 kb of subtelomere (or the available length); assigned if ANI >= 0.995 |
| hg_aln_ident, hg_aln_frac, hg_tier, hg_n_samples, hg_r2_freq | alignment check against that medoid (gap-compressed identity, aligned fraction); catalogue evidence for the haplogroup |
| novelty_flags | sa_novel_structure (no allele of the arm has this core structure), sa_novel_switch_position (structure known, no allele within 10 kb), sa_no_compatible_allele (partial window), hg_below_threshold, hg_close_second (second medoid also >= 0.995 and within 0.001), hg_alignment_check_failed (identity < 0.99 or aligned fraction < 0.90, about the 1st percentile of catalogue members), arm_unassigned |
| qc_flags | telomere_not_found, short_telomere (< 1 kb), junk_gt_1kb, partial_window, arm_inferred, arm_low_confidence, arm_conflict, complex_structure (>= 4 core blocks), sa_multiple_candidates, boundary_unstable, experimental_read_mode, mappy_backend_approximation |
| aligner, tool_version, catalogue_version | provenance of every row |

## Interpretation guide

* **Labels are relative to CHM13.** 'Foreign' telomere-adjacent sequence most resembles another CHM13 end (family). A switch is
  where that resemblance changes, inside shared paralogous subtelomeric sequence; it is not necessarily a structural breakpoint.
  The foreign fraction depends on the window: 17.4% of Release 2 ends are S20-foreign in revision 2, and other windows give
  other fractions (QC.md). Report the window used.
* **Family-level forms.** A form such as `20q>PF04` means that the outer block belongs to the 3q/19p shared block. `S20_best_end`
  says which member it most resembles (chr3q for SA-20q-01).
* **Matched alleles.** A matched structural allele means the same core block order as a catalogued allele, with every switch
  within 10 kb of its median. With several candidates (`sa_multiple_candidates`), treat the call as the set in `sa_candidates`.
  The tier is the catalogue's evidence for the allele, not for your end:
  * A: recurrent in >= 3 samples with orthogonal support, or >= 2 orthogonal lines;
  * B: recurrent only, or 2 samples with support;
  * C: all others, including singletons.
* **Boundary instability.** `boundary_unstable` marks calls that change when the boundary moves by 500 bp. These are typically
  ends whose first switch is near a block-length threshold. Prefer family-level or S20 summaries for such ends.
* **Novelty flags are relative to catalogue v0.1.**
  * `sa_novel_*`: no catalogued allele, singletons included, fits within tolerance.
  * `hg_below_threshold`: no *represented* haplogroup medoid reaches ANI 0.995. `representatives.fa.gz` holds medoids only for
    haplogroups seen in >= 3 samples, so an end of a rarer catalogued haplogroup is also reported as below threshold.
* **Haplogroups are compared by medoid.** The catalogue uses UPGMA clusters; the tool compares one medoid per haplogroup, which
  is less precise where haplogroups are diffuse (acroP, 10q, PAR1). See `hg_second_*` and `hg_close_second`.
* Frequencies are Release 2 haplotype-end frequencies; the samples are not a population-based sample.

## Accuracy (tests/validation/VALIDATION.md; minimap2 backend)

| test | ends | structural allele | haplogroup | mosaic, S5-S20, core structure, ±500 bp stability calls |
|---|---|---|---|---|
| (a) re-query of catalogue windows | 5,981 | 99.82% (100% with one compatible allele) | 94.82% (97.91% for represented haplogroups) | 100% |
| (b) leave-one-sample-out, 60 samples | 5,293 | 99.77% (removed alleles reported novel 74/76) | 97.71% (100/111) | as (a) |
| (c) CHM13 v2.0 whole assembly | 48 | 100% | 93.75% | 100% |
| (c) HG002 Q100 whole assembly | 86 | 100% | 96.51% | 100% |
| (c) WashU pedigree ends (never used for definitions) | 366 | 100% | 91.26% | 100% |
| (d) boundary moved ±500 bp, re-aligned (600 ends) | 600 | stable 97.17% / 96.83% (native 99.67 / 100%, foreign 94.84 / 93.65%) | stable 99.67% / 99.17% | S20 stable 95.67% / 96.50% |

In the pedigree, child and transmitting parent receive the same structural allele in 128/132 intact transmissions and the same
haplogroup in 129/132 (catalogue calls: 128/132 and 131/132).

## Limits

* **Arm inference without names or probes** uses the bin-based mosaic rule. It gives the correct catalogue group for 84% of named
  ends. Use `--chm13-genome` (the probe rule reproduces the catalogue probe calls for 600/600 tested ends) or `--arm-map`
  whenever possible.
* **Boundary sensitivity.** Moving the boundary by 500 bp changes 5-6% of foreign and within-family allele calls and < 0.5% of
  native ones (re-aligned test, VALIDATION.md (d)). The per-end ±500 bp flags are computed from shifted tracks without
  re-alignment: they agree with re-aligned calls for 99.2% of ends; they catch 26 of 32 unstable ends with no false alarms, and 6 of the 574 ends
  they call stable changed.
* **Tied alignments.** minimap2 breaks ties between equally scoring alignments using a hash of the query name. Queries are
  therefore named `<end_id>@0`, as in the catalogue build. Renaming an input can still change rare tied calls. For a given
  input, output is deterministic.
* **mappy backend** (no minimap2 executable). Under the track method it agrees with the minimap2 backend for S20 93.0%, core structure 92.0%, structural allele 91.5% and haplogroup 100% (400 ends), and runs at about 4 s per end per thread; the CLI prints a warning. Install minimap2 for catalogue-exact results. Rows are flagged `mappy_backend_approximation`.
* **Read mode is experimental.** It was validated only on HG002 read slices (`tests/validation/reads_hg002_summary.tsv`). Reads
  shorter than about 50 kb give partial structures, and read errors lower identities and ANI.
* Ends with less than 500 bp of clean telomere are not reported in contig mode; `--min-telomere` lowers the threshold.
* Results are versioned against catalogue v0.1 (revision 2). IDs are not stable across catalogue releases; see
  `data/catalogue/id_history.tsv`.

## Python API

```python
from subtel_atlas import query as Q
ref = Q.Reference()                                   # verified reference directory
ends = Q.contig_ends('chr20', seq)                    # or Q.window_end(name, header_comment, seq) / Q.read_end(name, seq)
records = Q.analyse(ends, ref, backend='minimap2', threads=4)
```

## Tests and validation

* `python -m pytest tests/` runs offline in about 15 s on small fixtures (`tests/fixtures`, 0.34 MB). Set
  `SUBTEL_ATLAS_REFS` to a fetch-refs directory to add a comparison with `ends.tsv`.
* After a catalogue rebuild, `python tests/regenerate_golden.py` refreshes the golden outputs.
* The full validation is `tests/validation/run_validation.py`, with stages select, query, compare, loso, shift, indep, probes,
  reads and mappy_backend; see `tests/validation/VALIDATION.md`.
