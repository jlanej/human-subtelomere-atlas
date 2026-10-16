# Validation of the subtel-atlas query tool (tool 0.1.0, catalogue v0.1 revision 2)

All numbers are computed by `tests/validation/run_validation.py` and read from the tables in this directory.
* **Catalogue.** Revision 2: track-based mosaic with paralog families and nativity classes (`data/catalogue/QC.md`, CHANGELOG).
  The files were verified unchanged from the start to the end of the validation (`SHA256SUMS` fingerprint).
* **Tool runs.** minimap2 2.31 backend; 4 threads for windows mode, 3 for whole-assembly runs.
* **References.** From `subtel-atlas fetch-refs`: the 48 CHM13 v2.0 ends range-read from the public FASTA, every file SHA-256-verified.

## 0. One implementation for catalogue and tool (`pipeline_equivalence.tsv`)

The package contains every method used by the catalogue build:
* the telomere caller and window rule;
* the identity-track mosaic and paralog families (`tracks.py`, imported by pipeline steps 04t, 06t and 08e);
* the structure and allele-matching rules (`alleles.py`, step 08);
* the arm rules (`arms.py`, step 05);
* the haplogroup sketch and alignment (`haplo.py`, `sketch.py`, step 07).

`tests/test_query_tool.py::test_pipeline_imports_shared_code` checks that the pipeline imports them.

Equivalence is shown end to end. The tool runs its own minimap2 alignment and reproduces the catalogue tables in sections (a)
and (c). The byte-level checks of the original refactoring, made on the bin-based first build, are listed as historical rows.

Review of the revision-2 changes made to the package by the catalogue build (`tracks.py`; edits to `query.py`, `refs.py`,
`alleles.py`):
* The track method and allele tie-break are consistent with the pipeline.
* Three gaps were fixed in `query.py`:
  1. The catalogue's window extension for boundaries beyond 20 kb was missing.
  2. Track queries are now named `<end_id>@0`, as in step 04t, because minimap2 breaks score ties with the query name.
  3. Ends without alignments now get a mosaic, as in step 06t.
* Added per-end outputs: the ±500 bp boundary-stability calls (same computation as step 08e), `window_extended`,
  `sa_is_native` and `sa_nativity_class`.
* The mappy fallback was re-validated; it is a weaker approximation under the track method (section (e)).

## (a) Self-consistency: re-querying catalogue windows (`self_consistency*.tsv`)

**Query set** (`query_set.tsv`; seed 20261015): 5,981 ends from 43 arm groups. It contains:
* every end of 60 Release 2 samples, 12 per superpopulation (5,293 ends);
* up to 3 ends per arm and structural-allele tier and 2 per arm and haplogroup tier (508);
* 40 ends with a boundary beyond 20 kb (76 extended windows in total);
* all CHM13 and HG002 Q100 ends (140).

The tool re-called each telomere on the window. The catalogue arm was supplied, so that mosaic and assignment are compared at
the same arm; the mosaic-rule arm was computed independently.

| metric | agreement (n = 5,981) |
|---|---|
| telomere boundary s0; window length (incl. 76 extended windows) | 100%; 100% |
| S5 / S10 / S20 state and family | 100% / 100% / 100% |
| form, 1 kb profile, core structure, core switch positions | 100% each |
| ±500 bp boundary-stability calls (`boundary_stable_allele`, `allele_minus500`, `allele_plus500`) | 100% |
| per-end nativity class | 100% |
| structural allele | 99.82% (5,970) |
| haplogroup | 94.82% (5,671) |
| mosaic-rule arm = catalogue mosaic call; mosaic-rule arm group = true arm (named ends) | 100%; 83.93% (3,175/3,783) |

Structural alleles in more detail:
* **By compatible alleles.** With one compatible allele the tool matches the catalogue for 5,547/5,547 ends. All 11
  differences are ends with 2-3 compatible alleles (17q, 17p, 20q, 5q). The catalogue's clustering placed these defining ends
  in an allele other than the one with the nearest median.
* **By tier.** A 99.83% (5,195/5,204), B 99.64% (554/556), C 100% (215/215).
* **By nativity.** Native 100% (4,230), within_family 99.13% (227/229), foreign 99.41% (1,507/1,516).
* **By catalogue boundary stability.** Boundary-stable alleles 99.98% (5,812/5,813); boundary-unstable alleles 94.05% (158/168).

Haplogroups in more detail:
* **Represented haplogroups** (with a medoid in `representatives.fa.gz`): 97.91% (5,666/5,787). The weakest arms are acroP
  (88.7%), 10q (90.5%) and PAR1 (91.2%).
* **Unrepresented haplogroups** (1-2 samples): 189 ends, all tier C. 170 are reported `hg_below_threshold`; 19 are absorbed
  into another haplogroup.

## (b) Leave-one-sample-out (`loso*.tsv`)

For each of the 60 samples (5,293 Release 2 ends), the sample's ends were removed from the references before assignment.
* **Structural alleles.** Re-clustered from the remaining defining ends with the catalogue rule (complete linkage, 10 kb
  Chebyshev cut); each cluster keeps the ID of its majority original allele.
* **Haplogroups.** A medoid from the held-out sample was replaced by the remaining clustered member with the highest mean ANI.
  Haplogroups left with fewer than 3 samples were dropped from the representatives. Haplogroups were not re-clustered (see
  Deviations).

Expected answer: the catalogue allele if it still exists without the sample, otherwise "novel".

| | n | accuracy | correct assignment | correct novel | false novel | wrong allele | assigned to other |
|---|---|---|---|---|---|---|---|
| structural allele, all | 5,293 | 99.77% | 5,203 | 78 | 0 | 10 | 2 |
| tier A / B / C | 4,767 / 418 / 104 | 99.83 / 99.28 / 99.04% | | | | | |
| allele removed by the hold-out | 76 | 97.37% | 0 | 74 | 0 | 0 | 2 |
| haplogroup, all | 5,293 | 97.71% | 5,067 | 105 | 18 | 92 | 11 |
| haplogroup removed by the hold-out | 111 | 90.09% | 0 | 100 | 0 | 0 | 11 |

* Hold-out errors largely coincide with the self-consistency errors: 8/10 structural-allele errors and 110/110 haplogroup
  errors are also wrong with the full reference.
* Accuracy by superpopulation: structural allele 99.44-100%, haplogroup 96.59-98.83% (AFR lowest).
* The two T2T genomes are themselves catalogue members, so they were also held out:
  * CHM13: structural allele 48/48, haplogroup 46/48;
  * HG002: structural allele 86/86 (2 correctly novel), haplogroup 83/86.

## (c) Independent assemblies, contig mode (`independent_*.tsv`, `runtime.tsv`)

Whole input files went through the CLI. This also tests telomeric-end detection, orientation, contig-name arms, the window rule
and coordinates. Input names differ from the catalogue IDs, so the minimap2 queries were named differently from the catalogue
build.

| dataset | telomeric ends found / catalogue ends | boundary, window | S5-S20, form, core structure | structural allele | boundary-stability calls | haplogroup |
|---|---|---|---|---|---|---|
| CHM13 v2.0 whole assembly | 48 / 48 | 100% | 100% | 100% | 100% | 93.75% |
| HG002 Q100 v1.1 whole diploid assembly | 86 / 92 | 100% | 100% | 100% | 100% | 96.51% |
| WashU pedigree v1.1 outer 250 kb ends (never in definitions) | 366 / 368 | 100% | 100% | 100% | 100% | 91.26% |

* **Ends not found.** The 8 missing ends (6 HG002, 2 pedigree) have a terminal (CCCTAA)n block of 247-413 bp, below the
  500 bp telomeric-end rule; the catalogue keeps T2T and pedigree ends regardless of telomere length.
* **Pedigree transmissions** (`independent_pedigree_transmissions.tsv`). In intact parent-to-child transmissions, child and
  transmitting parent receive:
  * the same structural allele in 128/132 (catalogue calls: 128/132);
  * the same haplogroup in 129/132 (131/132);
  * the same S20 call in 132/132.
* **No HG002 Release 2 comparison.** HG002 has no Release 2 ends in catalogue v0.1, so it is compared with the catalogue's
  HG002 Q100 entries (with the hold-out above).
* **Arm probe rule** (`probe_rule.tsv`; arm method unchanged since revision 1, table from that run): the tool reproduces the
  catalogue's probe arm and confidence for 600/600 ends with the minimap2 executable, and 560/600 with mappy.

## (d) Boundary-shift robustness of query results (`boundary_shift*.tsv`)

600 full-length ends from the query set: 300 with non-native alleles (252 foreign, 48 within_family) and 300 native.
Each perturbed input was re-queried with full re-alignment and compared with the unperturbed query.

**(A) Boundary moved by ±500 bp.** The header boundary override puts s0 500 bp outward (into the telomere) or inward:

| shift | n | S20 call | core structure | structural allele | haplogroup | allele = tool's ±500 prediction |
|---|---|---|---|---|---|---|
| -500 bp, all | 600 | 95.67% | 97.67% | 97.17% | 99.67% | 99.17% |
| -500 bp: native / within_family / foreign | 300 / 48 / 252 | 99.67 / 95.83 / 90.87% | | 99.67 / 93.75 / 94.84% | | |
| +500 bp, all | 600 | 96.50% | 97.33% | 96.83% | 99.17% | 99.17% |
| +500 bp: native / within_family / foreign | 300 / 48 / 252 | 99.67 / 91.67 / 93.65% | | 100 / 93.75 / 93.65% | | |

The per-end flag `boundary_stable_allele`, computed without re-alignment, marked 26 of the 32 ends whose allele changed under a
re-aligned ±500 bp shift (81% sensitivity), and all 26 ends it flagged unstable were unstable (no false alarms). Of the 574 ends
it called stable, 568 (99.0%) kept their allele and 6 changed.

**(B) Telomere lengthened or shortened.** (CCCTAA)n was inserted at the start of the telomere (+250, +500 bp), or telomere
sequence was trimmed (-250, -500 bp, ends with telomere >= 1.25 / 1.5 kb). The boundary was re-called by the tool, and the
window was kept at its length:

| perturbation | n | boundary moved by exactly the perturbation | S20 call | structural allele | native / within_family / foreign allele |
|---|---|---|---|---|---|
| +250 bp | 600 | 94.3% | 99.17% | 97.83% | 99.67 / 95.83 / 96.03% |
| +500 bp | 600 | 94.3% | 100% | 98.17% | 100 / 97.92 / 96.03% |
| -250 bp | 584 | 87.2% | 99.14% | 98.46% | 99.65 / 95.74 / 97.58% |
| -500 bp | 575 | 91.0% | 99.83% | 98.78% | 100 / 100 / 97.14% |

The residual changes have two causes:
* The caller's boundary moves by more or less than the inserted or trimmed amount (5.7-12.8% of ends).
* A fixed window length moves the inner end of the window: padding drops sequence at the far end, and trimming is N-padded.

±250 bp is not a multiple of the 100 bp grid, so it also changes the grid phase relative to the subtelomere. Native alleles
are essentially invariant (>= 99.65%).

## (e) mappy backend (`backend_mappy*.tsv`)

On 400 ends, mappy (all chains kept, `-P` mode; mappy cannot set `-p 0.01`) agrees with the minimap2 backend for:
* S5 86.5%, S10 89.8% and S20 93.0%;
* profile 60.3%, core structure 92.0% and structural allele 91.5%;
* haplogroup 100%.

It runs at 3.9 s per end single-threaded. On 80 of these ends, mappy's default chain filtering did worse (structural allele
77.5% vs 88.7%). Under revision 2 the mappy fallback is therefore a clearly weaker approximation than in revision 1. The CLI
prints a warning when it is used, and the minimap2 executable is required for catalogue-exact results.

## (f) Runtime (macOS arm64 machine shared with other jobs)

| input | ends | threads | wall time | peak memory |
|---|---|---|---|---|
| CHM13 v2.0, 3.1 Gb haploid (gzip FASTA) | 48 | 3 | 42.2 s | 2.23 GB |
| HG002 Q100 v1.1, 6.0 Gb diploid (gzip FASTA) | 86 | 3 | 77.6 s | 2.38 GB |
| pedigree ends, 368 × 250 kb | 366 | 3 | 207.7 s | 3.35 GB |
| windows mode, 5,981 catalogue windows | 5,981 | 4 | 0.478 s per end | |

The CHM13 and HG002 runs shared the machine with a 1-thread validation job. Per end, the tool runs one track alignment
against the 24 Mb CHM13 end set (up to 200 secondary chains) and the bin alignment used by the arm rule. A typical assembly with
about 90 telomeric ends takes 1-1.5 min.

## (g) Experimental read mode (`reads_hg002*.tsv`)

* **Reads.** Local slices of public HG002 alignments to HG002 Q100 v1.1 (Revio HiFi and ONT ultra-long). Up to 4 telomeric
  reads (>= 10 kb, >= 500 bp telomere) per end.
* **Truth.** The HG002 Q100 end the read is aligned to.

| reads | n reads (ends) | median window | arm | S20 | core structure | structural allele | haplogroup |
|---|---|---|---|---|---|---|---|
| HiFi, arm inferred | 354 (91) | 21.5 kb | 72.6% | 72.9% | 70.9% | 69.2% | 46.9% |
| HiFi, arm given | 354 (91) | 21.5 kb | 100% | 92.9% | 74.0% | 87.9% | 61.9% |
| ONT UL, arm inferred | 365 (92) | 100 kb | 81.4% | 81.1% | 81.1% | 78.9% | 28.5% |
| ONT UL, arm given | 365 (92) | 100 kb | 100% | 90.7% | 86.3% | 87.1% | 34.5% |

Read mode stays experimental. With a known arm, S20 calls and structural alleles are usable. Single-read haplogroups are not
reliable at ANI 0.995: sequencing errors lower k-mer ANI, and short reads cover only part of the 80 kb compared span.

## (h) Tests and clean install (`clean_install.tsv`)

* `python -m pytest tests/`: 23 passed and 1 skipped (the full-reference comparison, which passes when `SUBTEL_ATLAS_REFS`
  is set), in 15 s offline. The fixtures take 0.34 MB.
* A fresh conda env (python 3.11 + minimap2 2.31) took `pip install .`, pulling only mappy, numpy and pandas. `--help`,
  `fetch-refs` from the local checkout (39 s, all files verified) and a fixture query (2.2 s) all reproduced the catalogue calls.

## Failure modes

1. **Ambiguous structural alleles.** An end can lie within 10 kb of two same-structure allele medians. The query rule takes the
   nearest median, then the larger allele; the catalogue's clustering can place a defining end in the other allele (11/5,981).
   Flag: `sa_multiple_candidates`; read `sa_candidates`.
2. **Boundary sensitivity of non-native calls.** Moving the boundary by 500 bp changes 5-6% of foreign and within-family
   allele calls and fewer than 0.5% of native ones. Flag: `boundary_unstable` (81% sensitivity, no false alarms; 6 of the
   574 ends it calls stable changed under re-aligned shifts).
3. **Diffuse haplogroups.** Medoid assignment differs from UPGMA membership mostly in acroP, 10q and PAR1. Flags:
   `hg_close_second`, low `hg_ani`.
4. **Unrepresented rare haplogroups** (1-2 samples) have no medoid in `representatives.fa.gz`. Their ends are reported below
   threshold (170/189) or absorbed into a neighbour (19/189).
5. **mappy fallback.** About 92% agreement on structural alleles; use the minimap2 executable.
6. **Arm inference without names or probes** (mosaic rule): 84% correct at group level. Contig-name arms assume standard
   orientation.
7. **Telomeres < 500 bp** are not reported in contig mode, by definition; `--min-telomere` lowers the threshold.

## Deviations and scope

* Self-consistency used a stratified sample of 5,981 of the 20,284 non-pedigree ends (29%), not the whole catalogue.
* Leave-one-sample-out re-clustered structural alleles exactly. For haplogroups it re-selected medoids and representation
  without re-running UPGMA.
* The boundary-shift test used 600 ends (half with non-native alleles, i.e. non-native ends over-represented relative to the
  catalogue).
* The probe-rule table was carried over from the revision-1 run, because the arm method and arm calls did not change.
* The per-end `nativity_class` field was added to the tool after the self-consistency query had started. Its agreement
  (5,981/5,981) was computed from the recorded core structures with the same function.
* No HG002 Release 2 comparison (no such ends in v0.1). Read mode was tested on 719 HG002 reads from two technologies only.
* The test suite uses reduced references (6 CHM13 ends, telomeric 150 kb kept, rest N-masked). Its golden outputs test
  regression and determinism; agreement with the catalogue is covered here and by the optional full-reference test.
