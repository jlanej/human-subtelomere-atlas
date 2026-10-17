# Regenerating local intermediates (`work/`)

Everything under `work/` is git-ignored scratch: downloads, minimap2 indexes, alignments and per-run caches. None of it is
needed to use the atlas. The committed tables, figures, website and evidence showcase are complete without it. The 2026-10
release was built with a `work/` of about 77 GB, and that copy has been removed. Each component rebuilds into `work/` from
public data with the commands below, run from the repository root.

Set machine-specific paths first. `work/config.local.sh` is git-ignored and is read by the runners; it sets `ATLAS_THREADS`,
and `ATLAS_INPUT_DIR` for builds from an analysis home. Builds from public data only need `pipeline/INPUTS.md`.

| `work/` subdirectory (size at release) | What it held | Regenerate with |
|---|---|---|
| `ref/`, `catalogue/*.mmi` (3.7 GB, 5.5 GB) | CHM13 v2.0 and GRCh38 FASTA, minimap2 indexes | `bash pipeline/fetch_public_inputs.sh chm13 pedigree`; the `.mmi` files were optional speed-ups (`minimap2 -x asm5 -d` / `-x map-ont -d`; the evidence scripts take one with `--genome`) |
| `windows/` (3.5 GB) | 500 kb end windows of the 455 Release 2 haplotypes, HG002 Q100 and pedigree | `pipeline/catalogue_01_plan_windows.py`, `catalogue_02_fetch_windows.sh`, `catalogue_03_assemble_windows.py build` |
| `catalogue/` (11 GB) | bin and track alignments (PAF), mosaics, stability shifts | `bash pipeline/run_catalogue.sh` (steps in `pipeline/CATALOGUE.md`) |
| `tan2026/` (16 GB) | Tan et al. 2026 Zenodo sequences (10.5281/zenodo.21807156), their track alignments | `bash pipeline/validation/run_tan2026.sh` |
| `qtool/` (2.1 GB) | query-tool validation runs (self-consistency, leave-one-sample-out, HG002) | `python tests/validation/run_validation.py` (subcommands in its docstring) |
| `telomere_cis/`, `allele_ages/`, `junctions/` (0.9, 0.9, 10 GB) | per-analysis downloads (HG002 Q100 genome, ape ends, native 100 kb windows) | `bash analyses/<name>/run.sh` |
| `readevidence/` (8.4 GB) | calibration and population read slices, CHM13 map-ont index | `evidence/scripts/fetch_reads.py`, `run_calibration.py`, `site_error_rates.py`, `calibration_report.py`, `run_population.py all` then `report` (`evidence/calibration/CALIBRATION.md`, `evidence/population/README.md`) |
| `evidence/hg002/`, `evidence/trio/`, `evidence/population/` (2.9, 1.4, 5.4 GB) | HG002 Q100 multi-platform, HG003/HG004 and 1000 Genomes ONT read slices | `evidence/scripts/hg002_multiplatform.py`, `trio_reads.py fetch/eval/report` (`evidence/hg002/README.md`) |
| `evidence/release_build/`, `evidence/release/` | per-claim CRAM/BAM directories and the four release tarballs | `evidence/scripts/build_evidence_set.py` (`release` writes the tarballs, `SHA256SUMS` and `evidence/release_manifest.tsv`) |
| `literature/` (50 MB) | Crossref, OpenAlex and PubMed caches | `bash literature/reproduce.sh` (`literature/README.md`, section "Reproduce"; needs an OpenAlex API key) |
| `verify_release/`, `replicate_cache/`, `site_check/` | outputs of the checks | `bash verify_release.sh [--network]`, `python website/check_site.py` |
| `scratch_public/`, `scratch07/`, `test/`, `proto/`, `tmp*/` | clean-clone rebuild tests and prototypes | none needed (`pipeline/verify_public_rebuild.py` compares a subset rebuild with the released `ends.tsv`) |

Files that are too large for git but belong to the release are GitHub release assets, listed with SHA-256 in
[evidence/release_manifest.tsv](evidence/release_manifest.tsv). The tables asset also carries the full per-read population
table (`population_reads_all.tsv.gz`, 77,333 rows), which `evidence/scripts/population_report.py` writes to
`evidence/population/reads.tsv.gz`.

`verify_release.sh` needs none of `work/`. Without `--network` it checks the committed files. With `--network` it re-fetches
30 windows and re-derives one read-evidence claim from public reads. The claim check needs CHM13 v2.0: it takes it from
`SUBTEL_CHM13` or downloads it into `work/replicate_cache`.
