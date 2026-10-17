# Read-level evidence for non-native terminal blocks

Assembly-free test of single long reads against a catalogue claim that a chromosome end's telomere-adjacent block is closer
to another CHM13 v2.0 end (donor D) than to its own arm. A conclusive read:
* reaches the telomere;
* carries D alleles at paralog-diagnostic units distal to a maximum-likelihood switch;
* crosses the switch;
* continues into an anchor that aligns genome-wide with MAPQ 60 to the recipient chromosome.

**Start here:** [READ_EVIDENCE.md](READ_EVIDENCE.md) ([PDF](READ_EVIDENCE.pdf)) is the reviewer-facing report. The web
version is the *Read evidence* section of the site (`docs/evidence/`).

Reviewer resource (built by `scripts/build_evidence_set.py` and `scripts/evidence_figures.py` from the runs below):
* [claims.tsv](claims.tsv): one row per claim (HG002 Q100, HG003/HG004 trio tests, 1000 Genomes ONT). It gives the v1 verdict,
  the post-hoc **definitive** tier, read counts by tier and dataset, the model and the assembly switch with its CHM13 projection.
  The definitive tier is reporting only. It was defined after the population diagnostics and can only remove claims; the
  rule is in `scripts/evidence_lib.py`.
* [reads.tsv.gz](reads.tsv.gz): one row per cited read: public source URL, region, read group, units, purity, switch interval
  (read, subtelomere and CHM13 coordinates), anchor, own-assembly placement, tier, and the samtools command that fetches the
  read. [controls.tsv](controls.tsv): negative controls.
* [showcase/](showcase/README.md): eight exemplar claims with CRAM/BAM (CHM13 v2.0 and own assembly), window FASTA, read segments
  split at the switch, IGV sessions and "what you should see" notes. The full per-claim set is described in
  [release_manifest.tsv](release_manifest.tsv) (GitHub release assets, SHA-256).
* [replicate_claim.sh](replicate_claim.sh) `CLAIM_ID`: re-derives a verdict from public sources only.
  [REPLICATION_TEST.md](REPLICATION_TEST.md) records a clean-tree test ([replication_results.tsv](replication_results.tsv)).
* [figures/](figures/): read-through figures of the exemplars, the summary figure, and small multiples of every definitive claim.
* [SHA256SUMS](SHA256SUMS): checksums of the committed evidence files (checked by `verify_release.sh`).

Pipeline records:
* [CRITERIA.md](CRITERIA.md): pre-registered criteria, with calibration amendments A1 and C1-C3.
* [criteria_v1.json](criteria_v1.json): frozen thresholds (sha `a89c5eaf023ed032`). The registered version is
  [criteria_v1_preregistered.json](criteria_v1_preregistered.json).
* [calibration/CALIBRATION.md](calibration/CALIBRATION.md): sensitivity and specificity on HG002 Q100 (training and held-out
  halves) and on strictly native controls (HG002 and 1000 Genomes ONT).
* [hg002/README.md](hg002/README.md): all 24 HG002 Q100 claims re-tested on every public Q100 long-read alignment (16 files,
  12 independent datasets after molecule de-duplication), plus the trio transmission test in HG003/HG004 (GIAB).
* `scripts/`:
  * `hg002_multiplatform.py`: HG002 multi-platform spanning reads (fetch, evaluate, report, evidence BAMs);
  * `trio_reads.py`: HG003/HG004 transmission test;
  * `fetch_reads.py`: region fetch from the public alignments;
  * `run_calibration.py`: calibration runs;
  * `site_error_rates.py`: per-read error at diagnostic sites;
  * `calibration_report.py`: tables and figures;
  * `run_population.py`, `population_report.py`, `population_diagnostics.py`: 1000 Genomes ONT claims and diagnostics;
  * `build_evidence_set.py`: unified tables, own-assembly placement, CRAM/BAM set, showcase, release bundles, IGV sessions;
  * `evidence_figures.py`: read-through, summary and per-claim figures;
  * `replicate_claim.py`: the replication entry point behind `replicate_claim.sh`;
  * `evidence_lib.py`: shared helpers (placement, definitive tier, public-source commands).
* Command line: `subtel-atlas reads-evidence` (module `atlas/subtel_atlas/readevidence.py`; tests in
  `tests/test_readevidence.py`).
