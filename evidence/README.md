# Read-level evidence for non-native terminal blocks

Assembly-free test of single long reads against a catalogue claim that a chromosome end's telomere-adjacent block is closer
to another CHM13 v2.0 end (donor D) than to its own arm. A conclusive read:
* reaches the telomere;
* carries D alleles at paralog-diagnostic units distal to a maximum-likelihood switch;
* crosses the switch;
* continues into an anchor that aligns genome-wide with MAPQ 60 to the recipient chromosome.

* [CRITERIA.md](CRITERIA.md): pre-registered criteria, with calibration amendments A1 and C1-C3.
* [criteria_v1.json](criteria_v1.json): frozen thresholds (sha `a89c5eaf023ed032`). The registered version is
  [criteria_v1_preregistered.json](criteria_v1_preregistered.json).
* [calibration/CALIBRATION.md](calibration/CALIBRATION.md): sensitivity and specificity on HG002 Q100 (training and held-out
  halves) and on strictly native controls (HG002 and 1000 Genomes ONT).
* `scripts/`:
  * `fetch_reads.py`: region fetch from the public alignments;
  * `run_calibration.py`: calibration runs;
  * `site_error_rates.py`: per-read error at diagnostic sites;
  * `calibration_report.py`: tables and figures.
* Command line: `subtel-atlas reads-evidence` (module `atlas/subtel_atlas/readevidence.py`; tests in
  `tests/test_readevidence.py`).
