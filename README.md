# human-subtelomere-atlas

A validated, queryable catalogue of the terminal 100 kb of human chromosome ends (subtelomeres), built from the
HPRC Release 2 assemblies (230 individuals, 455 haplotype assemblies, 20,144 telomeric contig ends), with T2T-CHM13 v2.0
and HG002 Q100 v1.1 as references and WashU pedigree assemblies used for validation.

- **Browse:** https://jlanej.github.io/human-subtelomere-atlas/ (static site in [`docs/`](docs/))
- **Download:** [`data/catalogue/`](data/catalogue/) (tables, JSON, representative sequences; [SCHEMA](data/catalogue/SCHEMA.md), [QC](data/catalogue/QC.md))
- **Query your own assembly:** [`atlas/`](atlas/README.md) (`subtel-atlas query assembly.fa.gz`)

Release **v0.1** (pre-publication). Allele identifiers are stable within v0.1; see `QC.md` for the change log.

## What the catalogue contains

Each chromosome end is described as a **mosaic**: the terminal 100 kb, measured from the telomere boundary, is labelled
block by block with the CHM13 chromosome end it most resembles. Ends are grouped in two complementary ways:

| Unit | Definition | Count (tiers A / B / C) |
|---|---|---|
| Structural allele (`SA-<arm>-NN`) | Recurring block order and switch positions (10 kb tolerance) | 647 (129 / 166 / 352) |
| Haplogroup (`HG-<arm>-NN`) | Sequence cluster, ANI ≥ 0.995 over the first 80 kb | 593 (131 / 152 / 310) |

Near-identical CHM13 ends are grouped into six **paralog families** (PF01 4p + acrocentric p arms; PF02 1p/5q/6q;
PF03 11p/15q; PF04 3q/19p; PF05 Xp/Yp; PF06 Xq/Yq). Each end has a **nativity class**: of 20,144 Release 2 ends,
14,461 are native (own-arm sequence throughout), 747 within-family and 4,936 foreign (a terminal or internal block that
best matches an unrelated arm). Evidence tiers combine recurrence across unrelated samples with orthogonal support:
second-assembler (Verkko) concordance, read support across switch points, intact pedigree transmission, and presence in
a T2T assembly. Labels are relative to CHM13; a "switch" is a change in which CHM13 end the sequence most resembles,
usually inside shared paralogous sequence, and is not necessarily a structural breakpoint.

## Main findings

1. **Foreign terminal segments are common, old, recurrent alleles.** For example, SA-20q-01 (a 3q-like terminal block
   on 20q, median switch 31.5 kb from the telomere boundary) is carried by 146 individuals (Release 2 frequency 0.48).
   Shared switch points recur across unrelated individuals and across continental groups. Median-based clock ages of
   clustered forms are 0.2-1.5 Myr; these ages are model-dependent and order-of-magnitude only.
2. **They replicate in four independent cohorts.** Classifying 11,175 subtelomere sequences from HGSVC3, CPC, APC and JASA
   (Tan et al. 2026, CC BY 4.0) with the same method, all 88 common tier-A alleles recur in the pooled cohorts
   (76-88 per cohort), and non-native allele frequencies correlate with Release 2 at r = 0.97. No allele differs
   significantly from a matched 1000 Genomes superpopulation.
3. **Transfer is one-way and homology-mediated.** Carriers' donor ends look like non-carriers' (reciprocal products
   3.86% vs 3.91%, p = 0.88), and 85% of base-resolved switches lie within donor-recipient homology
   ([analyses/junctions](analyses/junctions/README.md)).
4. **No cis effect on telomere length.** Arm-specific telomere lengths reproduce published values (Spearman ρ = 0.73
   against Karimian et al. 2024), but foreign and native alleles of the same arm do not differ in telomere length once
   telomere-boundary artefacts are controlled ([analyses/telomere_cis](analyses/telomere_cis/README.md)).

Details, statistics and limitations: [analyses/](analyses/), the [Findings pages](https://jlanej.github.io/human-subtelomere-atlas/findings/),
and the literature comparison in [literature/NOVELTY.md](literature/NOVELTY.md).

## Relation to concurrent work

Two 2026 preprints analyse the same assemblies.
- Tan et al., *Widespread structural variations at human chromosome ends*, bioRxiv doi:10.64898/2026.09.03.748401:
  860 haplotypes, paralogy-block haplotypes per arm.
- Guarracino et al., *Concerted evolution and unorthodox recombination of human subtelomeres*, bioRxiv
  doi:10.64898/2026.07.10.737660.

This atlas is complementary. It adds an allele-level, evidence-tiered catalogue of terminal structures with shared
switch points, a query tool, a cross-check against the Tan et al. sequences, and tests of mechanism, age and cis
telomere-length effects. Multi-haplotype arms and the prevalence of shared blocks confirm those studies; they are not
claimed as new.

## Validation (summary)

| Check | Result | Source |
|---|---|---|
| Same haplotype, independent assembler (Tan et al. Verkko vs our hifiasm; 10,999 pairs) | structural allele 96.4%, form 98.7% | [data/validation](data/validation/README.md) |
| Robustness to ±500 bp shifts of the telomere boundary (1,043 ends) | non-native alleles 96.2% / 94.6% unchanged | [QC.md](data/catalogue/QC.md) |
| Query tool vs catalogue (5,981 re-queried ends; 60-sample leave-one-out) | structural allele 99.8% / 99.8% | [tests/validation](tests/validation/VALIDATION.md) |
| Independent T2T and pedigree assemblies (CHM13, HG002 Q100, WashU pedigree) | structural allele 100% | [tests/validation](tests/validation/VALIDATION.md) |

## Read-level evidence

Single long reads not used for the assemblies test the non-native terminal blocks directly. A read must start in the
telomere, carry the donor arm's paralog-diagnostic alleles up to the switch and the recipient's beyond it, and anchor
uniquely in the recipient chromosome; the criteria were frozen before population reads were inspected.

* HG002 Q100: all 15 v1-confirmed claims are *definitive* (12 with definitive reads from >= 2 independent datasets).
* Trio: the transmitting parent (HG003/HG004) confirms 4 of the HG002 claims.
* 1000 Genomes ONT: 24 of the 43 v1-confirmed claims in 42 individuals are definitive.
* Controls: 0 of 297 HG002 and 0 of 1,010 1KG-ONT negative-control tests were conclusive.

The definitive tier was defined after the population diagnostics and only removes claims. Report:
[evidence/READ_EVIDENCE.md](evidence/READ_EVIDENCE.md) ([PDF](evidence/READ_EVIDENCE.pdf)). Tables:
[claims.tsv](evidence/claims.tsv), [reads.tsv.gz](evidence/reads.tsv.gz). Exemplar reads with IGV sessions:
[evidence/showcase](evidence/showcase/README.md). Full CRAM/BAM set: release assets listed in
[evidence/release_manifest.tsv](evidence/release_manifest.tsv). Re-derive any verdict from public data with
`evidence/replicate_claim.sh CLAIM_ID` ([clean-tree test](evidence/REPLICATION_TEST.md): 9 of 9 runs identical).

## Quick start: query tool

```bash
conda create -n subtel -c conda-forge -c bioconda python=3.11 pip minimap2=2.31 && conda activate subtel
pip install .                  # from a checkout of this repository
subtel-atlas fetch-refs        # CHM13 v2.0 ends (range reads) + catalogue files, SHA-256 checked
subtel-atlas query my_assembly.fa.gz -t 4 > my_ends.tsv
```

minimap2 is required for results that match the catalogue exactly; the `mappy` fallback is approximate. Read mode
(`subtel-atlas reads`) is experimental. See [atlas/README.md](atlas/README.md) for output columns and interpretation.

## Reproducing the release

- `bash verify_release.sh [--network]`: verifies all release checksums, catalogue integrity, the test suite and the
  README headline counts and the read-evidence checksums and criteria hash; `--network` also re-fetches 30 random windows from
  the public assemblies and checks their SHA-256, and re-derives one read-evidence claim from public reads.
- `python -m pytest tests/`: offline tests (about 15 s).
- [`data/inputs/`](data/inputs/): small derived inputs (1.8 MB) that make the catalogue rebuildable from public assemblies;
  the tip-discovery scan and the earlier evidence analyses enter as these snapshots (see `pipeline/INPUTS.md`).
- [`report/NEW_BIOLOGY.pdf`](report/NEW_BIOLOGY.pdf): summary report of the findings (`report/build_report.py`).
- `pipeline/run_catalogue.sh`: catalogue build (see [pipeline/CATALOGUE.md](pipeline/CATALOGUE.md) and
  [pipeline/INPUTS.md](pipeline/INPUTS.md) for inputs and the rebuild-from-public-data recipe).
- `analyses/*/run.sh`, `pipeline/validation/run_tan2026.sh`, `literature/reproduce.sh`: analyses, cross-check and
  literature review.
- `website/build_site.py`: regenerates `docs/` from the data files.

## Repository layout

```
data/catalogue/    release tables, catalogue.json, representatives.fa.gz, SCHEMA, QC, checksums
data/validation/   cross-check against Tan et al. 2026 (derived tables only)
data/metadata/     sample populations (1000 Genomes / HPRC)
pipeline/          catalogue build, figures, validation pipelines
atlas/             query tool (Python package subtel_atlas, CLI subtel-atlas)
analyses/          telomere_cis, junctions, allele_ages (each with README, scripts, results, figures)
evidence/          read-level evidence: criteria, calibration, HG002/trio/1KG-ONT results, showcase reads, replication
literature/        systematic review, novelty table, references
figures/           atlas and validation figures (PNG, PDF, SVG)
docs/              GitHub Pages site (generated by website/build_site.py)
tests/             tests and query-tool validation
```

## Limitations

- Labels are relative to CHM13 v2.0. Arms whose CHM13 copy is itself a minor haplotype appear "foreign" in most people.
- Calls at 5 kb from the boundary (S5) are sensitive to the telomere-boundary call; S20 and structural alleles are the
  recommended units. The per-end `boundary_stable` flags mark sensitive ends.
- Acrocentric p arms and the pseudoautosomal regions are pooled (acroP, PAR1, PAR2). Duplicon-rich arms carry many rare
  structures, and haplogroups seen in fewer than 3 samples have no representative sequence.
- Assembled telomere lengths are lower bounds, and the source DNA is mostly from cell lines.
- Clock ages assume a genome-average mutation rate and ignore gene conversion.

## Licence, citation and acknowledgements

Code: MIT ([LICENSE](LICENSE)). Data, figures and text: CC BY 4.0 ([DATA_LICENSE.md](DATA_LICENSE.md), which also
lists upstream terms). Cite via [CITATION.cff](CITATION.cff).

We would like to acknowledge the Human Pangenome Reference Consortium (BioProject ID: PRJNA730823) and its funder, the
National Human Genome Research Institute (NHGRI). We thank the T2T Consortium, the Genome in a Bottle Consortium, the
1000 Genomes Project, the WashU pedigree project, and Tan et al. for releasing their subtelomere sequences.
