# Data licence and attribution

## Atlas data, figures and documentation

The catalogue tables, derived sequences, figures, website content and written documentation in this
repository (`data/`, `figures/`, `docs/`, `analyses/*/results/`, `analyses/*/figures/`, `literature/`
and Markdown files) are released under the
[Creative Commons Attribution 4.0 International licence (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
Please cite the atlas as described in [CITATION.cff](CITATION.cff) when you use these data.

Source code is released under the MIT licence ([LICENSE](LICENSE)).

## Upstream data

The atlas is derived from public data. Upstream terms continue to apply to the source data.

- **Human Pangenome Reference Consortium (HPRC) Release 2.** Assemblies and reads are in the public
  domain (https://humanpangenome.org/data-use/). Users are asked to follow the HPRC Data Use Best
  Practices and to include the following statement: *"We would like to acknowledge the Human
  Pangenome Reference Consortium (BioProject ID: PRJNA730823) and its funder, the National Human
  Genome Research Institute (NHGRI)."*
- **T2T-CHM13 v2.0** (doi:10.1126/science.abj6987) and **HG002 Q100 v1.1** (Genome in a Bottle /
  T2T consortium): public reference assemblies.
- **WashU pedigree (PAN010, PAN011, PAN027, PAN028)** assemblies (doi:10.64898/2025.12.14.693655):
  public assemblies, used here for validation only.
- **1000 Genomes Project / IGSR** sample and population tables (doi:10.1016/j.cell.2022.08.004;
  doi:10.1038/nature15393).
- **Tan et al. 2026**, *Widespread structural variations at human chromosome ends*, bioRxiv
  doi:10.64898/2026.09.03.748401. Subtelomere sequences from Zenodo doi:10.5281/zenodo.21807156 are
  licensed CC BY 4.0. Tables in `data/validation/` are derived from them; the sequences themselves
  are not redistributed here.
- **T2T great-ape assemblies** (Yoo et al. 2025, doi:10.1038/s41586-025-08816-3), used for the
  outgroup screen in `analyses/allele_ages`.

## Responsible use

Population labels (superpopulation, population) are those of the 1000 Genomes Project and HPRC and
describe sampling, not biological categories. Frequencies by superpopulation are given as descriptive
statistics. They should not be used to make claims about individuals or to stigmatise groups.
Following HPRC best practice, do not attempt to re-identify participants.
