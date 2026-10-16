"""subtel_atlas - query tool for the human subtelomere allele catalogue.

Assigns telomeric contig ends (or pre-cut telomere-first windows, or telomeric long reads) to the catalogue's arms, CHM13 mosaic
structures, structural alleles (SA-<arm>-NN) and sequence haplogroups (HG-<arm>-NN), using the same telomere caller, 5 kb bin
method, structure rules and haplogroup sequence span as the catalogue build (pipeline/catalogue_*.py imports these functions).
"""
__version__ = '0.1.0'
CATALOGUE_VERSION = 'v0.1'
