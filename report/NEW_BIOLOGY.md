# Common terminal alleles of human chromosome ends
*What the human subtelomere atlas adds: a validated allele catalogue of HPRC Release 2 chromosome ends, its replication in independent cohorts, and tests of mechanism, age and telomere length · October 2026 · atlas v0.1*

## Summary

- **Resource.** 20,144 telomeric contig ends from 230 HPRC Release 2 individuals were described as mosaics of CHM13 chromosome ends and grouped into 647 structural alleles (tiers A/B/C 129/166/352) and 593 sequence haplogroups, each with evidence and a query tool.
- **Foreign ends are alleles.** 24.5% of ends carry a block resembling an unrelated arm. Common forms are recurrent, old and shared across continents. SA-20q-01, a 3q-like terminal block on 20q, is carried by 146 individuals (frequency 0.48).
- **Replicated.** In 11,175 sequences from four independent cohorts (HGSVC3, CPC, APC, JASA; Tan et al. 2026), every common tier-A allele recurs in the pooled cohorts, and non-native allele frequencies correlate with Release 2 (r = 0.9708).
- **Mechanism.** Transfer is one-way and homology-mediated. Donor ends of carriers look like those of non-carriers (reciprocal 3.86% vs 3.91%), and 327 of 383 switch clusters lie within donor-recipient homology.
- **Telomere length.** Arm-specific lengths reproduce (Spearman ρ = 0.73 against Karimian et al. 2024). The terminal allele has no detectable cis effect once telomere-boundary artefacts are controlled.
- **Novelty.** Multi-haplotype arms and shared blocks confirm Tan et al. 2026 and Guarracino et al. 2026. The validated allele-level catalogue, cross-cohort replication, one-way homology-mediated mechanism, pre-out-of-Africa ages and the negative cis result are new.

## 1. The atlas

The atlas describes the outer 100 kb of every telomeric contig end in HPRC Release 2 (455 haplotype assemblies, 20,144 ends, including 8,146 on contigs without a chromosome name, whose arm was inferred). Windows start at the telomere/subtelomere boundary. One minimap2 pass against the 48 CHM13 v2.0 ends gives per-base identity tracks, and each 100 bp unit is labelled with the best-matching end. Near-identical CHM13 ends form six paralog families (PF01 = 4p/14p/15p/21p/22p; PF02 = 1p/5q/6q; PF03 = 11p/15q; PF04 = 3q/19p; PF05 = Xp/Yp; PF06 = Xq/Yq). A block is *home* if the own arm is within 0.5 identity points of the best match. Ends then fall into three nativity classes: native (14,461), within-family (747, closer to another member of a family containing the own arm) and foreign (4,936).

A **structural allele** is a recurring block order with switch positions within 10 kb. A **haplogroup** is a cluster at ANI ≥ 0.995 over the first 80 kb. Evidence tiers combine recurrence in unrelated samples with orthogonal lines of support: Verkko concordance, reads across switch points, intact pedigree transmission and presence in a T2T assembly. Labels are relative to CHM13. A switch marks a change in which CHM13 end the sequence most resembles, usually inside shared paralogous sequence, and need not be a structural breakpoint.

![F1_genome_overview.png](figures/atlas/F1_genome_overview.png)

**Figure 1. Terminal blocks shared between arms.** (a) Chords join the best-matching CHM13 end of each terminal non-home block to the recipient arm (tier A/B alleles), with width proportional to frequency on the recipient arm. The outer bars show the number of structural alleles per arm. (b) Nativity class of Release 2 ends per arm.

![F2_allele_mosaics.png](figures/atlas/F2_allele_mosaics.png)

**Figure 2. Allele mosaics of six arms.** Each row is one Release 2 end, with the telomere at left. Blocks are coloured by their best-matching CHM13 end, and rows are grouped by structural allele (frequency and tier at right).


## 2. Validation

| Check | n | Result |
|---|---|---|
| Same haplotype, independent assembly (Tan et al. Verkko vs our hifiasm) | 10,999 | structural allele 96.4%, form 98.7%, S20 97.8%; base-level switches identical after projection 95.9% |
| By nativity class (same comparison) |  | foreign 90.9%; native 99.6%; within_family 83.4% |
| Boundary shifted ±500 bp and re-aligned | 1043 | non-native allele unchanged 96.2% / 94.6%; native 97.7% / 99.0% |
| Query tool re-querying catalogue ends | 5,981 | mosaic, S-states and form 100%; structural allele 99.8% |
| Leave one sample out (60 samples) | 5,293 | structural allele 99.8%; held-out-only alleles reported novel 74/76 |
| Independent T2T / pedigree assemblies | CHM13, HG002 Q100, WashU | structural allele 100% in each (tests/validation/independent_summary.tsv) |

**Design change before release.** A first design tiled 5 kb bins from the telomere boundary. It was rejected because only 40% of non-native allele calls survived a 500 bp boundary shift. The track-and-family design used here raised agreement with independent assemblies from 93.6% to 96.4% (data/catalogue/QC.md, CHANGELOG). Calls 5 kb from the boundary (S5) remain boundary-sensitive and are not recommended as primary per-end calls.

![tan2026_concordance.png](figures/validation/tan2026_concordance.png)

**Figure 3. Cross-check against Tan et al. 2026.** (a) Agreement on the same haplotype. (b) Switch positions after alignment projection. (c, d) Allele frequencies in four independent cohorts compared with HPRC Release 2.


## 3. Common terminal alleles are old and replicate

All 88 common tier-A alleles (Release 2 frequency ≥ 5%) recur in the pooled independent cohorts, each in at least three individuals (11,175 sequences from 182 individuals). Per cohort, HGSVC3 88/88, CPC 78/88, APC 87/88, JASA 76/86 are observed. Of 5,573 allele-frequency tests, 14 reach q < 0.05 against the all-ancestry reference, and none against a matched superpopulation. In Release 2, 41 of 104 common forms differ among superpopulations (median FST 0.023), and 66 occur in all five. 34 of 39 switch clusters seen in ≥ 10 individuals are shared between continental African and non-African samples. Carriers of 23 of 34 testable forms cluster tightly (median nearest-neighbour carrier fraction 0.95 vs 0.14 expected), but only 13 of 34 meet the full single-origin criteria; 7 clustered forms are as diverse as native alleles, and 11 are not supported. None is strictly monophyletic, which is compatible with gene conversion at block edges. Median-based clock ages: 20q>3q 623 kyr, 3q>11p 527 kyr, 6p>11p 457 kyr. These ages assume a genome-average mutation rate and are order-of-magnitude estimates. No recurrent human switch configuration was found in five T2T great-ape genomes.

![fig_origin.png](analyses/allele_ages/figures/fig_origin.png)

**Figure 4. Single origin and age.** (a) Diversity of carriers compared with native donor-arm alleles. (b) Distance of carriers, natives and other recipients to the carrier consensus. (c) Clock ages. (d) Great-ape screen.


## 4. Transfer is one-way and homology-mediated

For 104 forms (3,917 carrier ends), the donor arm in the same haplotype is assessable for 2,464 carriers. Their donor ends look like those of non-carriers: reciprocal products 3.86% vs 3.91%, and native 47.0%. So the common terminal alleles arose by copying the donor terminal sequence, not by reciprocal exchange. At the 383 base-resolved switch clusters, 327 lie within donor-recipient homology, and 261 share more than 25 bp of identical parental sequence. Fewer than 1e-05 of random position pairs do. These are signatures of homology-directed copying (break-induced replication or gene conversion). Blunt or microhomology end-joining does not fit. Interstitial telomeric repeats are enriched only within 5 kb of the telomere.

![junctions_mechanism.png](analyses/junctions/figures/junctions_mechanism.png)

**Figure 5. Mechanism.** (a) Donor-end states of carriers and non-carriers. (b) Length of replaced recipient sequence. (c) Identical parental sequence at the switch. (d) Interstitial telomeric repeats near switches. (e) Three base-level exemplars, two in HG002 with read support.


## 5. No cis effect on telomere length

Assembled per-arm telomere length correlates with the long-read per-arm values of Karimian et al. 2024 (ρ = 0.73, 38 arms). Within individuals heterozygous for a foreign and a native allele of the same arm, an earlier analysis window suggested foreign alleles carry longer telomeres (+408 bp). That window shrinks as the telomere grows, so the association was an artefact. Windows anchored at the boundary give small negative differences (S20 -211 bp; S5 -515 bp). These arise because the telomere caller stops early at degenerate repeats next to the telomere. With a boundary-robust telomere length, the paired differences are -133 to +14 bp (all p ≥ 0.36). The minimum detectable effect is about 408 bp. Tan et al. 2026 report suggestive haplotype associations for 1p, 5q and 6q (FDR < 0.1); this within-person design neither supports nor excludes them.

![fig2_allele_test.png](analyses/telomere_cis/figures/fig2_allele_test.png)

**Figure 6. Telomere length and terminal allele.** (a) Paired differences by allele window. (b) Per-form differences. (c, d) Donor-likeness. (e) The window artefact.


## 6. Comparison with the literature

| Claim | Status before this work | Atlas evidence |
|---|---|---|
| Chromosome ends are mosaics of shared blocks; arms carry several haplotypes | Established (Linardopoulou 2005; Tan 2026; Guarracino 2026) | Confirmed; not claimed as new |
| Common foreign terminal blocks are real alleles, not assembly errors | Assumed in recent preprints, not validated per allele | Validated: second assembler, reads, pedigree, independent cohorts |
| Terminal replacements form by one-way, homology-directed copying | Contested (NHEJ inferred by Linardopoulou 2005) | Supported for common alleles |
| Terminal alleles are few, old and of single origin | Established only for 4q/10q (Lemmers 2010) | Carriers cluster for 23/34 testable forms; full single-origin criteria met by 13/34; switch clusters shared across Africa and non-Africa |
| Telomere length is set in cis by the terminal allele | Suggestive associations (Tan 2026); untested within persons | No detectable effect (≥ ~0.4 kb excluded) |
| Ectopic subtelomere exchange is frequent in pedigrees | Claimed (Guarracino 2026) | Not supported: 51 of 54 foreign ends inherited intact (3 transmissions) |

The full screened review (199 records) and the novelty table are in literature/REVIEW.md and literature/NOVELTY.md.


## 7. Limitations

Labels are relative to a single reference (CHM13), so an arm whose CHM13 copy is a minority haplotype appears foreign in most people. Acrocentric p arms and the pseudoautosomal regions are pooled. Duplicon-rich arms carry many rare structures with little orthogonal support (tier C). Most source DNA comes from cell lines, and assembled telomere lengths are lower bounds. Clock ages ignore gene conversion and assume a genome-average mutation rate. The same-haplotype cross-check uses the same sequencing data and tests assembly robustness, not sequencing error. Read-level assignment in the query tool is experimental.


## Data and code

Catalogue, figures, site and code: github.com/jlanej/human-subtelomere-atlas (data CC BY 4.0, code MIT). Browse at jlanej.github.io/human-subtelomere-atlas. Reproduction commands are in the repository README.


## References

- Guarracino et al. (2026). Concerted evolution and unorthodox recombination of human subtelomeres. *bioRxiv (Cold Spring Harbor Laboratory)*. doi:10.64898/2026.07.10.737660
- Tan et al. (2026). Widespread structural variations at human chromosome ends. *bioRxiv (Cold Spring Harbor Laboratory)*. doi:10.64898/2026.09.03.748401
- Karimian et al. (2024). Human telomere length is chromosome end–specific and conserved across individuals. *Science*. doi:10.1126/science.ado0431
- Linardopoulou et al. (2005). Human subtelomeres are hot spots of interchromosomal recombination and segmental duplication. *Nature*. doi:10.1038/nature04029
- Lemmers et al. (2010). Worldwide Population Analysis of the 4q and 10q Subtelomeres Identifies Only Four Discrete Interchromosomal Sequence Transfers in Human Evolution. *The American Journal of Human Genetics*. doi:10.1016/j.ajhg.2010.01.035
- Nurk et al. (2022). The complete sequence of a human genome. *Science*. doi:10.1126/science.abj6987
- Yoo et al. (2025). Complete sequencing of ape genomes. *Nature*. doi:10.1038/s41586-025-08816-3
- Tan et al. (2026) subtelomere sequences, Zenodo, CC BY 4.0. doi:10.5281/zenodo.21807156
- Human Pangenome Reference Consortium, Release 2 (BioProject PRJNA730823).
