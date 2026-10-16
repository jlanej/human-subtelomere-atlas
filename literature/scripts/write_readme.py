#!/usr/bin/env python3
"""Write literature/README.md from screening_flow.json, search_log.json and records.tsv."""
import csv, json, os
LIT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "literature")
flow = json.load(open(os.path.join(LIT, "screening_flow.json"))); log = json.load(open(os.path.join(LIT, "search_log.json")))
recs = list(csv.DictReader(open(os.path.join(LIT, "records.tsv")), delimiter="\t"))
npre = sum(r["preprint"] == "1" for r in recs)
txt = f"""# Literature review track

Systematic review supporting the human subtelomere atlas. Topics: subtelomeric duplication and inter-chromosomal exchange, telomere capture and break-induced replication, interstitial telomeric sequences, chromosome-end- and allele-specific telomere length, the 4q/10q, PAR1 and acrocentric exchange systems, complete genomes and pangenomes, existing subtelomere resources, long-read telomere methods and concerted-evolution theory.

## Deliverables

| File | Content |
|---|---|
| [REVIEW.md](REVIEW.md) | Structured narrative review by topic, citations as DOI links |
| [NOVELTY.md](NOVELTY.md) | Claim-by-claim table: status in the literature (established / contested / unknown), key citations, and how the final atlas results bear on each claim |
| [records.tsv](records.tsv) | {flow['included']} included records: DOI, title, first author, year, venue, OpenAlex ID, cited-by count, preprint flag, topic tags, relevance (1-3), one-line annotation, BibTeX key |
| [references.bib](references.bib) | BibTeX for all included records (Crossref metadata) |
| [doi_verification.tsv](doi_verification.tsv) | Crossref / doi.org resolution, Crossref update notices (corrections, retractions) per DOI |
| [curation.tsv](curation.tsv), [curation_excluded.tsv](curation_excluded.tsv) | Manual inclusion decisions with tags, relevance and annotation; exclusions after abstract reading with reasons |
| [screening.tsv](screening.tsv), [screening_flow.json](screening_flow.json) | Keyword pre-screen score and decision for every candidate; flow counts |
| [candidates.tsv](candidates.tsv), [seed_resolution.tsv](seed_resolution.tsv), [search_log.json](search_log.json) | All retrieved candidates with retrieval provenance (query, seed or citation step); seed DOI resolution; retrieval dates and endpoints |
| [evidence_numbers.json](evidence_numbers.json) | Atlas values quoted in NOVELTY.md and REVIEW.md, each with its source file; SHA-256 and modification time of every source |
| [cited_dois.txt](cited_dois.txt) | DOIs cited in the rendered documents |
| [SHA256SUMS](SHA256SUMS) | Checksums of the released files |

Topic tags in records.tsv: subtel_sd, subtel_exchange, telo_capture, bir, its, telo_length_cis, telo_methods, tvr, d4z4, par, acro, t2t, hprc, pedigree, concerted, resource, primate, yeast_model.

## Method

1. Retrieval (`scripts/fetch_literature.py`, configuration `scripts/search_config.json`): {log['n_queries']} OpenAlex relevance-ranked keyword queries (top {log['per_query']} each), {log['n_seeds']} seed papers pinned by DOI (resolved initially by Crossref free-text search and checked by title), and one citation-graph step (backward references and forward citations, most-cited and most-recent) from {log['n_expanded']} key papers. Raw responses were retrieved between {log['retrieved_utc_first']} and {log['retrieved_utc_last']} (UTC). {flow['candidates']} unique works.
2. Pre-screen (`scripts/screen_candidates.py`): deterministic title/abstract keyword score; seeds always retained; non-seeds retained at score 6 or more without off-topic organism flags. {flow['prescreen_retained']} retained, {flow['prescreen_excluded']} excluded.
3. Manual screening: titles of retained works read ({flow['excluded_at_title_screen']} excluded at title, {flow['retained_without_doi']} without DOI dropped); abstracts read for the shortlist (OpenAlex abstracts, PubMed abstracts via `scripts/fetch_pubmed_abstracts.py` where OpenAlex had none; the full text of the closest competitor preprint, doi:10.64898/2026.09.03.748401, was read for its telomere-length results); {flow['excluded_after_abstract']} excluded after abstract reading; {flow['included']} included. Relevance: 3 = bears directly on an atlas hypothesis or novelty claim, 2 = important context, 1 = background.
4. Verification and output (`scripts/build_records.py`): every DOI resolved against Crossref (DataCite DOIs via doi.org); Crossref `updated-by`, `update-to` and relation fields checked for retractions and corrections; BibTeX generated from Crossref metadata.
5. Rendering (`scripts/evidence_numbers.py`, `scripts/render_docs.py`, `scripts/write_readme.py`): `src/*.src.md` templates are rendered with author-year links taken from records.tsv and numbers taken from evidence_numbers.json, which reads the final atlas track outputs (`analyses/telomere_cis`, `analyses/junctions`, `analyses/allele_ages`), the catalogue (`data/catalogue`, v0.1 revision 2), the Tan et al. 2026 cross-check (`data/validation/tan2026_summary.json`, `tan2026_cohort_forms.tsv`), the query-tool validation tables (`tests/validation`, quoted if present) and the analysis home. Analysis-track values were computed on the earlier per-window form calls; catalogue values use the revision-2 paralog-family labels.

## Reproduce

From the repository root:

```
source work/config.local.sh          # git-ignored; sets ATLAS_INPUT_DIR (analysis home)
export OPENALEX_API_KEY=...          # OpenAlex API key
literature/reproduce.sh              # uses cached API responses in work/literature if present
```

Without the cache in `work/literature/` (git-ignored), a run queries the live services, and the candidate set can differ from this release as the indexes change. The included set is fixed by `curation.tsv`; a DOI in curation.tsv that is not retrieved stops the build. Abstracts are used only for screening and are not redistributed.

## Limitations

A single reviewer did the screening and annotation, with no duplicate screening. The keyword pre-screen (score 6 or more) excluded {flow['prescreen_excluded']} low-scoring candidates without title review. {npre} included records are preprints (flagged in records.tsv) and are not peer reviewed. One seed (Flint et al. 1994, Am J Hum Genet, telomere healing) has no registered DOI and is not included.
"""
open(os.path.join(LIT, "README.md"), "w").write(txt)
print("README.md written")
