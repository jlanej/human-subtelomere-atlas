#!/usr/bin/env python3
"""Render literature/REVIEW.md and literature/NOVELTY.md from literature/src/*.src.md.

[[doi]] tokens become markdown links '[First-author Year](https://doi.org/doi)' using author and
year from literature/records.tsv (Crossref metadata); a DOI absent from records.tsv is an error.
{{key}} tokens are filled from literature/evidence_numbers.json (or record counts n_records/n_rel3).
Also writes literature/cited_dois.txt (DOIs cited in the rendered documents).
"""
import csv, json, os, re, sys
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIT = os.path.join(REPO, "literature")
recs = {r["doi"]: r for r in csv.DictReader(open(os.path.join(LIT, "records.tsv")), delimiter="\t")}
ev = {k: v["value"] for k, v in json.load(open(os.path.join(LIT, "evidence_numbers.json"))).items() if not k.startswith("_")}
ev["n_records"] = len(recs); ev["n_rel3"] = sum(r["relevance"] == "3" for r in recs.values())
ev["n_candidates"] = f'{json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "screening_flow.json")))["candidates"]:,}'
cited = set(); missing = set()
def cite(m):
    doi = m.group(1).strip().lower()
    r = recs.get(doi)
    if not r:
        missing.add(doi); return m.group(0)
    cited.add(doi)
    href = "https://doi.org/" + doi.replace("(", "%28").replace(")", "%29")
    return f"[{r['first_author']} {r['year']}]({href})"
def num(m):
    k = m.group(1)
    if k not in ev:
        missing.add("number:" + k); return m.group(0)
    return str(ev[k])
for name in ("REVIEW", "NOVELTY"):
    src = open(os.path.join(LIT, "src", name + ".src.md")).read()
    out = re.sub(r"\{\{([^}]+)\}\}", num, re.sub(r"\[\[([^\]]+)\]\]", cite, src))
    open(os.path.join(LIT, name + ".md"), "w").write(out)
if missing:
    sys.exit("unresolved tokens: " + ", ".join(sorted(missing)))
open(os.path.join(LIT, "cited_dois.txt"), "w").write("\n".join(sorted(cited)) + "\n")
print(f"rendered REVIEW.md and NOVELTY.md; distinct DOIs cited: {len(cited)}")
