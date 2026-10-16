#!/usr/bin/env python3
"""Deterministic title/abstract pre-screen of literature candidates.

Reads work/literature/candidates.jsonl (from fetch_literature.py) and assigns each work a
keyword score. Works retrieved as seeds are always retained. Non-seed works are retained for
manual screening when score >= 6 and they are not flagged as off-topic organisms/fields.
Writes literature/screening.tsv (all candidates with score, flags and pre-screen decision).
Final inclusion and relevance grades are manual and recorded in literature/curation.tsv.
"""
import json, os, re
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "work", "literature", "candidates.jsonl")
OUT = os.path.join(REPO, "literature", "screening.tsv")
POS = {  # pattern: weight
    r"subtelomer": 4, r"telomere[- ]adjacent": 3, r"chromosome ends?\b": 2, r"telomere length": 2,
    r"segmental duplication": 2, r"interchromosomal|inter-chromosomal|non-homologous chromosomes": 2,
    r"telomere capture|telomere healing|chromosome healing": 4, r"break-induced replication": 3,
    r"interstitial telomer|intrachromosomal telomer|telomeric-like": 4, r"\bD4Z4\b|4q35|10q26|FSHD|facioscapulohumeral": 3,
    r"pseudoautosomal": 3, r"acrocentric": 3, r"telomere-to-telomere|\bT2T\b|complete (human )?genome|complete sequence": 2,
    r"pangenome": 2, r"\bHG002\b|CHM13": 2, r"long[- ]read|nanopore|HiFi|PacBio": 1, r"variant repeat|TVR|TTAGGG": 2,
    r"concerted evolution|molecular drive|gene conversion|unequal crossing|homogeni": 3,
    r"olfactory receptor": 1, r"\bhuman\b|homo sapiens": 1, r"pedigree|inherit": 1, r"haplotype|allele": 1,
    r"great ape|primate|chimpanzee": 1, r"STELA|single telomere length": 4, r"translocation": 1,
}
NEG = r"\b(plant|arabidopsis|rice|wheat|maize|soybean|cotton|tomato|grape|barley|watermelon|citrullus|fungal|fungus|candida|plasmodium|trypanosom|drosophila|nematode|caenorhabditis|bird|chicken|fish|insect|bacteri|viral|virus|rhinovirus|covid|sars-cov|cattle|horse|pig|dog|mouse genome|crop)\w*"
def main():
    rows = []
    for line in open(SRC):
        r = json.loads(line)
        text = (r["title"] + " " + r["abstract"]).lower()
        score = sum(w for p, w in POS.items() if re.search(p, text, re.I))
        neg = sorted(set(m.group(0) for m in re.finditer(NEG, (r["title"]).lower())))
        seed = any(p.startswith("seed:") for p in r["provenance"])
        keep = seed or (score >= 6 and not neg and r["type"] not in ("erratum", "paratext", "peer-review"))
        rows.append((r, score, neg, seed, keep))
    rows.sort(key=lambda x: (not x[4], -x[1], -x[0]["cited_by"], x[0]["openalex_id"]))
    with open(OUT, "w") as fh:
        fh.write("openalex_id\tdoi\tyear\tcited_by\tscore\tseed\toff_topic_flags\tprescreen\ttitle\n")
        for r, s, neg, seed, keep in rows:
            fh.write("\t".join(map(str, [r["openalex_id"], r["doi"], r["year"], r["cited_by"], s, int(seed), ",".join(neg),
                                         "retain" if keep else "exclude", r["title"].replace("\t", " ").replace("\n", " ")])) + "\n")
    print(f"candidates={len(rows)} retained={sum(x[4] for x in rows)} seeds={sum(x[3] for x in rows)}")
if __name__ == "__main__":
    main()
