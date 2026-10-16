#!/usr/bin/env python3
"""Fill abstracts missing from OpenAlex using NCBI PubMed E-utilities (anonymous, <=3 req/s).

Reads literature/screening.tsv (retained works) and work/literature/candidates.jsonl; for retained
works with a DOI and no OpenAlex abstract, looks up the PMID by DOI and fetches the abstract.
Writes work/literature/pubmed_abstracts.json {doi: {pmid, abstract}} (git-ignored; abstracts are
used only for screening/annotation and are not redistributed).
"""
import json, os, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET, csv
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.path.join(REPO, "work", "literature")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(url):
    for a in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                time.sleep(0.4); return r.read()
        except Exception:
            time.sleep(2 * (a + 1))
    return b""
def main():
    out_path = os.path.join(WORK, "pubmed_abstracts.json")
    out = json.load(open(out_path)) if os.path.exists(out_path) else {}
    keep = {r["openalex_id"] for r in csv.DictReader(open(os.path.join(REPO, "literature", "screening.tsv")), delimiter="\t") if r["prescreen"] == "retain"}
    todo = []
    for line in open(os.path.join(WORK, "candidates.jsonl")):
        r = json.loads(line)
        if r["openalex_id"] in keep and r["doi"] and not r["abstract"] and r["doi"] not in out:
            todo.append(r["doi"])
    for doi in todo:
        js = get(EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": f"{doi}[doi]", "retmode": "json"}))
        try:
            ids = json.loads(js)["esearchresult"]["idlist"]
        except Exception:
            ids = []
        if len(ids) != 1:
            out[doi] = {"pmid": "", "abstract": ""}; continue
        xml = get(EU + "efetch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ids[0], "retmode": "xml"}))
        try:
            root = ET.fromstring(xml)
            ab = " ".join("".join(e.itertext()) for e in root.iter("AbstractText"))
        except Exception:
            ab = ""
        out[doi] = {"pmid": ids[0], "abstract": ab}
    json.dump(out, open(out_path, "w"), indent=0)
    print(f"looked_up={len(todo)} with_abstract={sum(1 for v in out.values() if v['abstract'])} total={len(out)}")
if __name__ == "__main__":
    main()
