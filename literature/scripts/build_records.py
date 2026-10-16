#!/usr/bin/env python3
"""Build literature/records.tsv, literature/references.bib and literature/doi_verification.tsv.

Inputs : literature/curation.tsv (manual inclusion, tags, relevance, annotation)
         work/literature/candidates.jsonl (OpenAlex metadata from fetch_literature.py)
Every DOI is resolved against Crossref (https://api.crossref.org/works/<doi>); DOIs not registered
with Crossref (e.g. DataCite/Zenodo) are checked by a doi.org HEAD request. Crossref responses are
cached in work/literature/crossref_cache/ (git-ignored). Crossref 'update-to' / relation fields are
inspected for retractions. Standard library only.
"""
import csv, json, os, re, time, unicodedata, urllib.parse, urllib.request, hashlib
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIT = os.path.join(REPO, "literature"); WORK = os.path.join(REPO, "work", "literature")
CACHE = os.path.join(WORK, "crossref_cache"); os.makedirs(CACHE, exist_ok=True)
UA = "human-subtelomere-atlas-literature/1.0"

def crossref(doi):
    p = os.path.join(CACHE, hashlib.sha1(doi.encode()).hexdigest()[:16] + ".json")
    if os.path.exists(p):
        return json.load(open(p))
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    out = None
    for a in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60) as r:
                out = {"status": 200, "message": json.load(r)["message"]}; break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                out = {"status": 404}; break
            time.sleep(2 * (a + 1))
        except Exception:
            time.sleep(2 * (a + 1))
    if out is None:
        out = {"status": "error"}
    if out["status"] == 404:  # not a Crossref DOI: check the DOI resolver
        try:
            req = urllib.request.Request("https://doi.org/" + doi, method="HEAD", headers={"User-Agent": UA})
            class NoRedirect(urllib.request.HTTPRedirectHandler):
                def redirect_request(self, *a, **k): return None
            op = urllib.request.build_opener(NoRedirect)
            try:
                op.open(req, timeout=30); out["doi_org"] = 200
            except urllib.error.HTTPError as e:
                out["doi_org"] = e.code
        except Exception as e:
            out["doi_org"] = "error"
    json.dump(out, open(p, "w")); time.sleep(0.2)
    return out

def ascii_(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()

def bibesc(s):
    s = re.sub(r"<[^>]+>", "", s or "")
    return s.replace("&amp;", "\\&").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")

def main():
    cand = {}
    for line in open(os.path.join(WORK, "candidates.jsonl")):
        r = json.loads(line)
        if r["doi"]:
            cand.setdefault(r["doi"], r)
    cur = [r for r in csv.DictReader((l for l in open(os.path.join(LIT, "curation.tsv")) if not l.startswith("#")), delimiter="\t")]
    recs, ver, bib, keys = [], [], [], {}
    for c in cur:
        doi = c["doi"].strip().lower(); oa = cand[doi]
        cr = crossref(doi); m = cr.get("message") or {}
        retracted = None
        if cr["status"] == 200:
            # a retracted work carries 'updated-by' entries (type retraction/withdrawal) pointing at the notice;
            # 'update-to' is carried by the notice itself, so a DOI with it is a notice, not the work.
            ups = [u.get("type", "").lower() for u in m.get("updated-by", []) + m.get("update-to", [])]
            rel = m.get("relation", {})
            retracted = any(("retract" in u) or ("withdraw" in u) for u in ups) or ("is-retracted-by" in rel) or ("retraction" in (m.get("type") or ""))
            updates_seen = sorted(set(ups))
            ok, registry = True, "crossref"
        else:
            updates_seen = []
            ok = cr.get("doi_org") in (200, 301, 302, 303, 307, 308)
            registry = "non-crossref(doi.org)" if ok else "unresolved"
        title = ((m.get("title") or [None])[0] or oa["title"]).strip()
        title = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", title))
        year = None
        for k in ("published-print", "published-online", "issued", "posted"):
            dp = (m.get(k) or {}).get("date-parts")
            if dp and dp[0] and dp[0][0]:
                year = dp[0][0]; break
        year = year or oa["year"]
        venue = ((m.get("container-title") or [None])[0]) or oa["venue"] or (m.get("institution", [{}])[0].get("name") if m.get("institution") else "")
        if m.get("type") == "posted-content" and not venue:
            venue = "bioRxiv" if doi.startswith(("10.1101/", "10.64898/")) else "preprint"
        authors = [(a.get("family") or a.get("name") or "", a.get("given") or "") for a in m.get("author", [])] or \
                  [(n.split()[-1] if n else "", " ".join(n.split()[:-1])) for n in oa["authors"]]
        first = authors[0][0] if authors and authors[0][0] else "Anon"
        short = first + (" et al." if len(authors) > 2 else (" & " + authors[1][0] if len(authors) == 2 else ""))
        base = re.sub(r"[^A-Za-z]", "", ascii_(first)) + str(year)
        w = next((x for x in re.findall(r"[A-Za-z]+", ascii_(title)) if x.lower() not in {"a", "an", "the", "of", "on", "in", "and", "for", "to"}), "x")
        key = base + w.lower(); n = keys.get(key, 0); keys[key] = n + 1
        if n: key += "abcdefghij"[n]
        preprint = m.get("type") == "posted-content" or "biorxiv" in (venue or "").lower() or "medrxiv" in (venue or "").lower()
        recs.append({"doi": doi, "title": title, "first_author": short, "year": year, "venue": venue or "",
                     "openalex_id": oa["openalex_id"], "cited_by": oa["cited_by"], "preprint": int(bool(preprint)),
                     "topic_tags": c["topic_tags"], "relevance": c["relevance"], "annotation": c["annotation"], "bibkey": key})
        ver.append([doi, int(ok), registry, "" if retracted is None else int(retracted), ",".join(updates_seen), str(cr["status"]), title[:120]])
        typ = "article" if not preprint and m.get("type") in ("journal-article", None) and venue else "misc"
        if m.get("type") == "proceedings-article": typ = "inproceedings"
        f = [f"  title = {{{bibesc(title)}}}",
             "  author = {" + " and ".join(bibesc(f"{fa}, {gi}" if gi else fa) for fa, gi in authors) + "}",
             f"  year = {{{year}}}", f"  doi = {{{doi}}}"]
        if typ == "article":
            f.append(f"  journal = {{{bibesc(venue)}}}")
            for k, bk in (("volume", "volume"), ("issue", "number"), ("page", "pages")):
                if m.get(k): f.append(f"  {bk} = {{{bibesc(m[k]).replace('-', '--')}}}")
        else:
            f.append(f"  howpublished = {{{bibesc(venue or registry)}}}")
            if preprint: f.append("  note = {Preprint}")
        bib.append(f"@{typ}{{{key},\n" + ",\n".join(f) + "\n}\n")
    recs.sort(key=lambda r: (-int(r["relevance"]), r["topic_tags"].split(";")[0], r["year"] or 0, r["doi"]))
    cols = ["doi", "title", "first_author", "year", "venue", "openalex_id", "cited_by", "preprint", "topic_tags", "relevance", "annotation", "bibkey"]
    with open(os.path.join(LIT, "records.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(LIT, "doi_verification.tsv"), "w") as fh:
        fh.write("doi\tresolves\tregistry\tretracted\tcrossref_updates\tcrossref_status\tregistered_title\n")
        for v in sorted(ver): fh.write("\t".join(map(str, v)) + "\n")
    with open(os.path.join(LIT, "references.bib"), "w") as fh:
        fh.write("% Bibliography for the human subtelomere atlas literature review; generated by literature/scripts/build_records.py from Crossref metadata.\n\n")
        fh.write("\n".join(sorted(bib)))
    scr = list(csv.DictReader(open(os.path.join(LIT, "screening.tsv")), delimiter="\t"))
    inc = {r["doi"] for r in recs}
    exc = {r["doi"] for r in csv.DictReader(open(os.path.join(LIT, "curation_excluded.tsv")), delimiter="\t")}
    ret = [r for r in scr if r["prescreen"] == "retain"]
    flow = {"candidates": len(scr), "prescreen_excluded": len(scr) - len(ret), "prescreen_retained": len(ret),
            "retained_seeds": sum(r["seed"] == "1" for r in ret),
            "retained_without_doi": sum(1 for r in ret if not r["doi"]),
            "excluded_after_abstract": sum(1 for r in ret if r["doi"] in exc), "excluded_records_no_longer_retrieved_as_seed": len(exc - {r["doi"] for r in ret}),
            "excluded_at_title_screen": sum(1 for r in ret if r["doi"] and r["doi"] not in inc and r["doi"] not in exc),
            "included": len(inc), "included_not_prescreen_retained": len(inc - {r["doi"] for r in ret})}
    json.dump(flow, open(os.path.join(LIT, "screening_flow.json"), "w"), indent=1)
    print(json.dumps(flow))
    print(json.dumps({"records": len(recs), "resolved": sum(v[1] for v in ver), "crossref": sum(v[2] == "crossref" for v in ver),
                      "retracted": sum(1 for v in ver if v[3] == 1), "unresolved": [v[0] for v in ver if not v[1]],
                      "relevance": {k: sum(r["relevance"] == k for r in recs) for k in "123"}, "preprints": sum(r["preprint"] for r in recs)}))
if __name__ == "__main__":
    main()
