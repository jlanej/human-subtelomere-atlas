#!/usr/bin/env python3
"""Retrieve candidate literature records for the subtelomere atlas review.

Inputs : literature/scripts/search_config.json (queries, seed references, citation-expansion seeds)
Outputs: work/literature/openalex_raw/*.json   raw API responses (git-ignored snapshot)
         work/literature/candidates.jsonl       one record per unique work, with abstract (git-ignored)
         literature/candidates.tsv              metadata of all candidates (no abstracts)
         literature/seed_resolution.tsv         how each seed reference was resolved to a DOI

Requires the environment variable OPENALEX_API_KEY (an OpenAlex API key, passed as the
`api_key` query parameter). Crossref is queried anonymously. Standard library only.
Results depend on the live state of OpenAlex/Crossref; the retrieval date is recorded in
literature/candidates.tsv's companion file literature/search_log.json. Re-running uses the
cached raw responses in work/literature/openalex_raw unless --refresh is given.
"""
import json, os, sys, time, hashlib, urllib.parse, urllib.request, datetime, argparse

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIT = os.path.join(REPO, "literature")
WORK = os.path.join(REPO, "work", "literature")
RAW = os.path.join(WORK, "openalex_raw")
CFG = os.path.join(LIT, "scripts", "search_config.json")
OA = "https://api.openalex.org"
SELECT = "id,doi,title,publication_year,primary_location,cited_by_count,abstract_inverted_index,authorships,type,referenced_works"
UA = "human-subtelomere-atlas-literature/1.0"


def get_json(url, cache_name, refresh):
    path = os.path.join(RAW, cache_name + ".json")
    if os.path.exists(path) and not refresh:
        with open(path) as fh:
            return json.load(fh)
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
            with open(path, "w") as fh:
                json.dump(data, fh)
            time.sleep(0.12)
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                with open(path, "w") as fh:
                    json.dump({"_error": 404}, fh)
                return {"_error": 404}
            if e.code in (429, 500, 502, 503):
                time.sleep(2 * (attempt + 1)); continue
            raise
        except Exception:
            time.sleep(2 * (attempt + 1))
    return {"_error": "failed"}


def oa_url(path, params):
    p = dict(params); p["api_key"] = os.environ["OPENALEX_API_KEY"]
    return f"{OA}{path}?{urllib.parse.urlencode(p)}"


def key(s):
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def abstract(inv):
    if not inv:
        return ""
    pos = {}
    for w, idxs in inv.items():
        for i in idxs:
            pos[i] = w
    return " ".join(pos[i] for i in sorted(pos))


def compact(w):
    loc = w.get("primary_location") or {}
    src = (loc.get("source") or {})
    auths = [a.get("author", {}).get("display_name", "") for a in (w.get("authorships") or [])]
    doi = (w.get("doi") or "").replace("https://doi.org/", "").lower()
    return {"openalex_id": (w.get("id") or "").replace("https://openalex.org/", ""),
            "doi": doi, "title": w.get("title") or "", "year": w.get("publication_year"),
            "venue": src.get("display_name") or "", "type": w.get("type") or "",
            "cited_by": w.get("cited_by_count") or 0, "authors": auths,
            "abstract": abstract(w.get("abstract_inverted_index"))}


def crossref_resolve(ref, refresh):
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode({"query.bibliographic": ref, "rows": 3})
    d = get_json(url, "crossref_" + key(ref), refresh)
    items = (d.get("message") or {}).get("items") or []
    if not items:
        return None
    it = items[0]
    return {"doi": it.get("DOI", "").lower(), "title": (it.get("title") or [""])[0],
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0], "score": it.get("score")}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--refresh", action="store_true"); a = ap.parse_args()
    if "OPENALEX_API_KEY" not in os.environ:
        sys.exit("OPENALEX_API_KEY not set; OpenAlex requires an API key")
    os.makedirs(RAW, exist_ok=True)
    cfg = json.load(open(CFG))
    recs, hits = {}, {}   # openalex_id -> record ; openalex_id -> set of provenance strings

    def add(w, prov):
        if not w or "_error" in w or not w.get("id"):
            return None
        c = compact(w); oid = c["openalex_id"]
        recs.setdefault(oid, c); hits.setdefault(oid, set()).add(prov)
        return oid

    # 1. keyword queries
    for q in cfg["queries"]:
        url = oa_url("/works", {"search": q["q"], "per-page": cfg["per_query"], "select": SELECT})
        d = get_json(url, "q_" + key(q["q"]), a.refresh)
        for w in d.get("results", []):
            add(w, f"query:{q['tag']}:{q['q']}")
    # 2. seeds
    seedrows = []; seed_oids = {}
    for s in cfg["seeds"]:
        doi, how, cr = s.get("doi"), "pinned", None
        if not doi:
            cr = crossref_resolve(s["ref"], a.refresh)
            doi, how = (cr or {}).get("doi"), "crossref_top_hit"
        w = get_json(oa_url(f"/works/doi:{doi}", {"select": SELECT}), "doi_" + key(doi or "none"), a.refresh) if doi else {}
        oid = add(w, f"seed:{s['tag']}")
        seed_oids[s["ref"]] = oid
        seedrows.append([s["tag"], s["ref"], how, doi or "", (cr or {}).get("title", ""), str((cr or {}).get("score", "")), oid or "NOT_IN_OPENALEX"])
    # 3. one-step citation expansion (forward: most-cited citing works; backward: references)
    for pref in cfg["expand"]:
        oid = next((v for k, v in seed_oids.items() if k.startswith(pref) and v), None)
        if not oid:
            continue
        d = get_json(oa_url("/works", {"filter": f"cites:{oid}", "sort": "cited_by_count:desc", "per-page": 25, "select": SELECT}), "fwd_" + oid, a.refresh)
        for w in d.get("results", []):
            add(w, f"cites:{oid}")
        d = get_json(oa_url("/works", {"filter": f"cites:{oid}", "sort": "publication_date:desc", "per-page": 25, "select": SELECT}), "fwdrecent_" + oid, a.refresh)
        for w in d.get("results", []):
            add(w, f"cites_recent:{oid}")
        d = get_json(oa_url("/works", {"filter": f"cited_by:{oid}", "per-page": 50, "select": SELECT}), "bwd_" + oid, a.refresh)
        for w in d.get("results", []):
            add(w, f"cited_by:{oid}")

    with open(os.path.join(LIT, "seed_resolution.tsv"), "w") as fh:
        fh.write("tag\tseed_reference\tresolution\tdoi\tcrossref_title\tcrossref_score\topenalex_id\n")
        for r in seedrows:
            fh.write("\t".join(str(x).replace("\t", " ") for x in r) + "\n")
    order = sorted(recs, key=lambda o: (-len(hits[o]), -recs[o]["cited_by"], o))
    with open(os.path.join(WORK, "candidates.jsonl"), "w") as fh:
        for o in order:
            r = dict(recs[o]); r["provenance"] = sorted(hits[o]); fh.write(json.dumps(r) + "\n")
    with open(os.path.join(LIT, "candidates.tsv"), "w") as fh:
        fh.write("openalex_id\tdoi\tyear\tcited_by\ttype\tvenue\ttitle\tn_hits\tprovenance\n")
        for o in order:
            r = recs[o]
            fh.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in
                     [o, r["doi"], r["year"], r["cited_by"], r["type"], r["venue"], r["title"], len(hits[o]), ";".join(sorted(hits[o]))]) + "\n")
    mt = [os.path.getmtime(os.path.join(RAW, f)) for f in os.listdir(RAW)]
    fmt = lambda t: datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    log = {"retrieved_utc_first": fmt(min(mt)), "retrieved_utc_last": fmt(max(mt)), "n_raw_responses": len(mt),
           "openalex_endpoint": OA + "/works", "crossref_endpoint": "https://api.crossref.org/works",
           "n_queries": len(cfg["queries"]), "per_query": cfg["per_query"], "n_seeds": len(cfg["seeds"]),
           "n_expanded": len(cfg["expand"]), "n_unique_candidates": len(recs)}
    json.dump(log, open(os.path.join(LIT, "search_log.json"), "w"), indent=1)
    print(json.dumps(log))


if __name__ == "__main__":
    main()
