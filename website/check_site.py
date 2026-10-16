#!/usr/bin/env python3
"""QC for the generated site in docs/.

Usage (from the repository root):  python website/check_site.py [--out work/site_check/report.json]
Requires: Python 3.9+, html5lib.

Checks
  1. every HTML page parses with html5lib in strict mode (no parse errors), has <html lang>, <title>, one <h1>,
     unique ids, and alt text on every <img>;
  2. every internal href/src resolves to a file in docs/ (and #fragments to an id in the target page);
  3. links to the GitHub repository (blob/tree/raw main) point at files or folders that exist in this checkout;
  4. no external resources are loaded at runtime (script/link/img src must be relative);
  5. every JSON file and every data/*.js payload (window.X = {...};) parses;
  6. no absolute local paths or user names leak into the output;
  7. size report (total, per page with its data file) and the 15 MB budget.
Exit status 1 if any check fails.
"""
import argparse
import json
import os
import re
import sys
from html.parser import HTMLParser
from urllib.parse import urlparse, unquote

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DOCS = os.path.join(ROOT, "docs")
REPO = "https://github.com/jlanej/human-subtelomere-atlas/"
BUDGET = 15e6
_HOME_DIRS = "|".join(["/" + "Users" + "/[A-Za-z]", "/" + "home" + "/[a-z]+/", "[A-Z]:" + r"\\" + "Users"])
LEAK = re.compile("(" + _HOME_DIRS + "|file:///)")


class Collect(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links, self.ids, self.imgs, self.h1, self.title, self.lang = [], [], [], 0, False, None
        self.res = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.append(a["id"])
        if tag == "html":
            self.lang = a.get("lang")
        if tag == "h1":
            self.h1 += 1
        if tag == "title":
            self.title = True
        if tag == "a" and a.get("href"):
            self.links.append(("a", a["href"]))
        if tag in ("script", "img") and a.get("src"):
            self.links.append((tag, a["src"]))
            self.res.append((tag, a["src"]))
        if tag == "link" and a.get("href"):
            self.links.append(("link", a["href"]))
            self.res.append(("link", a["href"]))
        if tag == "img":
            self.imgs.append(a.get("alt"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "work", "site_check", "report.json"))
    args = ap.parse_args()
    import html5lib

    errors, warnings = [], []
    pages, ids_by_page = {}, {}
    html_files = []
    for dp, _, fs in os.walk(DOCS):
        for f in fs:
            if f.endswith(".html"):
                html_files.append(os.path.relpath(os.path.join(dp, f), DOCS).replace(os.sep, "/"))
    html_files.sort()
    parsed = {}
    for rel in html_files:
        text = open(os.path.join(DOCS, rel), encoding="utf-8").read()
        p = html5lib.HTMLParser(strict=True, namespaceHTMLElements=False)
        try:
            p.parse(text)
        except Exception as ex:  # strict mode raises on the first parse error
            errors.append(f"{rel}: HTML parse error: {ex}")
        c = Collect()
        c.feed(text)
        parsed[rel] = c
        ids_by_page[rel] = set(c.ids)
        dup = [i for i in set(c.ids) if c.ids.count(i) > 1]
        if dup:
            errors.append(f"{rel}: duplicate ids {dup[:5]}")
        if not c.lang:
            errors.append(f"{rel}: <html> without lang")
        if not c.title:
            errors.append(f"{rel}: no <title>")
        if c.h1 != 1:
            errors.append(f"{rel}: {c.h1} <h1> elements")
        miss_alt = sum(1 for a in c.imgs if not a)
        if miss_alt:
            errors.append(f"{rel}: {miss_alt} <img> without alt text")
        if LEAK.search(text):
            errors.append(f"{rel}: local path or user name pattern found: {LEAK.search(text).group(0)}")
        pages[rel] = {"bytes": len(text.encode("utf-8"))}
    n_int = n_repo = n_ext = 0
    ext_domains = {}
    for rel, c in parsed.items():
        base = os.path.dirname(rel)
        for tag, url in c.links:
            u = urlparse(url)
            if u.scheme in ("http", "https"):
                if tag in ("script", "img", "link"):
                    errors.append(f"{rel}: external runtime resource {url}")
                if url.startswith(REPO):
                    n_repo += 1
                    m = re.match(re.escape(REPO) + r"(blob|tree|raw)/main/([^#?]*)", url)
                    if m:
                        path = unquote(m.group(2))
                        if not os.path.exists(os.path.join(ROOT, path)):
                            warnings.append(f"{rel}: repository link to missing path {path}")
                else:
                    n_ext += 1
                    ext_domains[u.netloc] = ext_domains.get(u.netloc, 0) + 1
                continue
            if u.scheme in ("mailto", "data", "javascript"):
                if u.scheme == "javascript":
                    errors.append(f"{rel}: javascript: link")
                continue
            if url.startswith("/"):
                errors.append(f"{rel}: root-relative link {url} (breaks on project pages and file://)")
                continue
            n_int += 1
            target = os.path.normpath(os.path.join(base, unquote(u.path))).replace(os.sep, "/") if u.path else rel
            full = os.path.join(DOCS, target)
            if not os.path.exists(full):
                errors.append(f"{rel}: broken link {url}")
                continue
            if os.path.isdir(full):
                errors.append(f"{rel}: link to a directory {url} (needs index.html for file://)")
                continue
            if u.fragment and target.endswith(".html"):
                frag = unquote(u.fragment)
                js_anchor = target.startswith("arms/") and re.match(r"^(SA|HG)-", frag)
                if frag not in ids_by_page.get(target, set()) and not js_anchor:
                    errors.append(f"{rel}: missing anchor #{frag} in {target}")
    # JSON and JS payloads
    data_ok = 0
    for dp, _, fs in os.walk(DOCS):
        for f in fs:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, DOCS).replace(os.sep, "/")
            try:
                if f.endswith(".json"):
                    json.load(open(p))
                    data_ok += 1
                elif f.endswith(".js") and rel.startswith("data/"):
                    t = open(p, encoding="utf-8").read().strip()
                    m = re.match(r"^window\.(\w+) = (.*);$", t, re.S)
                    if not m:
                        raise ValueError("unexpected payload wrapper")
                    json.loads(m.group(2))
                    data_ok += 1
                if f.endswith((".js", ".css", ".json")) and LEAK.search(open(p, encoding="utf-8").read()):
                    errors.append(f"{rel}: local path or user name pattern found")
            except Exception as ex:
                errors.append(f"{rel}: data does not parse: {ex}")
    total = 0
    for dp, _, fs in os.walk(DOCS):
        for f in fs:
            total += os.path.getsize(os.path.join(dp, f))
    for rel in pages:
        if rel.startswith("arms/") and rel != "arms/index.html":
            arm = os.path.basename(rel)[:-5]
            dj = os.path.join(DOCS, "data", "arms", arm + ".js")
            pages[rel]["data_bytes"] = os.path.getsize(dj) if os.path.exists(dj) else None
    if total > BUDGET:
        errors.append(f"docs/ is {total / 1e6:.2f} MB > {BUDGET / 1e6:.0f} MB budget")
    rep = {"pages": len(html_files), "total_bytes": total, "internal_links": n_int, "repo_links": n_repo,
           "external_links": n_ext, "external_domains": ext_domains, "data_files_parsed": data_ok,
           "errors": errors, "warnings": warnings, "page_sizes": pages}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(rep, open(args.out, "w"), indent=1)
    print(f"pages {len(html_files)}; docs/ {total / 1e6:.2f} MB; internal links {n_int}; repository links {n_repo}; "
          f"other external links {n_ext} ({', '.join(sorted(ext_domains))}); data files parsed {data_ok}")
    print(f"errors {len(errors)}; warnings {len(warnings)}")
    for e in errors[:40]:
        print("ERROR", e)
    for w in warnings[:40]:
        print("WARN ", w)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
