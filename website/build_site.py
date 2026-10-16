#!/usr/bin/env python3
"""Build the static GitHub Pages site in docs/ from data/, analyses/, figures/ and literature/.

Usage (from the repository root):
    python pipeline/figures_atlas.py      # atlas figures + figures/atlas/PALETTE.tsv (needed first)
    python website/build_site.py          # writes docs/

The output is plain HTML/CSS/JS with relative links only: it can be opened from file:// or served from
https://jlanej.github.io/human-subtelomere-atlas/ without a build step. No external resources are loaded at runtime.
Requirements: Python 3.9+, pandas; Pillow (optional, used to make web-sized PNG copies of figures).
Files in docs/ are overwritten in place; nothing is deleted.
"""
import datetime
import hashlib
import html
import json
import os
import re
import shutil
import sys
from collections import Counter, OrderedDict

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DOCS = os.path.join(ROOT, "docs")
STATIC = os.path.join(ROOT, "website", "static")
CAT = os.path.join(ROOT, "data", "catalogue")
REPO_URL = "https://github.com/jlanej/human-subtelomere-atlas"
BLOB = REPO_URL + "/blob/main/"
TREE = REPO_URL + "/tree/main/"
RAW = REPO_URL + "/raw/main/"
SITE_URL = "https://jlanej.github.io/human-subtelomere-atlas/"
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POOLED_HOME = {"acroP": "15p", "PAR1": "Xp", "PAR2": "Xq"}
TAN = "10.64898/2026.09.03.748401"
GUARRACINO = "10.64898/2026.07.10.737660"
TAN_ZENODO = "10.5281/zenodo.21807156"

esc = html.escape


def rp(*p):
    return os.path.join(ROOT, *p)


def write(path, text):
    full = os.path.join(DOCS, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as fh:
        fh.write(text)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def human_bytes(n):
    n = float(n)
    for unit in ["B", "kB", "MB", "GB"]:
        if n < 1000 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1000


def load_json(path, default=None):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def fnum(x, d=0):
    """Comma-grouped integers for prose."""
    if d == 0:
        return f"{int(round(x)):,}"
    return f"{x:,.{d}f}"


# ------------------------------------------------------------------ images
def copy_image(src_rel, max_w=1600):
    """Copy a repository image into docs/assets/img/ (web-sized when Pillow is available). Returns docs-relative path."""
    src = rp(src_rel)
    name = src_rel.replace("/", "__")
    dst_rel = "assets/img/" + name
    dst = os.path.join(DOCS, dst_rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        from PIL import Image
        im = Image.open(src)
        if im.width > max_w:
            im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
        if im.mode in ("RGBA", "LA", "P"):
            bg = Image.new("RGB", im.size, (255, 255, 255))
            im2 = im.convert("RGBA")
            bg.paste(im2, mask=im2.split()[-1])
            im = bg
        im.save(dst, optimize=True)
        if os.path.getsize(dst) >= os.path.getsize(src):
            shutil.copyfile(src, dst)
    except Exception:
        shutil.copyfile(src, dst)
    return dst_rel


# ------------------------------------------------------------------ minimal markdown -> HTML
class MD:
    """Small Markdown renderer for the repository's own documents: headings, paragraphs, nested lists, pipe tables,
    fenced code, emphasis, inline code, links and images. Relative links are rewritten to GitHub URLs; relative images
    are copied into docs/assets/img/."""

    def __init__(self, src_rel, prefix, hshift=1, id_prefix=""):
        self.src_dir = os.path.dirname(src_rel)
        self.prefix = prefix
        self.hshift = hshift
        self.id_prefix = id_prefix

    def resolve(self, url):
        if re.match(r"^[a-z]+:", url) or url.startswith("#"):
            return url, None
        path, _, frag = url.partition("#")
        full = os.path.normpath(os.path.join(self.src_dir, path)).replace(os.sep, "/")
        return full, frag

    def link(self, url):
        full, frag = self.resolve(url)
        if frag is None:
            return full
        kind = "tree/main/" if os.path.isdir(rp(full)) else "blob/main/"
        return REPO_URL + "/" + kind + full + (("#" + frag) if frag else "")

    def image(self, alt, url):
        full, frag = self.resolve(url)
        if frag is None:
            return f'<img src="{esc(full)}" alt="{esc(alt)}">'
        if os.path.exists(rp(full)):
            dst = copy_image(full)
            return (f'<a href="{esc(self.prefix + dst)}"><img src="{esc(self.prefix + dst)}" alt="{esc(alt)}" '
                    f'loading="lazy"></a>')
        return f'<span class="muted">[image not found: {esc(full)}]</span>'

    def inline(self, s):
        codes = []

        def keep(m):
            codes.append("<code>" + esc(m.group(1)) + "</code>")
            return f"\x00{len(codes) - 1}\x00"

        s = re.sub(r"`([^`]+)`", keep, s)
        imgs = []

        def keep_img(m):
            imgs.append(self.image(m.group(1), m.group(2)))
            return f"\x01{len(imgs) - 1}\x01"

        s = re.sub(r"!\[([^\]]*)\]\(([^)\s]+)\)", keep_img, s)
        links = []

        def keep_link(m):
            links.append((m.group(1), self.link(m.group(2))))
            return f"\x02{len(links) - 1}\x02"

        s = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", keep_link, s)
        s = esc(s, quote=False)
        for t in ["sub", "sup", "br", "em", "strong"]:
            s = s.replace(f"&lt;{t}&gt;", f"<{t}>").replace(f"&lt;/{t}&gt;", f"</{t}>")
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"(?<![\w*])\*(?!\s)([^*]+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", s)
        s = re.sub(r"(?<![\w])_(?!\s)([^_]+?)(?<!\s)_(?![\w])", r"<em>\1</em>", s)
        s = re.sub(r"(https?://[^\s<)]+[^\s<).,;])", lambda m: f'<a href="{m.group(1)}">{m.group(1)}</a>'
                   if "\x02" not in m.group(1) else m.group(1), s)
        s = re.sub(r"\x02(\d+)\x02", lambda m: f'<a href="{esc(links[int(m.group(1))][1])}">'
                   f'{self.inline_text(links[int(m.group(1))][0])}</a>', s)
        s = re.sub(r"\x01(\d+)\x01", lambda m: imgs[int(m.group(1))], s)
        s = re.sub(r"\x00(\d+)\x00", lambda m: codes[int(m.group(1))], s)
        return s

    def inline_text(self, t):
        t = esc(t, quote=False)
        t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
        return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)

    def slug(self, text):
        return self.id_prefix + re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:60]

    def render(self, text):
        lines = text.split("\n")
        out, i, n = [], 0, len(lines)
        while i < n:
            ln = lines[i]
            if not ln.strip():
                i += 1
                continue
            if ln.strip().startswith("```"):
                j = i + 1
                buf = []
                while j < n and not lines[j].strip().startswith("```"):
                    buf.append(lines[j])
                    j += 1
                out.append("<pre><code>" + esc("\n".join(buf)) + "</code></pre>")
                i = j + 1
                continue
            m = re.match(r"^(#{1,6})\s+(.*)$", ln)
            if m:
                lev = min(6, len(m.group(1)) + self.hshift)
                txt = m.group(2).strip()
                out.append(f'<h{lev} id="{self.slug(re.sub(r"[`*]", "", txt))}">{self.inline(txt)}</h{lev}>')
                i += 1
                continue
            if ln.lstrip().startswith("|") and i + 1 < n and re.match(r"^\s*\|?\s*:?-{2,}", lines[i + 1]):
                head = self.cells(ln)
                j = i + 2
                body = []
                while j < n and lines[j].lstrip().startswith("|"):
                    body.append(self.cells(lines[j]))
                    j += 1
                num = [all(re.match(r"^[-+]?[\d.,]+(e-?\d+)?%?$", r[k].strip() or "0") for r in body if k < len(r))
                       for k in range(len(head))]
                t = ['<div class="table-wrap free"><table><thead><tr>']
                t += [f'<th scope="col"{" class=num" if num[k] else ""}>{self.inline(h)}</th>'
                      for k, h in enumerate(head)]
                t.append("</tr></thead><tbody>")
                for r in body:
                    t.append("<tr>" + "".join(f'<td{" class=num" if k < len(num) and num[k] else ""}>'
                                              f'{self.inline(c)}</td>' for k, c in enumerate(r)) + "</tr>")
                t.append("</tbody></table></div>")
                out.append("".join(t))
                i = j
                continue
            if re.match(r"^\s*([-*]|\d+\.)\s+", ln):
                j, html_list = self.list_block(lines, i)
                out.append(html_list)
                i = j
                continue
            if ln.startswith(">"):
                buf = []
                while i < n and lines[i].startswith(">"):
                    buf.append(lines[i].lstrip("> "))
                    i += 1
                out.append("<blockquote>" + self.inline(" ".join(buf)) + "</blockquote>")
                continue
            buf = []
            while i < n and lines[i].strip() and not re.match(r"^(#{1,6})\s", lines[i]) \
                    and not lines[i].lstrip().startswith("|") and not lines[i].strip().startswith("```") \
                    and not re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                buf.append(lines[i].strip())
                i += 1
            if buf:
                out.append("<p>" + self.inline(" ".join(buf)) + "</p>")
            else:
                i += 1
        return "\n".join(out)

    @staticmethod
    def cells(line):
        s = line.strip()
        if s.startswith("|"):
            s = s[1:]
        if s.endswith("|"):
            s = s[:-1]
        parts, cur, inc = [], "", False
        for ch in s:
            if ch == "`":
                inc = not inc
            if ch == "|" and not inc:
                parts.append(cur.strip())
                cur = ""
            else:
                cur += ch
        parts.append(cur.strip())
        return parts

    def list_block(self, lines, i):
        items = []  # (indent, ordered, text)
        n = len(lines)
        while i < n:
            ln = lines[i]
            m = re.match(r"^(\s*)([-*]|\d+\.)\s+(.*)$", ln)
            if m:
                items.append([len(m.group(1).replace("\t", "    ")), m.group(2)[-1] == ".", m.group(3)])
                i += 1
                continue
            if ln.strip() and items and (ln.startswith("  ") or ln.startswith("\t")) and \
                    not ln.lstrip().startswith("|"):
                items[-1][2] += " " + ln.strip()
                i += 1
                continue
            break
        html_out, stack = [], []
        for ind, ordered, txt in items:
            while stack and ind < stack[-1][0]:
                html_out.append("</li></" + stack.pop()[1] + ">")
            if not stack or ind > stack[-1][0]:
                tag = "ol" if ordered else "ul"
                html_out.append(f"<{tag}>")
                stack.append((ind, tag))
            else:
                html_out.append("</li>")
            html_out.append("<li>" + self.inline(txt))
        while stack:
            html_out.append("</li></" + stack.pop()[1] + ">")
        return i, "".join(html_out)


def md_file(src_rel, prefix, hshift=1, id_prefix=""):
    with open(rp(src_rel), encoding="utf-8") as fh:
        return MD(src_rel, prefix, hshift, id_prefix).render(fh.read())


def md_section(src_rel, heading_regex, prefix, hshift=1, stop_level=None, id_prefix=""):
    """Render the section of a Markdown file starting at the first heading matching heading_regex, up to the next
    heading of the same or higher level."""
    with open(rp(src_rel), encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    start, lev = None, None
    for k, ln in enumerate(lines):
        m = re.match(r"^(#{1,6})\s+(.*)$", ln)
        if m and re.search(heading_regex, m.group(2)):
            start, lev = k, len(m.group(1))
            break
    if start is None:
        return ""
    end = len(lines)
    for k in range(start + 1, len(lines)):
        m = re.match(r"^(#{1,6})\s+", lines[k])
        if m and len(m.group(1)) <= (stop_level or lev):
            end = k
            break
    return MD(src_rel, prefix, hshift, id_prefix).render("\n".join(lines[start:end]))


def md_paragraphs_starting(src_rel, start_regex, prefix):
    """Render every paragraph of a Markdown file whose first line matches start_regex."""
    with open(rp(src_rel), encoding="utf-8") as fh:
        text = fh.read()
    paras = re.split(r"\n\s*\n", text)
    keep = [p for p in paras if re.match(start_regex, p.strip())]
    return MD(src_rel, prefix).render("\n\n".join(keep))


# ------------------------------------------------------------------ page template
NAV = [("index.html", "Home"), ("arms/index.html", "Arms"), ("alleles.html", "Alleles"),
       ("findings/index.html", "Findings"), ("methods.html", "Methods"), ("validation.html", "Validation"),
       ("query.html", "Query tool"), ("downloads.html", "Downloads"), ("literature.html", "Literature")]


def page(path, title, body, active=None, scripts=(), description=None):
    depth = path.count("/")
    pre = "../" * depth
    nav = []
    for href, lab in NAV:
        cur = ' aria-current="page"' if active == href else ""
        nav.append(f'<li><a href="{pre}{href}"{cur}>{lab}</a></li>')
    sc = "\n".join(f'<script src="{pre}{s}"></script>' for s in scripts)
    desc = esc(description or "Human subtelomere atlas: structural alleles of human chromosome ends from HPRC "
                              "Release 2, HG002 Q100 and CHM13.")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)} | Human subtelomere atlas</title>
<meta name="description" content="{desc}">
<link rel="stylesheet" href="{pre}assets/css/site.css">
</head>
<body>
<a class="skip" href="#main">Skip to content</a>
<header class="site"><div class="wrap">
<a class="brand" href="{pre}index.html">Human subtelomere atlas <small>catalogue {esc(VERSION)}</small></a>
<nav class="main" aria-label="Main"><ul>{''.join(nav)}</ul></nav>
</div></header>
<main id="main" class="wrap">
{body}
</main>
<footer class="site"><div class="wrap">
<p>Human subtelomere atlas, catalogue {esc(VERSION)}. Static pages generated by <code>website/build_site.py</code> from the
release tables in <code>data/</code> and the analyses in <code>analyses/</code>.
Source: <a href="{REPO_URL}">{REPO_URL.replace('https://', '')}</a>. No cookies, trackers or external resources.</p>
</div></footer>
<script src="{pre}assets/js/common.js"></script>
{sc}
</body>
</html>
"""


def figure(src_rel, alt, caption_html, prefix, full_link=True):
    dst = copy_image(src_rel)
    links = []
    if full_link:
        links.append(f'<a href="{BLOB}{src_rel}">full resolution</a>')
        pdf = src_rel[:-4] + ".pdf"
        if os.path.exists(rp(pdf)):
            links.append(f'<a href="{BLOB}{pdf}">PDF</a>')
        svg = src_rel[:-4] + ".svg"
        if os.path.exists(rp(svg)):
            links.append(f'<a href="{BLOB}{svg}">SVG</a>')
    extra = f' <span class="muted">({", ".join(links)})</span>' if links else ""
    return (f'<figure><a href="{prefix}{dst}"><img src="{prefix}{dst}" alt="{esc(alt)}" loading="lazy"></a>'
            f'<figcaption>{caption_html}{extra}</figcaption></figure>')


def md_inline(text, src_rel="figures/atlas/CAPTIONS.md", prefix=""):
    return MD(src_rel, prefix).inline(text)


def tier_badge(t):
    return f'<span class="tier tier-{esc(str(t))}">{esc(str(t))}</span>'


def table_html(df, cols, num=(), fmt=None, id_col=None, caption=None, wrap_class="table-wrap"):
    fmt = fmt or {}
    h = [f'<div class="{wrap_class}"><table>']
    if caption:
        h.append(f'<caption class="muted">{caption}</caption>')
    h.append("<thead><tr>" + "".join(f'<th scope="col"{" class=num" if c in num else ""}>{esc(lab)}</th>'
                                      for c, lab in cols) + "</tr></thead><tbody>")
    for _, r in df.iterrows():
        rid = f' id="{esc(str(r[id_col]))}"' if id_col else ""
        cells = []
        for c, _ in cols:
            v = r[c]
            if c in fmt:
                s = fmt[c](v, r)
            elif pd.isna(v) if not isinstance(v, (list, dict)) else False:
                s = ""
            else:
                s = esc(str(v))
            cells.append(f'<td{" class=num" if c in num else ""}>{s}</td>')
        h.append(f"<tr{rid}>" + "".join(cells) + "</tr>")
    h.append("</tbody></table></div>")
    return "".join(h)


# ------------------------------------------------------------------ data
def parse_captions(path):
    out = {}
    if not os.path.exists(path):
        return out
    cur = None
    for ln in open(path, encoding="utf-8"):
        m = re.match(r"^## (\S+)", ln)
        if m:
            cur = m.group(1)
            continue
        if cur and ln.startswith("**"):
            out[cur] = ln.strip()
    return out


def caption_title(cap):
    m = re.match(r"^\*\*(.+?)\*\*", cap or "")
    return re.sub(r"[*_]", "", m.group(1)) if m else ""


VERSION = "v0.1"


# ------------------------------------------------------------------ derived numbers
def ages(D):
    p = rp("analyses", "allele_ages", "results", "origin_verdicts.tsv")
    if not os.path.exists(p):
        return None
    v = pd.read_csv(p, sep="\t")
    col = "age_kyr_median_mu" if "age_kyr_median_mu" in v.columns else None
    vcol = "single_origin" if "single_origin" in v.columns else None
    if not col or not vcol:
        return None
    cl = v[v[vcol].astype(str).str.startswith(("consistent", "clustered"))]
    a = cl[col].dropna()
    return {"n_clustered": int(len(cl)), "n_in_range": int(((a >= 200) & (a <= 1500)).sum()),
            "median": float(a.median()) if len(a) else None}


# ------------------------------------------------------------------ data files for the site
def js_payload(var, obj):
    return f"window.{var} = " + json.dumps(obj, separators=(",", ":"), allow_nan=False) + ";\n"


# ------------------------------------------------------------------ pages
def arm_order(arms):
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    key = {}
    for i, c in enumerate(chroms):
        key[c + "p"], key[c + "q"] = 2 * i, 2 * i + 1
    key.update({"acroP": key["13p"], "PAR1": key["Xp"], "PAR2": key["Xq"]})
    return sorted(arms, key=lambda a: key.get(a, 999))


def arm_stats(D):
    sa, hg, e = D["sa"], D["hg"], D["ends"]
    r2 = e[e.set == "HPRC_R2"]
    rows = []
    for arm in arm_order(sa.arm.unique()):
        g = sa[sa.arm == arm]
        nat = g[g.is_native]
        non = g[~g.is_native].sort_values("r2_freq", ascending=False)
        top = non.iloc[0] if len(non) else None
        rr = r2[r2.catalogue_arm == arm]
        rows.append({"arm": arm, "r2_ends": len(rr), "r2_samples": rr["sample"].nunique(), "n_sa": len(g),
                     "A": int((g.tier == "A").sum()), "B": int((g.tier == "B").sum()), "C": int((g.tier == "C").sum()),
                     "native_freq": float(nat.r2_freq.sum()) if len(nat) else 0.0,
                     "n_hg": int((hg.arm == arm).sum()),
                     "top": f"{top.allele_id} ({top.structure}, {top.r2_freq:.3f})" if top is not None else ""})
    return pd.DataFrame(rows)


def lit_title(doi):
    p = rp("literature", "records.tsv")
    try:
        r = pd.read_csv(p, sep="\t")
        m = r[r.doi.astype(str).str.lower() == doi.lower()]
        if len(m):
            return str(m.iloc[0]["title"])
    except Exception:
        pass
    return ""


def page_arms_index(D, AS):
    pre = "../"
    caps = D["captions"]

    def armlink(v, r):
        return f'<a href="{esc(v)}.html">{esc(v)}</a>'

    tbl = table_html(AS, [("arm", "Arm"), ("r2_ends", "Release 2 ends"), ("r2_samples", "Samples"),
                          ("n_sa", "Structural alleles"), ("A", "Tier A"), ("B", "Tier B"), ("C", "Tier C"),
                          ("native_freq", "Native frequency"), ("n_hg", "Haplogroups"),
                          ("top", "Most frequent non-native allele")],
                     num=("r2_ends", "r2_samples", "n_sa", "A", "B", "C", "native_freq", "n_hg"),
                     fmt={"arm": armlink, "native_freq": lambda v, r: f"{v:.3f}"}, wrap_class="table-wrap free")
    body = f"""
<h1>Arms</h1>
<p>One page per catalogue arm. Ends whose distal sequence is shared between chromosomes are pooled: <code>acroP</code> =
13p, 14p, 15p, 21p, 22p; <code>PAR1</code> = Xp, Yp; <code>PAR2</code> = Xq, Yq. Each page has an interactive mosaic of
the arm's structural alleles, its haplogroups and a filterable table of member ends.</p>
<nav class="arm-nav" aria-label="Arms">{''.join(f'<a href="{a}.html">{a}</a>' for a in AS.arm)}</nav>
{tbl}
<p class="muted">Native frequency = summed Release 2 frequency of the arm's native (all-home) structural alleles.</p>
{figure('figures/atlas/F2_allele_mosaics.png', caption_title(caps.get('F2_allele_mosaics')) +
        ': stacked mosaics of real ends for 20q, 3q, 6p, 8p, 16p and 11p.', md_inline(caps.get('F2_allele_mosaics', '')), pre)}
"""
    write("arms/index.html", page("arms/index.html", "Arms", body, active="arms/index.html"))


ANALYSES = OrderedDict([
    ("telomere-cis", {"dir": "analyses/telomere_cis", "title": "Telomere length is not set in cis by the terminal allele",
                      "short": "Telomere length and the terminal allele (H1)"}),
    ("junctions", {"dir": "analyses/junctions", "title": "Foreign terminal blocks arise by one-way, homology-mediated transfer",
                   "short": "Mechanism at junctions (H2)"}),
    ("allele-ages", {"dir": "analyses/allele_ages", "title": "Common forms are old and shared across populations",
                     "short": "Population distribution, origin and age"}),
])


def analysis_figs(d, pre):
    figdir = rp(d, "figures")
    out = []
    capmd = {}
    for f in sorted(os.listdir(figdir)) if os.path.isdir(figdir) else []:
        if f.endswith("_caption.md"):
            capmd[f[:-11]] = open(os.path.join(figdir, f), encoding="utf-8").read()
    for f in sorted(os.listdir(figdir)) if os.path.isdir(figdir) else []:
        if f.endswith(".png"):
            stem = f[:-4]
            cap = capmd.get(stem)
            cap_html = MD(f"{d}/figures/{f}", pre).render(cap) if cap else \
                f'<code>{esc(d)}/figures/{esc(f)}</code>; caption in the analysis README.'
            alt = f"Figure {stem.replace('_', ' ')} of the {d.split('/')[-1].replace('_', ' ')} analysis"
            out.append(figure(f"{d}/figures/{f}", alt, cap_html, pre))
    return "\n".join(out)


def render_json_generic(obj, depth=0):
    if isinstance(obj, dict):
        rows = "".join(f"<tr><th scope=row>{esc(str(k))}</th><td>{render_json_generic(v, depth + 1)}</td></tr>"
                       for k, v in obj.items())
        return f'<table class="kv"><tbody>{rows}</tbody></table>'
    if isinstance(obj, list):
        if obj and all(isinstance(x, dict) for x in obj):
            keys = []
            for x in obj:
                for k in x:
                    if k not in keys:
                        keys.append(k)
            head = "".join(f"<th scope=col>{esc(str(k))}</th>" for k in keys)
            body = "".join("<tr>" + "".join(f"<td>{render_json_generic(x.get(k, ''), depth + 1)}</td>" for k in keys)
                           + "</tr>" for x in obj)
            return f'<div class="table-wrap free"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
        return esc(", ".join(str(x) for x in obj))
    if isinstance(obj, float):
        return esc(f"{obj:.4g}")
    return esc(str(obj))


def page_query():
    p = rp("atlas", "README.md")
    if os.path.exists(p):
        content = md_file("atlas/README.md", "", hshift=1)
        status = ""
    else:
        content = ""
        status = ('<p class="warn">The query tool documentation (<code>atlas/README.md</code>) is not yet in the '
                  'repository at the time this site was built; rebuild the site after it is added.</p>')
    body = f"""
<h1>Query tool</h1>
<p>The command-line query tool in <a href="{TREE}atlas">atlas/</a> looks up ends, alleles and haplogroups in the
release tables and assigns new sequences to catalogue alleles. Documentation: <a href="{BLOB}atlas/README.md">atlas/README.md</a>.
For browsing, use the <a href="arms/index.html">arm pages</a> and the <a href="alleles.html">allele search</a>.</p>
{status}
{content}
"""
    write("query.html", page("query.html", "Query tool", body, active="query.html"))


def find_licence():
    for f in ["LICENSE", "LICENSE.md", "LICENSE.txt", "LICENCE", "COPYING"]:
        if os.path.exists(rp(f)):
            return f
    return None


def page_downloads(D):
    man = D["manifest"].copy()

    def link(v, r):
        return f'<a href="{RAW}{esc(v)}">{esc(v)}</a>'

    man["size"] = man["bytes"].map(human_bytes)
    man["rows"] = man["rows"].map(lambda v: "" if pd.isna(v) else fnum(v))
    cat_tbl = table_html(man, [("path", "File"), ("size", "Size"), ("rows", "Rows"), ("sha256", "SHA-256"),
                               ("description", "Description")], num=("size", "rows"),
                         fmt={"path": link, "sha256": lambda v, r: f"<code>{esc(str(v))}</code>"},
                         wrap_class="table-wrap free")
    figrows = []
    fd = rp("figures", "atlas")
    for f in sorted(os.listdir(fd)):
        p = os.path.join(fd, f)
        figrows.append({"path": f"figures/atlas/{f}", "size": human_bytes(os.path.getsize(p)), "sha256": sha256(p)})
    fig_tbl = table_html(pd.DataFrame(figrows), [("path", "File"), ("size", "Size"), ("sha256", "SHA-256")],
                         num=("size",), fmt={"path": link, "sha256": lambda v, r: f"<code>{esc(str(v))}</code>"},
                         wrap_class="table-wrap free")
    valrows = []
    for vd in ["data/validation", "tests/validation"]:
        if os.path.isdir(rp(vd)):
            for f in sorted(os.listdir(rp(vd))):
                pth = rp(vd, f)
                if os.path.isfile(pth) and not f.endswith((".pyc", ".DS_Store")):
                    valrows.append({"path": f"{vd}/{f}", "size": human_bytes(os.path.getsize(pth)), "sha256": sha256(pth)})
    val_tbl = table_html(pd.DataFrame(valrows), [("path", "File"), ("size", "Size"), ("sha256", "SHA-256")],
                         num=("size",), fmt={"path": link, "sha256": lambda v, r: f"<code>{esc(str(v))}</code>"},
                         wrap_class="table-wrap free") if valrows else \
        '<p class="warn">Validation outputs were not present when this site was built.</p>'
    lic = find_licence()
    if lic:
        lic_html = f'<p>Licence: see <a href="{BLOB}{lic}">{lic}</a>.</p>'
    else:
        lic_html = ('<p class="warn">No licence file is present in this release of the repository. Until one is added, '
                    'reuse terms for the catalogue files are not declared here. Upstream data keep their own terms: HPRC '
                    'Release 2 and T2T assemblies (public, see the providers), Tan et al. 2026 sequences (CC BY 4.0).</p>')
    prov = D["prov"]
    tools = ", ".join(f"{k} {v}" for k, v in prov.get("tools", {}).items())
    analyses = "".join(f'<li><a href="{TREE}{v["dir"]}/results">{v["dir"]}/results/</a> '
                       f'({esc(v["short"])})</li>' for v in ANALYSES.values())
    body = f"""
<h1>Downloads</h1>
<p>Catalogue {esc(VERSION)}. Files are served from the GitHub repository; verify them with
<code>sha256sum -c data/catalogue/SHA256SUMS</code> from the repository root. Column definitions:
<a href="methods.html#schema">data dictionary</a>.</p>
{lic_html}
<h2 id="catalogue">Catalogue release files</h2>
{cat_tbl}
<p>Any 100 kb end window not in <code>representatives.fa.gz</code> can be fetched from the public assembly and
checksum-verified with <code>python pipeline/fetch_end_window.py &lt;end_id&gt;</code>
(<a href="{BLOB}pipeline/fetch_end_window.py">script</a>). Provenance, input checksums and parameters:
<a href="{BLOB}data/catalogue/PROVENANCE.json">PROVENANCE.json</a>.</p>
<h2 id="figures">Atlas figures</h2>
{fig_tbl}
<h2 id="validation-files">Validation outputs</h2>
<p>Tan et al. 2026 cross-check (<code>data/validation/</code>) and query-tool validation (<code>tests/validation/</code>);
checksums computed when this site was built.</p>
{val_tbl}
<h2 id="site-data">Site data</h2>
<p>The arm pages read compact per-arm data files (<code>data/arms/&lt;arm&gt;.js</code>, derived from the release tables);
every arm page can export its filtered member table as TSV. The allele search index is
<a href="data/alleles_index.json">data/alleles_index.json</a>.</p>
<h2 id="analyses">Analysis results</h2>
<ul>{analyses}</ul>
<h2 id="versions">Versions</h2>
<p>Catalogue release {esc(VERSION)}; software: {esc(tools)}.</p>
"""
    write("downloads.html", page("downloads.html", "Downloads", body, active="downloads.html"))


def page_literature():
    p = rp("literature", "records.tsv")
    n = None
    top = ""
    try:
        r = pd.read_csv(p, sep="\t")
        n = len(r)
        if "relevance" in r.columns:
            t = r[pd.to_numeric(r.relevance, errors="coerce") >= 3].copy()
            t["doi_link"] = t.doi
            top = table_html(t.sort_values("year", ascending=False),
                             [("first_author", "First author"), ("year", "Year"), ("title", "Title"), ("venue", "Venue"),
                              ("doi_link", "DOI")], num=("year",),
                             fmt={"doi_link": lambda v, rr: f'<a href="https://doi.org/{esc(str(v))}">{esc(str(v))}</a>',
                                  "year": lambda v, rr: esc(str(int(v))) if pd.notna(v) else ""},
                             wrap_class="table-wrap free")
    except Exception:
        pass
    body = f"""
<h1>Literature</h1>
<p>An annotated review of prior and concurrent work on human subtelomere structure, exchange and evolution, and a
novelty assessment of the atlas findings against it:</p>
<ul>
<li><a href="{BLOB}literature/REVIEW.md">literature/REVIEW.md</a>: annotated literature review</li>
<li><a href="{BLOB}literature/NOVELTY.md">literature/NOVELTY.md</a>: novelty assessment of each finding</li>
<li><a href="{BLOB}literature/references.bib">literature/references.bib</a> and
<a href="{BLOB}literature/records.tsv">literature/records.tsv</a>{f' ({n} curated records)' if n else ''}</li>
</ul>
<h2 id="concurrent">Concurrent work</h2>
<p>Two 2026 preprints analyse the same class of data and should be read alongside the atlas: Tan et al.
(<a href="https://doi.org/{TAN}">doi:{TAN}</a>) and Guarracino et al. (<a href="https://doi.org/{GUARRACINO}">doi:{GUARRACINO}</a>).
The atlas is complementary to both; see the <a href="index.html#related">overview</a>.</p>
<h2 id="core">Core references (highest relevance)</h2>
{top or '<p class="muted">records.tsv not available at build time.</p>'}
<p class="muted">This list is generated from literature/records.tsv at build time.</p>
"""
    write("literature.html", page("literature.html", "Literature", body, active="literature.html"))


# ------------------------------------------------------------------ revision-2 pages and data
def load():
    global VERSION
    D = {}
    D["sa"] = pd.read_csv(os.path.join(CAT, "alleles_structural.tsv"), sep="\t")
    D["hg"] = pd.read_csv(os.path.join(CAT, "haplogroups.tsv"), sep="\t")
    D["fr"] = pd.read_csv(os.path.join(CAT, "allele_frequencies_by_superpop.tsv"), sep="\t")
    D["ends"] = pd.read_csv(os.path.join(CAT, "ends.tsv"), sep="\t", low_memory=False)
    D["mb"] = pd.read_csv(os.path.join(CAT, "mosaic_blocks.tsv"), sep="\t", low_memory=False,
                          usecols=["end_id", "level", "block", "label", "family", "source", "best_end", "sub_start",
                                   "sub_end"])
    D["pf"] = pd.read_csv(os.path.join(CAT, "paralog_families.tsv"), sep="\t")
    D["families"] = OrderedDict((r.family, r.members.split(",")) for r in D["pf"].itertuples())
    D["cat"] = load_json(os.path.join(CAT, "catalogue.json"), {})
    D["prov"] = load_json(os.path.join(CAT, "PROVENANCE.json"), {})
    D["manifest"] = pd.read_csv(os.path.join(CAT, "MANIFEST.tsv"), sep="\t")
    D["pal"] = pd.read_csv(rp("figures", "atlas", "PALETTE.tsv"), sep="\t")
    D["palqc"] = pd.read_csv(rp("figures", "atlas", "palette_qc.tsv"), sep="\t") \
        if os.path.exists(rp("figures", "atlas", "palette_qc.tsv")) else None
    D["captions"] = parse_captions(rp("figures", "atlas", "CAPTIONS.md"))
    D["rec"] = load_json(rp("analyses", "junctions", "results", "reciprocal_summary.json"), {})
    D["jun"] = load_json(rp("analyses", "junctions", "results", "junction_summary.json"), {})
    D["origin"] = load_json(rp("analyses", "allele_ages", "results", "origin_summary.json"), {})
    D["pop"] = load_json(rp("analyses", "allele_ages", "results", "pop_summary.json"), {})
    D["base"] = load_json(rp("analyses", "telomere_cis", "results", "baseline_stats.json"), {})
    VERSION = D["cat"].get("version") or D["prov"].get("release") or "v0.1"
    return D


def palette(D):
    pal = D["pal"]
    fam = D["families"]
    out = {"names": list(pal.chm13_end), "hex": list(pal.hex), "arm": list(pal.arm),
           "families": {f: ms for f, ms in fam.items()},
           "family_of": {m: f for f, ms in fam.items() for m in ms}}
    return out, {n: i for i, n in enumerate(out["names"])}


def fam_text(D, label):
    """'PF04' -> 'PF04 (3q, 19p)'; other labels unchanged."""
    ms = D["families"].get(label)
    return f"{label} ({', '.join(ms)})" if ms else label


def headline(D):
    e, sa, hg = D["ends"], D["sa"], D["hg"]
    H = {}
    H["n_ends"] = len(e)
    H["by_set"] = e.set.value_counts().to_dict()
    r2 = e[e.set == "HPRC_R2"]
    H["n_r2"] = len(r2)
    H["r2_samples"] = r2["sample"].nunique()
    H["n_sa"] = len(sa)
    H["sa_tiers"] = sa.tier.value_counts().to_dict()
    H["sa_nat"] = sa.nativity_class.value_counts().to_dict()
    H["n_hg"] = len(hg)
    H["hg_tiers"] = hg.tier.value_counts().to_dict()
    H["n_arms"] = sa.arm.nunique()
    H["r2_nat"] = r2.nativity_class.value_counts().to_dict()
    H["r2_S20_foreign"] = float((r2.S20_state == "foreign").mean())
    named = r2[r2.arm_source == "contig_name"]
    H["n_named"] = len(named)
    for w in ["S5", "S10", "S20"]:
        c = named[f"cis_{w}_state"].dropna() if f"cis_{w}_state" in named else pd.Series(dtype=str)
        H[f"cis_{w}"] = float((c == "foreign").mean()) if len(c) else float("nan")
        H["cis_n"] = int(len(c))
    s20q = sa[sa.allele_id == "SA-20q-01"]
    H["sa20q"] = s20q.iloc[0].to_dict() if len(s20q) else None
    return H


def representative_blocks(D):
    """Core blocks of each representative end with the best single CHM13 end of each block (taken from the raw block
    with the same label and the largest overlap; raw blocks carry best_end)."""
    mb = D["mb"]
    reps = set(D["sa"].representative_end_id)
    core = mb[(mb.level == "core") & mb.end_id.isin(reps)].sort_values(["end_id", "block"])
    raw = mb[(mb.level == "raw") & mb.end_id.isin(reps)]
    rawg = {k: d for k, d in raw.groupby("end_id")}
    out = {}
    for e, d in core.groupby("end_id"):
        blocks = []
        rr = rawg.get(e)
        for b in d.itertuples():
            be = None
            if b.label != "none" and rr is not None:
                cand = rr[(rr.label == b.label) & (rr.sub_end > b.sub_start) & (rr.sub_start < b.sub_end)]
                if len(cand):
                    ov = cand.apply(lambda r: min(r.sub_end, b.sub_end) - max(r.sub_start, b.sub_start), axis=1)
                    w = ov.groupby(cand.best_end).sum()
                    be = w.idxmax()
            blocks.append((int(b.sub_start), int(b.sub_end), b.label, be if isinstance(be, str) else None))
        out[e] = blocks
    return out


def build_data(D, pal, pidx):
    sa, hg, fr, e, mb = D["sa"], D["hg"], D["fr"], D["ends"], D["mb"]
    freq = fr[fr.superpopulation.isin(POPS)].pivot_table(index="allele_id", columns="superpopulation",
                                                         values="frequency")
    raw = mb[mb.level == "raw"]
    rawg = {k: d.sort_values("block") for k, d in raw.groupby("end_id")}
    repb = representative_blocks(D)
    none_i = pidx["none"]

    group = {"acroP": ["13p", "14p", "15p", "21p", "22p"], "PAR1": ["Xp", "Yp"], "PAR2": ["Xq", "Yq"]}

    def idx_of(best_end, label, arm):
        """home blocks: own arm (pooled groups: best end if a member); others: best single end."""
        if label == "none":
            return none_i
        if label == "home":
            be = best_end[3:] if isinstance(best_end, str) else ""
            if arm in group:
                return pidx.get("chr" + (be if be in group[arm] else POOLED_HOME[arm]), none_i)
            return pidx.get("chr" + arm, none_i)
        if isinstance(best_end, str) and best_end in pidx:
            return pidx[best_end]
        if isinstance(label, str) and label in D["families"]:
            return pidx.get("chr" + D["families"][label][0], none_i)
        if isinstance(label, str):
            return pidx.get("chr" + label, none_i)
        return none_i

    def clamp(v):
        return int(max(0, min(100000, v)))

    def ints(v):
        return [int(x) for x in str(v).split(",")] if pd.notna(v) else []

    arms_data = {}
    for arm, g in sa.groupby("arm"):
        alleles = []
        for r in g.itertuples():
            blk = [[clamp(s), clamp(en), idx_of(be, lab, arm), lab,
                    (be[3:] if isinstance(be, str) and lab not in ("home", "none") else "")]
                   for s, en, lab, be in repb.get(r.representative_end_id, [])]
            fq = {p: (round(float(freq.loc[r.allele_id, p]), 4) if r.allele_id in freq.index and p in freq.columns
                      and pd.notna(freq.loc[r.allele_id, p]) else None) for p in POPS}
            alleles.append({
                "id": r.allele_id, "st": r.structure, "form": r.form if isinstance(r.form, str) else "",
                "nat": bool(r.is_native), "nc": r.nativity_class if isinstance(r.nativity_class, str) else "",
                "tier": r.tier, "ev": r.evidence if isinstance(r.evidence, str) else "",
                "ne": int(r.n_ends), "nr2": int(r.n_r2_ends), "ns": int(r.n_samples),
                "f": round(float(r.r2_freq), 4) if pd.notna(r.r2_freq) else None, "fq": fq,
                "sw": ints(r.switches_median), "swmin": ints(r.switches_min), "swmax": ints(r.switches_max),
                "blk": blk, "rep": r.representative_end_id,
                "hgs": r.haplogroups_top if isinstance(r.haplogroups_top, str) else "",
                "evc": {"vkc": int(r.n_vk_concordant), "vkd": int(r.n_vk_discordant), "lay": int(r.n_layout_supported),
                        "rd": int(r.n_hg002_reads_supported), "ped": int(r.n_pedigree_intact),
                        "t2t": int(r.n_chm13 + r.n_hg002_q100)},
            })
        hgs = [{"id": r.haplogroup_id, "tier": r.tier, "ne": int(r.n_ends), "ns": int(r.n_samples),
                "f": round(float(r.r2_freq), 4) if pd.notna(r.r2_freq) else None}
               for r in hg[hg.arm == arm].itertuples()]
        mcols = ["end_id", "sample", "superpopulation", "set", "assembly_accession", "structural_allele_id",
                 "haplogroup_id", "nativity_class", "S20_state", "S20_best_end", "form", "sa_assignment", "flags"]
        m = e[e.catalogue_arm == arm]
        members = []
        for row in m[mcols].itertuples(index=False):
            vals = ["" if (not isinstance(v, str) and pd.isna(v)) else v for v in row]
            rb = rawg.get(row[0])
            blocks = [[clamp(b.sub_start), clamp(b.sub_end), idx_of(b.best_end, b.label, arm)] for b in rb.itertuples()] \
                if rb is not None else []
            members.append(vals + [blocks])
        arms_data[arm] = {"arm": arm, "version": VERSION,
                          "n_r2_ends": int(D["cat"].get("arms", {}).get(arm, {}).get("n_r2_ends", 0)),
                          "alleles": alleles, "haplogroups": hgs, "mcols": mcols, "members": members}
    cols = ["id", "type", "arm", "structure", "form", "nativity_class", "families", "tier", "evidence", "n_ends",
            "n_samples", "r2_freq"] + POPS
    rows = []
    famre = re.compile(r"PF\d+")
    for r in sa.itertuples():
        fq = [round(float(freq.loc[r.allele_id, p]), 4) if r.allele_id in freq.index and pd.notna(freq.loc[r.allele_id, p])
              else None for p in POPS]
        rows.append([r.allele_id, "structural allele", r.arm, r.structure, r.form if isinstance(r.form, str) else "",
                     r.nativity_class if isinstance(r.nativity_class, str) else "",
                     " ".join(sorted(set(famre.findall(str(r.structure))))), r.tier,
                     r.evidence if isinstance(r.evidence, str) else "", int(r.n_ends), int(r.n_samples),
                     round(float(r.r2_freq), 4) if pd.notna(r.r2_freq) else None] + fq)
    for r in hg.itertuples():
        fq = [round(float(freq.loc[r.haplogroup_id, p]), 4) if r.haplogroup_id in freq.index and
              pd.notna(freq.loc[r.haplogroup_id, p]) else None for p in POPS]
        top = r.structural_alleles_top if isinstance(r.structural_alleles_top, str) else ""
        rows.append([r.haplogroup_id, "haplogroup", r.arm, top, "", "", "", r.tier,
                     r.evidence if isinstance(r.evidence, str) else "", int(r.n_ends), int(r.n_samples),
                     round(float(r.r2_freq), 4) if pd.notna(r.r2_freq) else None] + fq)
    index = {"version": VERSION, "cols": cols, "rows": rows, "families": {f: ms for f, ms in D["families"].items()}}
    return arms_data, index, repb


def definitions_html(D):
    fams = "; ".join(f"<strong>{esc(f)}</strong> = {esc(', '.join(ms))}" for f, ms in D["families"].items())
    pf = D["pf"]
    tau = pf.tau.iloc[0] if "tau" in pf else None
    margin = pf.margin.iloc[0] if "margin" in pf else None
    return f"""
<div class="note" id="definitions-box"><p><strong>How blocks are labelled (catalogue revision 2).</strong> At every 100 bp
point of an end's subtelomere, the identity of the surrounding {int(pf.window_bp.iloc[0]) // 1000 if 'window_bp' in pf else 5} kb
to each of the 48 CHM13 v2.0 ends is read from one alignment of the whole window.</p>
<ul>
<li><strong>Home rule.</strong> A point is <em>home</em> if a CHM13 copy of the end's own arm (or of an arm in its pooled
group) is within {margin * 100 if margin is not None else 0.5:g} percentage points of the best identity, so home wins ties.</li>
<li><strong>Paralog families.</strong> Some CHM13 ends share long, near-identical duplicated blocks, so the single best end
can flip between them. Ends that are tied with each other at &ge; {tau if tau is not None else 0.75:g} of the points where
the less often tied one is tied form a family; a non-home block is labelled with the family of its best end, otherwise
with that end's arm. Families: {fams}. The best single CHM13 end of each block is kept as <code>best_end</code>.</li>
<li><strong>none.</strong> No CHM13 end within 0.80 identity: sequence without a close CHM13 counterpart (neither home nor
foreign evidence).</li>
<li><strong>Nativity class</strong> of a core structure: <em>native</em> = only home (and none) blocks; <em>within
family</em> = every non-home block is in a family that contains the end's own arm; <em>foreign</em> = at least one block
from another family or arm.</li>
</ul></div>"""


def page_index(D, H, AG, repb):
    pre = ""
    caps = D["captions"]
    rec, jun, org, pop, base = D["rec"], D["jun"], D["origin"], D["pop"], D["base"]
    s = H["sa20q"]
    be20 = ""
    if s:
        b = repb.get(s["representative_end_id"], [])
        if b and b[0][3]:
            be20 = b[0][3][3:]
    tan_title, gua_title = lit_title(TAN), lit_title(GUARRACINO)
    bs = H["by_set"]
    noncar = rec.get("pooled", {}).get("frac", {}).get("noncarriers", {}).get("reciprocal")
    chi_p = rec.get("pooled", {}).get("chi2_donor_status", {}).get("p")
    ape = org.get("ape", {})
    age_txt = ""
    if AG:
        age_txt = (f" Median-based clock ages of {AG['n_in_range']} of {AG['n_clustered']} clustered forms fall between "
                   f"0.2 and 1.5 Myr (median {AG['median']:.0f} kyr); these ages depend on the mutation-rate and "
                   f"genealogy model.")
    nat = H["r2_nat"]
    st20 = ""
    if s:
        fam = str(s["structure"]).split("|")[0]
        st20 = (f'<a href="arms/20q.html#SA-20q-01">SA-20q-01</a> (structure <code>{esc(str(s["structure"]))}</code>: '
                f'20q carrying a terminal block of paralog family {esc(fam_text(D, fam))}'
                f'{", best single end " + esc(be20) if be20 else ""}; median switch '
                f'{float(s["switches_median"]) / 1000:.1f} kb from the telomere boundary) is carried by {s["n_samples"]} '
                f'samples (Release 2 frequency {s["r2_freq"]}; tier {s["tier"]}).')
    body = f"""
<h1>Human subtelomere atlas</h1>
<p class="lead">A validated, queryable catalogue of the structural alleles at human chromosome ends. It describes the
terminal 100 kb next to every telomere in {fnum(bs.get('HPRC_R2', 0))} telomeric contig ends from HPRC Release 2
({H['r2_samples']} samples), the {bs.get('HG002_Q100', 0)} ends of the HG002 Q100 assembly and the {bs.get('CHM13', 0)}
ends of T2T-CHM13 v2.0, as mosaics of blocks labelled home, by paralog family or by the CHM13 arm they resemble most.</p>
<p>Concurrent work on the same question: Tan et al. 2026 (<a href="https://doi.org/{TAN}">doi:{TAN}</a>) and Guarracino
et al. 2026 (<a href="https://doi.org/{GUARRACINO}">doi:{GUARRACINO}</a>); this atlas is complementary to both
(<a href="#related">details</a>).</p>
{definitions_html(D)}
<div class="note"><strong>Terminology.</strong> A <em>switch</em> (junction) is a change in which CHM13 end the sequence
most resembles; it lies inside shared, paralogous sequence and is not necessarily a structural breakpoint. All labels are
relative to the CHM13 reference.</div>
<div class="grid">
<div class="card"><span class="stat">{fnum(H['n_ends'])}</span>chromosome ends ({fnum(bs.get('HPRC_R2', 0))} HPRC
Release 2, {bs.get('CHM13', 0)} CHM13, {bs.get('HG002_Q100', 0)} HG002 Q100, {bs.get('pedigree_aux', 0)} pedigree ends used
only for validation)</div>
<div class="card"><span class="stat">{fnum(H['n_sa'])}</span>structural alleles on {H['n_arms']} catalogue arms
(tier A {H['sa_tiers'].get('A', 0)}, B {H['sa_tiers'].get('B', 0)}, C {H['sa_tiers'].get('C', 0)}; native
{H['sa_nat'].get('native', 0)}, within family {H['sa_nat'].get('within_family', 0)}, foreign
{H['sa_nat'].get('foreign', 0)})</div>
<div class="card"><span class="stat">{fnum(H['n_hg'])}</span>sequence haplogroups (k-mer ANI &ge;
{D['cat'].get('haplogroup_threshold_ani', 0.995)} over the first 80 kb; tier A {H['hg_tiers'].get('A', 0)},
B {H['hg_tiers'].get('B', 0)}, C {H['hg_tiers'].get('C', 0)})</div>
<div class="card"><span class="stat">{fnum(nat.get('foreign', 0))}</span>Release 2 ends with a foreign core structure
(native {fnum(nat.get('native', 0))}, within family {fnum(nat.get('within_family', 0))})</div>
</div>
{figure('figures/atlas/F1_genome_overview.png', caption_title(caps.get('F1_genome_overview')) +
        ': arms on a ring, chords from the best CHM13 end of terminal blocks to the recipient arm, and nativity classes per arm.',
        md_inline(caps.get('F1_genome_overview', '')), pre)}
<h2 id="findings">Key findings</h2>
<p class="muted">Finding 1 uses catalogue revision 2. Findings 2-4 come from the three analyses, which used the earlier
per-window best-end form calls (e.g. <code>20q&gt;3q</code>); their counts are not directly comparable with revision-2
S20 or nativity counts (<a href="findings/index.html#definitions">definitions and concordance table</a>).</p>
<ol>
<li><strong>Non-native terminal sequence is common.</strong> Of {fnum(H['n_r2'])} Release 2 ends,
{fnum(nat.get('foreign', 0))} ({nat.get('foreign', 0) / max(1, H['n_r2']):.1%}) have a foreign and
{fnum(nat.get('within_family', 0))} ({nat.get('within_family', 0) / max(1, H['n_r2']):.1%}) a within-family core
structure; the first 20 kb are mostly foreign in {H['r2_S20_foreign']:.1%} (S20 state). The fraction depends on the
window and on the labelling rule: the earlier per-window best-end calls gave {H.get('cis_S5', 0):.1%},
{H.get('cis_S10', 0):.1%} and {H.get('cis_S20', 0):.1%} for the first 5, 10 and 20 kb of
{fnum(H.get('cis_n', 0))} chromosome-named tips.</li>
<li><strong>Common forms are old, recurrent structural alleles with independent validation.</strong> For example, {st20}
In the analyses, unrelated people share base-identical switch points:
{pop.get('junction_clusters_ge10_shared_AFRcont_OOA', '')} of {pop.get('junction_clusters_ge10_samples', '')} junction
clusters with &ge; 10 samples contain both continental African and non-African samples.{age_txt} None of the
{ape.get('junction_clusters_tested', '')} recurrent junctions tested was found in the T2T great-ape genomes
({ape.get('candidate_shared_junctions', '')} candidates after all {ape.get('carrier_probe_spanning_hits', '')} probe hits
were diagnosed as shared duplicons); absence is weak evidence, because orthology of ape chromosome ends could be assigned
only in part.</li>
<li><strong>The transfer is one-way and homology-mediated.</strong> The donor arm's end in carrier haplotypes looks like
that of non-carriers (reciprocal configuration {rec.get('reciprocal_frac', 0):.2%} vs {noncar or 0:.2%}; chi-square
p = {chi_p if chi_p is not None else float('nan'):.2f}), and {jun.get('junction_class', {}).get('homology_spanning', '')} of
{jun.get('analysed', '')} junction clusters
({jun.get('junction_class', {}).get('homology_spanning', 0) / max(1, jun.get('analysed', 1)):.0%}) lie within
donor-recipient homology. Interstitial telomeric sequence is not a general feature of junctions.</li>
<li><strong>The terminal allele does not set telomere length in cis.</strong> Arm-specific telomere lengths in the assemblies
reproduce published nanopore measurements (Spearman &rho; = {base.get('mean_rel_tl', {}).get('spearman_rho', float('nan')):.2f},
{base.get('n_arms_compared', '')} arms), but foreign and own terminal alleles in the same individual do not differ in
telomere length; an earlier positive association was an artefact of a telomere-length-dependent classification window.</li>
</ol>
<p><a href="findings/index.html">Findings in detail</a> &middot; <a href="validation.html">Validation</a> &middot;
<a href="arms/index.html">Browse arms</a> &middot; <a href="alleles.html">Search alleles</a></p>
<h2 id="related">Relation to concurrent work</h2>
<p>This atlas is not the first population catalogue of human subtelomeres. It complements two concurrent preprints:</p>
<ul>
<li>Tan et al. 2026, <em>{esc(tan_title)}</em>, bioRxiv, <a href="https://doi.org/{TAN}">doi:{TAN}</a>: complete
subtelomeres from 860 haplotypes over 500 kb, organised into paralogy-block haplotypes; sequences released under CC BY 4.0
(<a href="https://doi.org/{TAN_ZENODO}">doi:{TAN_ZENODO}</a>).</li>
<li>Guarracino et al. 2026, <em>{esc(gua_title)}</em>, bioRxiv, <a href="https://doi.org/{GUARRACINO}">doi:{GUARRACINO}</a>:
all-vs-all comparison of chromosome ends in near-complete assemblies, grouping pseudo-homologous regions into sequence
communities and reporting candidate subtelomere recombination in a pedigree.</li>
</ul>
<p>The atlas contributes per-end membership for every HPRC Release 2 telomeric tip, telomere-boundary-anchored structural
alleles labelled by paralog family with nativity classes, sequence haplogroups with explicit evidence tiers (Verkko, read
and pedigree support), boundary-shift stability tests, representative sequences, frequencies by superpopulation, a query
tool, and hypothesis tests of mechanism, age and telomere-length effects. A same-haplotype cross-check against the
independent Tan et al. sequences is on the <a href="validation.html#tan2026">validation page</a>; see the
<a href="literature.html">literature review</a> for a fuller comparison.</p>
<h2 id="cite">How to cite</h2>
<p>Human subtelomere atlas, allele catalogue {esc(VERSION)} (2026). {REPO_URL}. Please cite the concurrent studies above
where relevant, and the underlying data: HPRC (Liao et al. 2023,
<a href="https://doi.org/10.1038/s41586-023-05896-x">doi:10.1038/s41586-023-05896-x</a>; Release 2 data at
<a href="https://humanpangenome.org/">humanpangenome.org</a>), T2T-CHM13 (Nurk et al. 2022,
<a href="https://doi.org/10.1126/science.abj6987">doi:10.1126/science.abj6987</a>), and the HG002 Q100 assembly. Allele
identifiers are specific to release {esc(VERSION)} (revision 2); the mapping from earlier identifiers is in
<a href="{BLOB}data/catalogue/id_history.tsv">id_history.tsv</a>.</p>
"""
    write("index.html", page("index.html", "Overview", body, active="index.html"))


def page_arm(D, arm, data, AS, pal, repb):
    pre = "../"
    st = AS[AS.arm == arm].iloc[0]
    sa = D["sa"][D["sa"].arm == arm].sort_values("n_ends", ascending=False).copy()
    hg = D["hg"][D["hg"].arm == arm].sort_values("n_ends", ascending=False).copy()
    pooled = {"acroP": "13p, 14p, 15p, 21p and 22p", "PAR1": "Xp and Yp", "PAR2": "Xq and Yq"}.get(arm)
    pooled_txt = f" This group pools the ends of {pooled}." if pooled else ""
    fr = D["fr"]
    frp = fr[fr.superpopulation.isin(POPS)].pivot_table(index="allele_id", columns="superpopulation", values="frequency")
    for df, key in [(sa, "allele_id"), (hg, "haplogroup_id")]:
        for p in POPS:
            df[p] = [frp.loc[i, p] if i in frp.index and p in frp.columns else None for i in df[key]]
    sa["best_ends"] = ["|".join(("home" if b[2] == "home" else "none" if b[2] == "none" else (b[3][3:] if b[3] else "-"))
                                for b in repb.get(r, [])) for r in sa.representative_end_id]
    f3 = lambda v, r: "" if pd.isna(v) else f"{v:.3f}"
    f2 = lambda v, r: "" if pd.isna(v) else f"{v:.2f}"
    num_cols = ("n_ends", "n_samples", "r2_freq") + tuple(POPS)
    famre = re.compile(r"PF\d+")
    sa_tbl = table_html(sa, [("allele_id", "Allele"), ("structure", "Structure"), ("best_ends", "Best end per block"),
                             ("nativity_class", "Nativity"), ("switches_median", "Median switch(es), bp"),
                             ("tier", "Tier"), ("evidence", "Evidence"), ("n_ends", "Ends"), ("n_samples", "Samples"),
                             ("r2_freq", "R2 freq.")] + [(p, p) for p in POPS] + [("haplogroups_top", "Haplogroups")],
                        num=num_cols, id_col="allele_id",
                        fmt={"tier": lambda v, r: tier_badge(v), "r2_freq": f3, **{p: f2 for p in POPS},
                             "structure": lambda v, r: famre.sub(lambda m: f'<abbr title="{esc(fam_text(D, m.group(0)))}">'
                                                                      f'{m.group(0)}</abbr>', esc(str(v))),
                             "nativity_class": lambda v, r: esc(str(v).replace("_", " ")) if isinstance(v, str) else ""})
    hg_tbl = table_html(hg, [("haplogroup_id", "Haplogroup"), ("tier", "Tier"), ("evidence", "Evidence"),
                             ("n_ends", "Ends"), ("n_samples", "Samples"), ("r2_freq", "R2 freq.")] +
                        [(p, p) for p in POPS] + [("medoid_end_id", "Medoid end"),
                                                  ("structural_alleles_top", "Structural alleles"),
                                                  ("mean_within_ani", "Mean within ANI")],
                        num=num_cols + ("mean_within_ani",), id_col="haplogroup_id",
                        fmt={"tier": lambda v, r: tier_badge(v), "r2_freq": f3, **{p: f2 for p in POPS},
                             "mean_within_ani": lambda v, r: "" if pd.isna(v) else f"{v:.4f}"})
    # families relevant to this arm
    own = [POOLED_HOME.get(arm, arm)] if arm not in ("acroP", "PAR1", "PAR2") else \
        {"acroP": ["13p", "14p", "15p", "21p", "22p"], "PAR1": ["Xp", "Yp"], "PAR2": ["Xq", "Yq"]}[arm]
    own_f = sorted({f for f, ms in D["families"].items() if any(o in ms for o in own)})
    used_f = sorted(set(famre.findall(" ".join(sa.structure.astype(str)))))
    fam_line = []
    if own_f:
        fam_line.append("this arm belongs to " + ", ".join(esc(fam_text(D, f)) for f in own_f) +
                        " (blocks of that family on this arm are <em>within family</em>)")
    other = [f for f in used_f if f not in own_f]
    if other:
        fam_line.append("its alleles contain blocks of " + ", ".join(esc(fam_text(D, f)) for f in other))
    fam_txt = ("Paralog families: " + "; ".join(fam_line) + ".") if fam_line else \
        "No paralog-family blocks occur in this arm's alleles."
    nat = D["ends"][(D["ends"].set == "HPRC_R2") & (D["ends"].catalogue_arm == arm)].nativity_class.value_counts()
    arms = list(AS.arm)
    navs = "".join(f'<a href="{a}.html"{" aria-current=page" if a == arm else ""}>{a}</a>' for a in arms)
    sets = sorted(set(m[3] for m in data["members"]))
    body = f"""
<nav class="arm-nav" aria-label="Arms">{navs}</nav>
<h1>Arm {esc(arm)}</h1>
<p>{fnum(st.r2_ends)} Release 2 ends from {st.r2_samples} samples (native {int(nat.get('native', 0))}, within family
{int(nat.get('within_family', 0))}, foreign {int(nat.get('foreign', 0))}); {st.n_sa} structural alleles (tier A {st.A},
B {st.B}, C {st.C}) and {st.n_hg} sequence haplogroups.{pooled_txt}</p>
<p>{fam_txt} <a href="../methods.html#labels">How blocks are labelled</a>.</p>
<h2 id="mosaic">Structural alleles</h2>
<p>One row per structural allele: the representative end's core blocks (&ge; 10 kb) from the telomere boundary (left,
dark bar) inward. Home blocks carry the arm's own colour; other blocks are coloured by their best single CHM13 v2.0 end;
members of a paralog family share a hue and differ in lightness; dashed outline = block labelled by paralog family. Triangles mark median core switch positions; the
thin bar under a row spans the members' minimum to maximum switch. Hover over a row, or focus it with the keyboard, for
details (structure, families, best ends, nativity class, frequencies, evidence); click or press Enter to filter the member
table to that allele.</p>
<div class="controls" role="group" aria-label="Allele filters">
<label class="inline"><input type="checkbox" id="f-tier-A" checked> tier A</label>
<label class="inline"><input type="checkbox" id="f-tier-B" checked> tier B</label>
<label class="inline"><input type="checkbox" id="f-tier-C" checked> tier C</label>
<label>Nativity <select id="f-nat"><option value="">all</option><option value="native">native</option>
<option value="within_family">within family</option><option value="foreign">foreign</option></select></label>
<label>Minimum samples <input type="number" id="f-minsamples" min="1" value="1" style="width:5rem"></label>
<label>Order <select id="f-order"><option value="members">by member count</option><option value="switch">by first switch</option>
<option value="id">by ID</option></select></label>
<button type="button" class="btn" id="toggle-members-mosaic" aria-expanded="false" aria-controls="member-mosaic">Show member-end mosaic</button>
<span class="count" id="allele-count" aria-live="polite"></span>
</div>
<div class="mosaic-wrap" id="allele-mosaic"><noscript><p>The interactive mosaic needs JavaScript; the allele table below has
the same information.</p></noscript></div>
<div class="legend" id="mosaic-legend" aria-label="Colour key"></div>
<div class="mosaic-wrap" id="member-mosaic" hidden></div>
<p class="muted">Labels are relative to CHM13: a colour change is a change in which CHM13 end is most similar, not
necessarily a structural breakpoint. Palette: <a href="../methods.html#palette">methods</a>.</p>
<h3 id="allele-table">Allele table</h3>
{sa_tbl}
<h2 id="haplogroups">Sequence haplogroups</h2>
<p>Clusters of ends with k-mer ANI &ge; {D['cat'].get('haplogroup_threshold_ani', 0.995)} over the first 80 kb of
subtelomere (UPGMA).</p>
{hg_tbl}
<h2 id="members">Member ends</h2>
<p id="selection-status" class="muted" aria-live="polite"></p>
<div class="controls" role="group" aria-label="Member filters">
<label>Search <input type="search" id="m-search" placeholder="end, sample, accession, flag"></label>
<label>Structural allele <select id="m-allele"><option value="">all</option></select></label>
<label>Haplogroup <select id="m-hg"><option value="">all</option></select></label>
<label>Nativity <select id="m-nat"><option value="">all</option><option value="native">native</option>
<option value="within_family">within family</option><option value="foreign">foreign</option></select></label>
<label>Superpopulation <select id="m-pop"><option value="">all</option>{''.join(f'<option>{p}</option>' for p in POPS)}<option value="NA">none (CHM13, pedigree)</option></select></label>
<label>Set <select id="m-set"><option value="">all</option>{''.join(f'<option>{esc(x)}</option>' for x in sets)}</select></label>
<button type="button" class="btn" id="m-download">Download shown rows (TSV)</button>
<span class="count" id="member-count" aria-live="polite"></span>
</div>
<div class="table-wrap" id="member-table"><noscript><p>The member table needs JavaScript. All ends are in
<a href="{BLOB}data/catalogue/ends.tsv">data/catalogue/ends.tsv</a>.</p></noscript></div>
<p class="muted">Nativity is the class of the end's own core structure; S20 state and S20 best end describe the first
20 kb of subtelomere (majority label; best single CHM13 end). Pedigree ends (set <code>pedigree_aux</code>) are validation
evidence only and are not counted in allele sizes or frequencies. Full tables:
<a href="{BLOB}data/catalogue/ends.tsv">ends.tsv</a>, <a href="{BLOB}data/catalogue/alleles_structural.tsv">alleles_structural.tsv</a>,
<a href="{BLOB}data/catalogue/haplogroups.tsv">haplogroups.tsv</a>, <a href="{BLOB}data/catalogue/mosaic_blocks.tsv">mosaic_blocks.tsv</a>;
representative sequences in <a href="{BLOB}data/catalogue/representatives.fa.gz">representatives.fa.gz</a>; any other
window with <code>python pipeline/fetch_end_window.py &lt;end_id&gt;</code>; new sequences with the
<a href="../query.html">query tool</a>.</p>
"""
    write(f"arms/{arm}.html", page(f"arms/{arm}.html", f"Arm {arm}", body, active="arms/index.html",
                                   scripts=("data/palette.js", f"data/arms/{arm}.js", "assets/js/arm.js"),
                                   description=f"Structural alleles and haplogroups of chromosome end {arm}."))


def page_alleles(D, index):
    fam_opts = "".join(f'<option value="{esc(f)}">{esc(fam_text(D, f))}</option>' for f in D["families"])
    body = f"""
<h1>Alleles</h1>
<p>All {fnum(len(D['sa']))} structural alleles and {fnum(len(D['hg']))} sequence haplogroups of catalogue
{esc(VERSION)}. Search matches every word against ID, arm, structure, form, family and evidence (for example
<code>20q PF04</code>, <code>16q|home</code> or <code>pedigree</code>). Click an ID to open its arm page. Frequencies are
among Release 2 ends of the arm in each superpopulation. Structures read telomere &rarr; centromere: <code>home</code> =
own arm, <code>PFnn</code> = paralog family, an arm name = that CHM13 arm, <code>none</code> = no CHM13 end within 0.80
identity (<a href="methods.html#labels">definitions</a>).</p>
<div class="controls" role="search" aria-label="Allele search">
<label>Search <input type="search" id="q" placeholder="e.g. 20q PF04" style="min-width:14rem"></label>
<label>Type <select id="q-type"><option value="">all</option><option>structural allele</option><option>haplogroup</option></select></label>
<label>Tier <select id="q-tier"><option value="">all</option><option>A</option><option>B</option><option>C</option></select></label>
<label>Arm <select id="q-arm"><option value="">all</option></select></label>
<label>Nativity <select id="q-nat"><option value="">all</option><option value="native">native</option>
<option value="within_family">within family</option><option value="foreign">foreign</option></select></label>
<label>Paralog family <select id="q-fam"><option value="">any</option>{fam_opts}</select></label>
<label>Minimum samples <input type="number" id="q-minsamples" min="0" value="0" style="width:5rem"></label>
<button type="button" class="btn" id="q-download">Download shown rows (TSV)</button>
<span class="count" id="allele-count" aria-live="polite"></span>
</div>
<div class="table-wrap" id="allele-table"><noscript><p>The search table needs JavaScript. The same data are in
<a href="data/alleles_index.json">data/alleles_index.json</a> and in
<a href="{BLOB}data/catalogue/alleles_structural.tsv">alleles_structural.tsv</a> /
<a href="{BLOB}data/catalogue/haplogroups.tsv">haplogroups.tsv</a>.</p></noscript></div>
<p class="muted">Search index: <a href="data/alleles_index.json">data/alleles_index.json</a> (columns
{esc(', '.join(index['cols']))}). Nativity and family filters apply to structural alleles; haplogroups are sequence
clusters and have no structure of their own (their main structural alleles are listed).</p>
"""
    write("alleles.html", page("alleles.html", "Alleles", body, active="alleles.html",
                               scripts=("data/alleles_index.js", "assets/js/alleles.js")))


def concordance(D, n_forms=20):
    """Analysis forms (ends.tsv prior_form: the earlier per-window best-end call used by the three analyses) vs
    catalogue revision 2 (nativity class, form, structural allele) on the same Release 2 ends."""
    e = D["ends"]
    r2 = e[(e.set == "HPRC_R2") & e.prior_form.notna()]
    sa = D["sa"].set_index("allele_id")
    top = r2.prior_form.value_counts()
    named = ["20q>3q", "3q>11p", "6p>11p", "1p>11p", "8p>1p", "16p>9q"]
    forms = list(top.index[:n_forms]) + [f for f in named if f in top.index and f not in list(top.index[:n_forms])]
    rows = []
    for f in forms:
        d = r2[r2.prior_form == f]
        n = len(d)
        nc = d.nativity_class.fillna("unassigned").value_counts()
        fm = d.form.fillna("(none: S20 not foreign)").value_counts()
        al = d.structural_allele_id.dropna().value_counts()
        a0 = al.index[0] if len(al) else ""
        rows.append({"analysis_form": f, "ends": n, "headline_extra": f not in list(top.index[:n_forms]),
                     "native": nc.get("native", 0) / n, "within_family": nc.get("within_family", 0) / n,
                     "foreign": nc.get("foreign", 0) / n,
                     "catalogue_form": f"{fm.index[0]} ({fm.iloc[0] / n:.0%})",
                     "main_allele": (f"{a0} {sa.loc[a0, 'structure']} ({al.iloc[0] / n:.0%})" if a0 else ""),
                     "n_alleles": int(len(al))})
    return pd.DataFrame(rows)


def page_findings(D):
    pre = "../"
    caps = D["captions"]
    conc = concordance(D)
    pct = lambda v, r: f"{v:.0%}"
    conc_tbl = table_html(conc, [("analysis_form", "Analysis form (earlier call)"), ("ends", "Release 2 ends"),
                                 ("native", "native"), ("within_family", "within family"), ("foreign", "foreign"),
                                 ("catalogue_form", "Main revision-2 form"), ("main_allele", "Main structural allele"),
                                 ("n_alleles", "Alleles")],
                          num=("ends", "native", "within_family", "foreign", "n_alleles"),
                          fmt={"native": pct, "within_family": pct, "foreign": pct}, wrap_class="table-wrap free")
    defs = f"""
<h2 id="definitions">Definitions: analysis forms vs catalogue revision 2</h2>
<p class="warn">The three analyses (telomere_cis, junctions, allele_ages) used the earlier per-window best-end form calls
(a terminal block called by the single best CHM13 end of a 20 kb window, e.g. <code>20q&gt;3q</code>; column
<code>prior_form</code> of ends.tsv). Catalogue revision 2 labels blocks from per-base identity tracks, applies the home
rule and groups near-identical CHM13 ends into paralog families (e.g. <code>20q&gt;PF04</code>, PF04 = 3q/19p), so
analysis counts are not directly comparable with revision-2 S20, form or nativity counts. The table maps the
{int((~conc.headline_extra).sum())} most frequent analysis forms (plus {int(conc.headline_extra.sum())} further headline
forms of the analyses) onto revision 2 for the same Release 2
ends: share of ends in each nativity class, the most common revision-2 form, and the most common structural allele.</p>
{conc_tbl}
<p class="muted">Computed at build time from data/catalogue/ends.tsv (<code>prior_form</code>, <code>nativity_class</code>,
<code>form</code>, <code>structural_allele_id</code>) and alleles_structural.tsv. "(none: S20 not foreign)" = the
revision-2 S20 window is not mostly foreign (for example a terminal block shorter than 20 kb or a home tie).</p>
"""
    cards = "".join(f'<div class="card"><h3><a href="{k}.html">{esc(v["short"])}</a></h3><p>{esc(v["title"])}.</p>'
                    f'<p><a href="{BLOB}{v["dir"]}/README.md">README on GitHub</a></p></div>'
                    for k, v in ANALYSES.items())
    figs = []
    for key, alt in [("F1_genome_overview", "Genome overview ring and nativity class per arm"),
                     ("F2_allele_mosaics", "Allele mosaics of six arms"),
                     ("F3_population_frequencies", "Heatmap of common non-native allele frequencies by superpopulation"),
                     ("F4_validation", "Evidence tiers and evidence-line combinations")]:
        figs.append(f'<h3 id="{key}">{esc(caption_title(caps.get(key)))}</h3>' +
                    figure(f"figures/atlas/{key}.png", f"{caption_title(caps.get(key))}: {alt}.",
                           md_inline(caps.get(key, "")), pre))
    body = f"""
<h1>Findings</h1>
<p>Three hypothesis-driven analyses test what the common terminal forms are, how they arose and whether they affect
telomere length. Each page summarises the verdicts; full methods, statistics and caveats are in the analysis READMEs.</p>
<div class="grid">{cards}</div>
{defs}
<h2 id="atlas-figures">Atlas figures (catalogue revision 2)</h2>
<p>Generated by <a href="{BLOB}pipeline/figures_atlas.py"><code>pipeline/figures_atlas.py</code></a> from the release
tables; PNG (300 dpi), PDF and SVG in <a href="{TREE}figures/atlas">figures/atlas/</a>; captions in
<a href="{BLOB}figures/atlas/CAPTIONS.md">CAPTIONS.md</a>.</p>
{''.join(figs)}
"""
    write("findings/index.html", page("findings/index.html", "Findings", body, active="findings/index.html"))
    note = ('<p class="warn">Definitions: this analysis used the earlier per-window best-end form calls (e.g. '
            '<code>20q&gt;3q</code>); its counts are not directly comparable with catalogue revision-2 S20, form or '
            'nativity counts. See the <a href="index.html#definitions">concordance table</a>.</p>')
    for key, v in ANALYSES.items():
        d = v["dir"]
        readme = f"{d}/README.md"
        if key == "telomere-cis":
            summ = md_section(readme, r"Summary of verdicts", pre, hshift=0)
        elif key == "junctions":
            summ = "<h2 id=\"verdicts\">Verdicts</h2>" + md_paragraphs_starting(readme, r"^\*\*Verdict", pre)
        else:
            summ = (md_section(readme, r"^1\. Population distribution", pre, hshift=0) +
                    "<h2 id=\"origin-age\">Single origin and relative age: results</h2>" +
                    md_section(readme, r"^Results$", pre, hshift=0).split("</h3>", 1)[-1] +
                    md_section(readme, r"^Interpretation$", pre, hshift=0) +
                    md_section(readme, r"^3\. Great apes", pre, hshift=0))
        cav = md_section(readme, r"^Caveats", pre, hshift=0)
        body = f"""
<p><a href="index.html">Findings</a> &rsaquo; {esc(v['short'])}</p>
<h1>{esc(v['title'])}</h1>
{note}
<p>Full methods, statistics and reproduction instructions: <a href="{BLOB}{readme}">{esc(readme)}</a>
(results tables in <a href="{TREE}{d}/results">{esc(d)}/results/</a>; reproduce with <code>bash {esc(d)}/run.sh</code>).</p>
{summ}
<h2 id="figures">Figures</h2>
{analysis_figs(d, pre)}
{cav}
"""
        write(f"findings/{key}.html", page(f"findings/{key}.html", v["short"], body, active="findings/index.html"))


def page_methods(D, pal):
    prov = D["prov"]
    params = prov.get("parameters", {})
    tools = prov.get("tools", {})
    src = prov.get("sources", {})
    ptab = "".join(f"<tr><th scope=row><code>{esc(k)}</code></th><td>{esc(str(v))}</td></tr>" for k, v in params.items())
    ttab = "".join(f"<tr><th scope=row>{esc(k)}</th><td>{esc(str(v))}</td></tr>" for k, v in tools.items())
    stab = []
    for k, v in src.items():
        if isinstance(v, dict):
            parts = []
            for kk, vv in v.items():
                vv = str(vv)
                if vv.startswith("http"):
                    vv = f'<a href="{esc(vv)}">{esc(vv)}</a>'
                elif kk == "doi":
                    vv = f'<a href="https://doi.org/{esc(vv)}">doi:{esc(vv)}</a>'
                else:
                    vv = esc(vv)
                parts.append(f"{esc(kk)}: {vv}")
            stab.append(f"<tr><th scope=row>{esc(k)}</th><td>{'<br>'.join(parts)}</td></tr>")
        else:
            stab.append(f"<tr><th scope=row>{esc(k)}</th><td>{esc(str(v))}</td></tr>")
    pfdf = D["pf"][[c for c in ["family", "n_members", "members", "min_cotie_rate", "median_cotie_rate", "tied_points"]
                    if c in D["pf"].columns]]
    pf_tbl = table_html(pfdf, [(c, c) for c in pfdf.columns], num=("n_members", "min_cotie_rate", "median_cotie_rate",
                                                                    "tied_points"), wrap_class="table-wrap free")
    swatches = "".join(f'<span><span class="sw" style="background:{h}"></span>{esc(n)}'
                       f'{(" (" + esc(pal["family_of"][a]) + ")") if a in pal["family_of"] else ""} <code>{h}</code></span>'
                       for n, h, a in zip(pal["names"], pal["hex"], pal["arm"]))
    q = D["palqc"]
    qtxt = ""
    if q is not None:
        w = q[(q.cooccurrence_weight >= 50) & ~q.same_family]
        nf = q[~q.same_family]
        sf = q[q.same_family]
        qtxt = (f"For the {len(w)} pairs of different families or ends that co-occur in the mosaics (weight &ge; 50), the "
                f"smallest CIEDE2000 difference is {w.dE2000.min():.1f} in normal vision and "
                f"{w.dE2000_deuteranopia.min():.1f} in simulated deuteranopia; over all {len(nf)} such pairs it is "
                f"{nf.dE2000.min():.1f}; shades within a family differ by &ge; {sf.dE2000.min():.1f}.")
    caps = D["captions"]
    qc10 = md_section("data/catalogue/QC.md", r"^10\. Paralog families", "", hshift=2, id_prefix="qc10-")
    body = f"""
<h1>Methods</h1>
<h2 id="overview">Overview</h2>
<ol>
<li><strong>Windows.</strong> For every telomeric contig end of the HPRC Release 2 assemblies, the outermost
{params.get('window_bp', 100000) // 1000} kb were range-read from the public assemblies and oriented telomere &rarr; centromere
(tip frame); CHM13 v2.0, HG002 Q100 and WashU pedigree ends were added (pedigree ends are used only as validation
evidence). The terminal (CCCTAA)n block is masked and all coordinates are measured from the telomere/subtelomere boundary.</li>
<li><strong>Home arm.</strong> From the contig name; otherwise from 10 kb probes at 250 kb and 1 Mb mapped to CHM13
(calibrated on named tips); otherwise from the mosaic. Acrocentric short arms and the PARs are pooled
(<code>acroP</code>, <code>PAR1</code>, <code>PAR2</code>).</li>
<li><strong>Mosaics (revision 2).</strong> One alignment of each window to the outer 500 kb of the 48 CHM13 ends gives
per-base identity tracks; a sliding identity on a 100 bp grid anchored at the contig tip is labelled home, paralog family,
CHM13 arm or none, short runs are merged and boundaries refined at base level. The grid does not depend on where the
telomere boundary is called, which makes alleles stable under boundary shifts (<a href="validation.html#stability">validation</a>).</li>
<li><strong>Structural alleles.</strong> Core structures (blocks &ge; 10 kb) with switch positions agreeing within 10 kb
(complete linkage) define <code>SA-&lt;arm&gt;-NN</code>; partial and low-confidence ends are assigned to the compatible
allele with the nearest median switches. Identifiers changed in revision 2; see <a href="{BLOB}data/catalogue/id_history.tsv">id_history.tsv</a>.</li>
<li><strong>Haplogroups.</strong> Within-arm all-vs-all FracMinHash ANI over the first
{params.get('haplogroup_compared_bp', 80000) // 1000} kb of subtelomere, UPGMA cut at ANI {params.get('haplogroup_ani', 0.995)},
validated by alignment to the medoid (<code>HG-&lt;arm&gt;-NN</code>).</li>
<li><strong>Evidence and tiers.</strong> Recurrence in &ge; 3 unrelated samples plus orthogonal lines (Verkko assemblies of
the same haplotypes, HG002 Q100 reads or hifiasm read layouts across switches, intact pedigree transmission, presence in a
T2T assembly) give tiers A/B/C; see <a href="validation.html">validation</a>.</li>
</ol>
<h2 id="labels">Paralog families, the home rule and nativity classes</h2>
{definitions_html(D)}
<p>Paralog families (from <a href="{BLOB}data/catalogue/paralog_families.tsv">paralog_families.tsv</a>):</p>
{pf_tbl}
<h3 id="labels-details">Details and worked example (QC report section 10)</h3>
{qc10}
<h2 id="parameters">Parameters and software</h2>
<div class="grid"><div class="table-wrap free"><table><caption class="muted">Build parameters (PROVENANCE.json)</caption>
<tbody>{ptab}</tbody></table></div>
<div class="table-wrap free"><table><caption class="muted">Software versions</caption><tbody>{ttab}</tbody></table></div></div>
<h3 id="sources">Data sources</h3>
<div class="table-wrap free"><table><tbody>{''.join(stab)}</tbody></table></div>
<h2 id="pipeline">Pipeline steps</h2>
{md_file('pipeline/CATALOGUE.md', '', hshift=2, id_prefix='pipe-').split('</h3>', 1)[-1]}
<h2 id="palette">Colour palette</h2>
<p>Non-home blocks are coloured by their best single CHM13 v2.0 end and home blocks by the arm's own end (a paralog
tied within the home margin is not drawn as foreign). Members of a paralog family share one hue and differ in
lightness, so a family-labelled block reads as its family whichever member is its best end; every other end has its own
colour. Colours were chosen once by a seeded optimisation
(<a href="{BLOB}pipeline/atlas_palette.py">pipeline/atlas_palette.py</a>) that maximises the CIEDE2000 colour difference,
in normal vision and in simulated deuteranopia (Machado et al. 2009, doi:10.1109/TVCG.2009.113), between CHM13 ends that are adjacent or co-occur in
the catalogue mosaics. {qtxt} Pairwise values: <a href="{BLOB}figures/atlas/palette_qc.tsv">palette_qc.tsv</a>; hex
codes: <a href="{BLOB}figures/atlas/PALETTE.tsv">PALETTE.tsv</a>. With 48 classes some pairs that never co-occur remain
similar; the arm pages name the family and best end of every block on hover.</p>
{figure('figures/atlas/palette_chm13_ends.png', 'Palette key: 24 chromosomes by p and q arm, with the paralog families.', md_inline(caps.get('palette_chm13_ends', '')), '')}
<div class="legend">{swatches}</div>
<h2 id="schema">Data dictionary</h2>
{md_file('data/catalogue/SCHEMA.md', '', hshift=2, id_prefix='schema-').split('</h3>', 1)[-1]}
"""
    write("methods.html", page("methods.html", "Methods", body, active="methods.html"))


TAN_KEYS = [("n_pairs", "pairs compared"), ("n_samples", "samples"), ("identity_median", "median sequence identity"),
            ("frac_boundary_exact", "telomere boundary identical"), ("frac_boundary_le_100bp", "boundary within 100 bp"),
            ("S20_state_agree", "S20 state agrees"), ("form_agree", "form agrees"),
            ("core_structure_agree", "core structure agrees"), ("sa_agree", "structural allele agrees"),
            ("sa_compatible", "structural allele compatible"), ("frac_tan_unassigned", "Tan record unassigned"),
            ("switch_base_frac_proj_le_10bp", "base-level switches within 10 bp")]


def tan_metrics(tan):
    sh = tan.get("same_haplotype", {}) if isinstance(tan, dict) else {}
    blocks = [(k, sh[k]) for k in ["primary", "all_pairs"] if isinstance(sh.get(k), dict)]
    if not blocks:
        return ""
    head = "".join(f"<th scope=col class=num>{esc(k.replace('_', ' '))}</th>" for k, _ in blocks)
    rows = []
    for key, lab in TAN_KEYS:
        if not any(key in b for _, b in blocks):
            continue
        cells = []
        for _, b in blocks:
            v = b.get(key, sh.get(key) if key == "n_samples" else None)
            cells.append(f"<td class=num>{'' if v is None else (fnum(v) if isinstance(v, int) else f'{v:.4g}')}</td>")
        rows.append(f"<tr><th scope=row>{esc(lab)} <code>{esc(key)}</code></th>{''.join(cells)}</tr>")
    return (f'<div class="table-wrap free"><table><caption class="muted">Same-haplotype comparison, key metrics from '
            f'tan2026_summary.json</caption><thead><tr><th scope=col>metric</th>{head}</tr></thead><tbody>'
            f'{"".join(rows)}</tbody></table></div>')


def page_validation(D):
    e = D["ends"]
    caps = D["captions"]
    tiers_par = ""
    with open(os.path.join(CAT, "SCHEMA.md"), encoding="utf-8") as fh:
        txt = fh.read()
    m = re.search(r"(Tiers\. Evidence lines:.*?)(?:\n\n|\n##)", txt, re.S)
    if m:
        tiers_par = MD("data/catalogue/SCHEMA.md", "").render(m.group(1))
    ev_rows = []
    for col, lab in [("vk_status", "Verkko concordance"), ("lay_status", "hifiasm read layout"),
                     ("reads_status", "HG002 Q100 reads across switch"), ("ped_status", "Pedigree transmission")]:
        vc = e[col].fillna("(empty)").value_counts()
        ev_rows.append(f"<tr><th scope=row>{lab}<br><code>{col}</code></th><td>" +
                       ", ".join(f"{esc(str(k))}: {fnum(v)}" for k, v in vc.items()) + "</td></tr>")
    stab = md_section("data/catalogue/QC.md", r"^9\. Track-based", "", hshift=1, id_prefix="stab-")
    stab = stab.split("</h2>", 1)[-1] if stab.startswith("<h2") else stab
    if not stab:
        stab = '<p class="warn">The stability section of QC.md was not found at build time.</p>'
    # query-tool validation
    vpath = "tests/validation/VALIDATION.md"
    if os.path.exists(rp(vpath)):
        qv = md_file(vpath, "", hshift=2, id_prefix="qv-")
        qv = qv.split("</h3>", 1)[-1] if qv.startswith("<h3") else qv
        qv = (f'<p>Rendered from <a href="{BLOB}{vpath}">{vpath}</a> (per-test tables in '
              f'<a href="{TREE}tests/validation">tests/validation/</a>).</p>' + qv)
    else:
        qv = '<p class="warn">Query-tool validation in progress: tests/validation/VALIDATION.md is not present.</p>'
    # Tan et al. cross-check
    tan = load_json(rp("data", "validation", "tan2026_summary.json"))
    if tan is not None or os.path.exists(rp("data", "validation", "README.md")):
        parts = [f'<p>Same-haplotype comparison of catalogue ends with the independent Tan et al. 2026 records '
                 f'(<a href="https://doi.org/{TAN}">doi:{TAN}</a>; sequences CC BY 4.0, '
                 f'<a href="https://doi.org/{TAN_ZENODO}">doi:{TAN_ZENODO}</a>), from '
                 f'<a href="{TREE}data/validation">data/validation/</a>.</p>']
        if tan is not None:
            parts.append(tan_metrics(tan))
        if os.path.exists(rp("data", "validation", "README.md")):
            parts.append(md_section("data/validation/README.md", r"^Results", "", hshift=1, id_prefix="tan-"))
            parts.append(md_section("data/validation/README.md", r"^Limits", "", hshift=1, id_prefix="tanl-"))
        if tan is not None:
            parts.append('<details><summary>All values in tan2026_summary.json</summary>' + render_json_generic(tan) +
                         '</details>')
        tan_html = "".join(parts)
    else:
        tan_html = ("<p class=\"warn\">Cross-check in progress. When <code>data/validation/tan2026_summary.json</code> is "
                    "present, this section reports how the atlas's structural alleles correspond to the independent "
                    f"Tan et al. 2026 sequences (<a href=\"https://doi.org/{TAN}\">doi:{TAN}</a>).</p>")
    qc = md_file("data/catalogue/QC.md", "", hshift=1, id_prefix="qc-")
    qc = qc.split("</h2>", 1)[-1] if qc.startswith("<h2") else qc
    body = f"""
<h1>Validation</h1>
<p>Contents: <a href="#tiers">evidence tiers</a> &middot; <a href="#stability">boundary-shift stability</a> &middot;
<a href="#query-tool">query-tool validation</a> &middot; <a href="#tan2026">cross-check with Tan et al. 2026</a> &middot;
<a href="#evidence">per-end evidence</a> &middot; <a href="#qc">full QC report</a></p>
<h2 id="tiers">Evidence tiers</h2>
<p>Every structural allele and haplogroup carries an evidence tier (from <a href="{BLOB}data/catalogue/SCHEMA.md">SCHEMA.md</a>):</p>
{tiers_par}
{figure('figures/atlas/F4_validation.png', caption_title(caps.get('F4_validation')) + ': tier counts, UpSet of evidence lines, and per-allele evidence for common alleles.', md_inline(caps.get('F4_validation', '')), '')}
<h2 id="stability">Catalogue stability under telomere-boundary shifts</h2>
<p>Rendered from section 9 of <a href="{BLOB}data/catalogue/QC.md">QC.md</a>: allele, structure and S20 calls after
re-aligning windows from boundaries shifted by &plusmn;100 to &plusmn;1000 bp, compared with the unshifted calls.</p>
{stab}
<h2 id="query-tool">Query-tool validation</h2>
{qv}
<h2 id="tan2026">Cross-check with Tan et al. 2026</h2>
{tan_html}
<h2 id="evidence">Per-end evidence status</h2>
<p>Counts over all {fnum(len(e))} ends in <code>ends.tsv</code> (statuses are defined in the data dictionary).</p>
<div class="table-wrap free"><table><tbody>{''.join(ev_rows)}</tbody></table></div>
<h2 id="qc">Catalogue QC report</h2>
<p>Rendered from <a href="{BLOB}data/catalogue/QC.md">data/catalogue/QC.md</a> (generated by the catalogue pipeline).</p>
{qc}
"""
    write("validation.html", page("validation.html", "Validation", body, active="validation.html"))


# ------------------------------------------------------------------ main
def main():
    os.makedirs(DOCS, exist_ok=True)
    D = load()
    pal, pidx = palette(D)
    H = headline(D)
    AG = ages(D)
    AS = arm_stats(D)
    # static assets
    for src, dst in [("site.css", "assets/css/site.css"), ("common.js", "assets/js/common.js"),
                     ("arm.js", "assets/js/arm.js"), ("alleles.js", "assets/js/alleles.js")]:
        os.makedirs(os.path.dirname(os.path.join(DOCS, dst)), exist_ok=True)
        shutil.copyfile(os.path.join(STATIC, src), os.path.join(DOCS, dst))
    write(".nojekyll", "")
    arms_data, index, repb = build_data(D, pal, pidx)
    write("data/palette.js", js_payload("ATLAS_PALETTE", pal))
    write("data/alleles_index.js", js_payload("ATLAS_INDEX", index))
    write("data/alleles_index.json", json.dumps(index, separators=(",", ":"), allow_nan=False))
    for arm, data in arms_data.items():
        write(f"data/arms/{arm}.js", js_payload("ATLAS_ARM", data))
    page_index(D, H, AG, repb)
    page_arms_index(D, AS)
    for arm in AS.arm:
        page_arm(D, arm, arms_data[arm], AS, pal, repb)
    page_alleles(D, index)
    page_findings(D)
    page_methods(D, pal)
    page_validation(D)
    page_query()
    page_downloads(D)
    page_literature()
    write("README.md", "# docs/ (generated)\n\nStatic GitHub Pages site generated by `website/build_site.py`; do not edit "
                       "by hand. Serve from the `docs/` folder of the main branch.\n")
    total = 0
    for dp, _, fs in os.walk(DOCS):
        for f in fs:
            total += os.path.getsize(os.path.join(dp, f))
    print(f"docs/ written: {total / 1e6:.2f} MB; arms {len(arms_data)}; headline {json.dumps({k: H[k] for k in ['n_ends', 'n_sa', 'n_hg']})}; ages {AG}")


if __name__ == "__main__":
    main()
