#!/usr/bin/env python3
"""Atlas overview figures (F1-F4 + palette key) from the v0.1 catalogue release tables.

Usage (from the repository root):
    python pipeline/figures_atlas.py            # writes figures/atlas/*.{png,pdf,svg}, PALETTE.tsv, CAPTIONS.md

Inputs (read only): data/catalogue/{alleles_structural,haplogroups,allele_frequencies_by_superpop,ends,mosaic_blocks}.tsv.
Every number drawn or written in a caption is computed here from those files.

Figures
  F1 genome-wide overview: arms on a ring, structural-allele counts per arm, terminal foreign forms (tier A/B) as chords
     from the donor-like arm to the recipient arm, chord width = summed Release 2 frequency; per-arm allele composition.
  F2 allele mosaics of real ends for 20q, 3q, 6p, 8p, 16p, 11p, grouped by structural allele.
  F3 common forms by superpopulation (tier A non-native alleles with Release 2 frequency >= 0.05).
  F4 validation: tiers, evidence-line combinations and per-allele evidence counts for the common forms.
  palette: the 48-colour CHM13 end palette (also written to PALETTE.tsv; reused by the website).
"""
import json
import math
import os
import sys
from collections import Counter, OrderedDict, defaultdict

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CAT = os.path.join(ROOT, "data", "catalogue")
OUT = os.path.join(ROOT, "figures", "atlas")

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POOLED = {"13p": "acroP", "14p": "acroP", "15p": "acroP", "21p": "acroP", "22p": "acroP",
          "Xp": "PAR1", "Yp": "PAR1", "Xq": "PAR2", "Yq": "PAR2"}
# CHM13 end drawn for the 'home' sequence of a pooled group
POOLED_HOME = {"acroP": "15p", "PAR1": "Xp", "PAR2": "Xq"}
NONE_COLOUR = "#ececec"
TELO_COLOUR = "#3b3b3b"
TIER_AB = "#4d4d4d"
TIER_C = "#b3b3b3"
CLASS_COLOURS = {  # allele classes (greyscale, so they never collide with the chromosome hues)
    "native": "#e3e3e3",
    "tier C (any structure)": "#b3b3b3",
    "internal foreign block (tier A/B)": "#6b6b6b",
    "terminal foreign block (tier A/B)": "#1c1c1c",
    "unassigned": "#ffffff",
}


# ---------------------------------------------------------------- style
def style():
    """Publication rcParams: 8/7/6 pt role ladder, open frame, outward ticks, frameless legends, 300 dpi, Type-42 fonts."""
    rc = {
        "font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": False, "legend.frameon": False, "figure.dpi": 200,
        "savefig.dpi": 300, "savefig.bbox": "tight", "axes.titleweight": "normal",
        "axes.titlelocation": "left", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    }
    matplotlib.rcParams.update(rc)


def panel_letter(ax, letter, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")


# ---------------------------------------------------------------- palette (family-aware; see pipeline/atlas_palette.py)
def _ensure_palette():
    p = os.path.join(OUT, "PALETTE.tsv")
    ok = os.path.exists(p) and "family" in pd.read_csv(p, sep="\t", nrows=1).columns
    if not ok:
        sys.path.insert(0, os.path.dirname(__file__))
        import atlas_palette
        atlas_palette.main(force=True)
    return pd.read_csv(p, sep="\t")


def _families():
    pf = pd.read_csv(os.path.join(CAT, "paralog_families.tsv"), sep="\t")
    return {r.family: r.members.split(",") for r in pf.itertuples()}


PALDF = _ensure_palette()


def _palette_qc_text():
    p = os.path.join(OUT, "palette_qc.tsv")
    if not os.path.exists(p):
        return ""
    q = pd.read_csv(p, sep="\t")
    w = q[(q.cooccurrence_weight >= 50) & ~q.same_family]
    nf = q[~q.same_family]
    sf = q[q.same_family]
    return (f"For the {len(w)} pairs of different families or ends with co-occurrence weight >= 50, minimum dE2000 is "
            f"{w.dE2000.min():.1f} (simulated deuteranopia {w.dE2000_deuteranopia.min():.1f}); over all "
            f"{len(nf)} such pairs {nf.dE2000.min():.1f}; between shades of one family >= {sf.dE2000.min():.1f}.")


PALETTE_QC_TEXT = _palette_qc_text()
PAL = {r.arm: r.hex for r in PALDF.itertuples() if r.arm != "none"}
FAMILIES = _families()
FAMILY_OF = {m: f for f, ms in FAMILIES.items() for m in ms}
with open(os.path.join(CAT, "catalogue.json")) as _fh:
    VERSION = json.load(_fh).get("version", "v0.1")


GROUP_ARMS = {"acroP": ["13p", "14p", "15p", "21p", "22p"], "PAR1": ["Xp", "Yp"], "PAR2": ["Xq", "Yq"]}


def home_end(cat_arm, best_end):
    """Colour key for a home-labelled block: the arm's own CHM13 end (pooled groups: the best end if it is a member of
    the group, else the group's default end). Under the home rule a home block's best single end can be a tied paralog,
    which must not be drawn as foreign colour."""
    grp = GROUP_ARMS.get(cat_arm)
    if grp:
        be = best_end[3:] if isinstance(best_end, str) and best_end.startswith("chr") else best_end
        return "chr" + (be if be in grp else POOLED_HOME[cat_arm])
    return "chr" + cat_arm


def arm_colour(token, cat_arm=None):
    """Colour for a CHM13 end ('chr3q'), an arm token ('3q'), 'home', 'none' or a family label ('PF04', drawn with the
    shade of its first member)."""
    if token is None or (isinstance(token, float) and np.isnan(token)) or token == "none":
        return NONE_COLOUR
    t = token[3:] if token.startswith("chr") else token
    if t == "home":
        t = POOLED_HOME.get(cat_arm, cat_arm)
    if t in FAMILIES:
        t = FAMILIES[t][0]
    return PAL.get(t, NONE_COLOUR)


# ---------------------------------------------------------------- data
def load():
    sa = pd.read_csv(os.path.join(CAT, "alleles_structural.tsv"), sep="\t")
    hg = pd.read_csv(os.path.join(CAT, "haplogroups.tsv"), sep="\t")
    fr = pd.read_csv(os.path.join(CAT, "allele_frequencies_by_superpop.tsv"), sep="\t")
    ends = pd.read_csv(os.path.join(CAT, "ends.tsv"), sep="\t", low_memory=False,
                       usecols=["end_id", "set", "sample", "superpopulation", "catalogue_arm", "structural_allele_id",
                                "sa_assignment", "haplogroup_id", "sub_len", "nativity_class"])
    mb = pd.read_csv(os.path.join(CAT, "mosaic_blocks.tsv"), sep="\t", low_memory=False,
                     usecols=["end_id", "level", "block", "label", "family", "source", "best_end", "is_home",
                              "sub_start", "sub_end"])
    return sa, hg, fr, ends, mb


def arm_order(arms):
    """Genome order; pooled groups placed at the position of their first member."""
    key = {}
    for i, c in enumerate(CHROMS):
        key[c + "p"] = 2 * i
        key[c + "q"] = 2 * i + 1
    key["acroP"] = key["13p"]
    key["PAR1"] = key["Xp"]
    key["PAR2"] = key["Xq"]
    return sorted(arms, key=lambda a: key[a])


def arm_chrom(a):
    return {"acroP": "13-15,21,22", "PAR1": "X/Y", "PAR2": "X/Y"}.get(a, a[:-1])


def allele_class(row):
    if row is None:
        return "unassigned"
    if row["is_native"]:
        return "native"
    if row["tier"] == "C":
        return "tier C (any structure)"
    return "terminal foreign block (tier A/B)" if isinstance(row["form"], str) else "internal foreign block (tier A/B)"


def allele_label(r):
    st = r["structure"]
    return f'{r["allele_id"]}  {st}'


# ---------------------------------------------------------------- render checks
def bbox_overlaps(fig, ignore=()):
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible() and t not in ignore]
    over = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                over.append((a.get_text()[:30], b.get_text()[:30]))
    return over


def save(fig, name, svg=True):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    if svg:
        fig.savefig(os.path.join(OUT, name + ".svg"))
    over = bbox_overlaps(fig)
    plt.close(fig)
    return over


# ---------------------------------------------------------------- palette key (shared helper)
def draw_palette_key(ax, used=None, title="Colour = CHM13 v2.0 end the block resembles most (best end)"):
    """Top: 24 x 2 grid (columns chromosomes, rows p / q). Below: the paralog families, whose members share a hue and
    differ in lightness (outlined in the grid). Ends not drawn in the figure are faded."""
    for i, c in enumerate(CHROMS):
        for j, s in enumerate("pq"):
            e = c + s
            alpha = 1.0 if (used is None or e in used) else 0.18
            ax.add_patch(Rectangle((i, 1 - j), 0.86, 0.86, facecolor=PAL[e], edgecolor="none", alpha=alpha))
            if e in FAMILY_OF:
                ax.add_patch(Rectangle((i, 1 - j), 0.86, 0.86, facecolor="none", edgecolor="#1c1c1c", lw=0.6,
                                       alpha=max(alpha, 0.5)))
        ax.text(i + 0.43, -0.12, c, ha="center", va="top", fontsize=6)
    ax.text(-0.25, 1.43, "p", ha="right", va="center", fontsize=6)
    ax.text(-0.25, 0.43, "q", ha="right", va="center", fontsize=6)
    ax.add_patch(Rectangle((24.8, 0.5), 0.86, 0.86, facecolor=NONE_COLOUR, edgecolor="#9a9a9a", lw=0.4))
    ax.text(25.9, 0.93, "none", ha="left", va="center", fontsize=6)
    ax.text(-0.6, 2.15, title, ha="left", va="bottom", fontsize=7)
    ax.text(-0.6, -1.85, "Paralog families\n(outlined; shared hue):", ha="left", va="center", fontsize=6)
    for k, (f, ms) in enumerate(sorted(FAMILIES.items())):
        row, col = divmod(k, 2)
        x = 8.0 + col * 14.5
        y = -1.0 - row * 0.85
        ax.text(x, y, f, ha="left", va="center", fontsize=6)
        for m_i, m in enumerate(ms):
            ax.add_patch(Rectangle((x + 1.9 + m_i * 0.6, y - 0.32), 0.54, 0.64, facecolor=PAL[m], edgecolor="none"))
        ax.text(x + 2.05 + len(ms) * 0.6, y, ",".join(ms), ha="left", va="center", fontsize=6)
    ax.set_xlim(-1.0, 38.0)
    ax.set_ylim(-3.3, 2.9)
    ax.set_axis_off()


# ---------------------------------------------------------------- F1
NATIVITY_COLOURS = OrderedDict([  # greyscale so they never collide with the CHM13 end colours
    ("native", "#e3e3e3"), ("within_family", "#8c8c8c"), ("foreign", "#1c1c1c"), ("unassigned", "#ffffff")])


def terminal_blocks(ends, mb):
    """Per end: label of the outermost core block and the best single CHM13 end of the first raw block carrying that
    label (raw blocks keep best_end; core blocks keep the label)."""
    core1 = mb[(mb.level == "core") & (mb.block == 1)].set_index("end_id").label.to_dict()
    raw = mb[mb.level == "raw"].sort_values(["end_id", "block"])
    out = {}
    for e, lab in core1.items():
        out[e] = [lab, None]
    for e, lab, be in zip(raw.end_id, raw.label, raw.best_end):
        t = out.get(e)
        if t is not None and t[1] is None and lab == t[0] and isinstance(be, str):
            t[1] = be
    return out


def fig1(sa, ends, mb, captions):
    arms = arm_order(sa.arm.unique().tolist())
    r2 = ends[ends.set == "HPRC_R2"]
    sai = sa.set_index("allele_id")
    # (b) nativity class of Release 2 ends (class of the end's own core structure)
    comp = {}
    for a, d in r2.groupby("catalogue_arm"):
        if a not in arms:
            continue
        cls = d.nativity_class.fillna("unassigned")
        n = len(d)
        comp[a] = {k: float((cls == k).sum()) / n for k in NATIVITY_COLOURS}
        comp[a]["n"] = n
    tiers = sa.groupby(["arm", "tier"]).size().unstack(fill_value=0).reindex(arms).fillna(0)
    # chords: Release 2 ends in tier A/B alleles whose outermost core block is not home/none
    tb = terminal_blocks(r2, mb)
    tierAB = set(sa.loc[sa.tier.isin(["A", "B"]), "allele_id"])
    rows = []
    for e, arm, aid, cls in zip(r2.end_id, r2.catalogue_arm, r2.structural_allele_id, r2.nativity_class):
        if aid not in tierAB or e not in tb:
            continue
        lab, be = tb[e]
        if lab in ("home", "none") or not isinstance(be, str):
            continue
        rows.append((arm, be[3:], lab, cls))
    ch = pd.DataFrame(rows, columns=["arm", "donor", "label", "nativity"])
    denom = r2[r2.structural_allele_id.notna()].groupby("catalogue_arm").size()
    chords = ch.groupby(["arm", "donor"]).agg(n=("label", "size"), label=("label", lambda x: x.mode().iloc[0]),
                                              within=("nativity", lambda x: (x == "within_family").mean())).reset_index()
    chords["freq"] = chords.n / chords.arm.map(denom)
    chords["donor_pos"] = chords.donor.map(lambda t: POOLED.get(t, t))
    chords = chords[(chords.arm != chords.donor_pos) & chords.donor_pos.isin(arms)].sort_values("freq")

    fig = plt.figure(figsize=(7.2, 5.1))
    axc = fig.add_axes([0.0, 0.15, 0.62, 0.81])
    axb = fig.add_axes([0.71, 0.25, 0.27, 0.67])
    axk = fig.add_axes([0.0, -0.035, 0.68, 0.15])

    gaps, chrom_prev = [], None
    for a in arms:
        chm = arm_chrom(a)
        gaps.append(0.0 if chm == chrom_prev else 0.35)
        chrom_prev = chm
    unit = 2 * math.pi / (len(arms) + sum(gaps))
    pos, cur = {}, 0.0
    for a, g in zip(arms, gaps):
        cur += g * unit
        pos[a] = (cur, cur + unit * 0.92)
        cur += unit

    def ang(t):
        return math.pi / 2 - t

    def pt(t, r):
        return r * math.cos(ang(t)), r * math.sin(ang(t))

    R0, R1, RL, RB0, RBMAX = 1.0, 1.07, 1.10, 1.30, 1.62
    maxn = tiers.sum(axis=1).max()
    for a in arms:
        t0, t1 = pos[a]
        ts = np.linspace(t0, t1, 12)
        home = POOLED_HOME.get(a, a)
        xy = [pt(t, R0) for t in ts] + [pt(t, R1) for t in ts[::-1]]
        axc.add_patch(plt.Polygon(xy, closed=True, facecolor=PAL[home], edgecolor="none"))
        r = RB0
        for tr, colr in [("AB", TIER_AB), ("C", TIER_C)]:
            v = sum(tiers.loc[a, x] for x in tr if x in tiers.columns)
            h = (RBMAX - RB0) * v / maxn
            if h > 0:
                xy = [pt(t, r) for t in ts] + [pt(t, r + h) for t in ts[::-1]]
                axc.add_patch(plt.Polygon(xy, closed=True, facecolor=colr, edgecolor="none"))
            r += h
        tm = (t0 + t1) / 2
        x, y = pt(tm, RL)
        deg = math.degrees(ang(tm)) % 360
        right = math.cos(ang(tm)) >= 0
        axc.text(x, y, a, rotation=deg if right else deg - 180, rotation_mode="anchor",
                 ha="left" if right else "right", va="center", fontsize=6)
    wmax = 7.0
    for _, c in chords.iterrows():
        tr = sum(pos[c.arm]) / 2
        td = sum(pos[c.donor_pos]) / 2
        p0, p1 = pt(td, R0 - 0.01), pt(tr, R0 - 0.01)
        d = abs(((tr - td + math.pi) % (2 * math.pi)) - math.pi)
        k = 0.55 * (1 - d / math.pi) ** 1.3
        ctrl = pt((td + tr) / 2 if abs(tr - td) <= math.pi else (td + tr) / 2 + math.pi, k)
        path = Path([p0, ctrl, p1], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        within = c.within > 0.5
        axc.add_patch(PathPatch(path, facecolor="none", edgecolor=arm_colour(c.donor), lw=0.35 + wmax * c.freq,
                                alpha=0.9 if c.freq >= 0.05 else 0.6, capstyle="butt", zorder=2 + c.freq,
                                linestyle=(0, (2.0, 1.2)) if within else "-"))
        axc.plot(*pt(tr, R0 - 0.035), marker="o", ms=1.2 + 6 * c.freq, color="#1c1c1c", mec="white", mew=0.3,
                 zorder=5)
    top = chords.sort_values("freq", ascending=False).head(4)
    axc.set_xlim(-1.80, 1.80)
    axc.set_ylim(-1.80, 1.80)
    axc.set_aspect("equal")
    axc.set_axis_off()
    axc.text(-1.62, 1.80, "Terminal blocks resembling another arm (tier A/B alleles)", fontsize=8, ha="left",
             va="top")
    panel_letter(axc, "a", x=0.04, y=0.975)
    lx, ly = 1.18, -1.32
    axc.text(lx, ly + 0.34, "Frequency on\nrecipient arm", fontsize=6, ha="left", va="bottom")
    for i, f in enumerate([0.05, 0.2, 0.45]):
        yy = ly + 0.24 - i * 0.13
        axc.plot([lx, lx + 0.18], [yy, yy], color="#555555", lw=0.35 + wmax * f, solid_capstyle="butt")
        axc.text(lx + 0.24, yy, f"{f:g}", fontsize=6, va="center")
    axc.plot([lx, lx + 0.18], [ly - 0.17, ly - 0.17], color="#555555", lw=1.2, linestyle=(0, (2.0, 1.2)))
    axc.text(lx + 0.24, ly - 0.17, "within family", fontsize=6, va="center")
    axc.plot(lx + 0.05, ly - 0.32, marker="o", ms=3, color="#1c1c1c", mec="white", mew=0.3)
    axc.text(lx + 0.13, ly - 0.32, "recipient arm", fontsize=6, va="center")
    lx2, ly2 = -1.78, -1.32
    axc.text(lx2, ly2 + 0.34, "Structural alleles\nper arm (outer bars)", fontsize=6, ha="left", va="bottom")
    for i, (lab, colr) in enumerate([("tier A/B", TIER_AB), ("tier C", TIER_C)]):
        yy = ly2 + 0.24 - i * 0.13
        axc.add_patch(Rectangle((lx2, yy - 0.045), 0.12, 0.09, facecolor=colr, edgecolor="none"))
        axc.text(lx2 + 0.17, yy, lab, fontsize=6, va="center")
    amax = tiers.sum(axis=1).idxmax()
    axc.text(lx2, ly2 - 0.17, f"longest bar: {amax}, {int(maxn)}", fontsize=6, va="center")

    order = sorted(comp, key=lambda a: comp[a]["native"])
    y = np.arange(len(order))
    left = np.zeros(len(order))
    for k, col in NATIVITY_COLOURS.items():
        v = np.array([comp[a][k] for a in order])
        axb.barh(y, v, left=left, color=col, edgecolor="#9a9a9a" if k == "unassigned" else "none", lw=0.3,
                 height=0.82, label=k)
        left += v
    axb.set_yticks(y)
    axb.set_yticklabels(order)
    axb.set_ylim(-0.7, len(order) - 0.3)
    axb.set_xlim(0, 1.0)
    axb.set_xticks([0, 0.5, 1.0])
    axb.set_xticklabels(["0", "0.5", "1"])
    axb.set_xlabel("Fraction of Release 2 ends", labelpad=2)
    axb.tick_params(axis="y", length=0, pad=1)
    mostly = [a for a in order if comp[a]["native"] < 0.5]
    axb.set_title(f"Native ends are a minority\non {len(mostly)} of {len(order)} arms", fontsize=8, x=-0.16)
    panel_letter(axb, "b", x=-0.20, y=1.02)
    lab_txt = {"native": "native (home, none)", "within_family": "within family", "foreign": "foreign",
               "unassigned": "no core structure"}
    handles = [Patch(facecolor=c, edgecolor="#9a9a9a" if k in ("unassigned", "native") else "none", lw=0.3,
                     label=lab_txt[k]) for k, c in NATIVITY_COLOURS.items()]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.705, 0.155), ncol=2, fontsize=6,
               handlelength=1.0, handletextpad=0.4, borderaxespad=0, labelspacing=0.3, columnspacing=0.8)
    used = set([POOLED_HOME.get(a, a) for a in arms]) | set(chords.donor)
    draw_palette_key(axk, used=used, title="Colour = CHM13 v2.0 end (ring: own arm; chord: best end of the "
                                           "terminal block; faded: not drawn)")
    over = save(fig, "F1_genome_overview")
    t4 = "; ".join(f"{r.arm}>{r.donor} {r.freq:.2f}" for r in top.itertuples())
    nwf = int((chords.within > 0.5).sum())
    nat_tot = r2.nativity_class.fillna("unassigned").value_counts()
    captions["F1_genome_overview"] = (
        f"**Figure 1 | Genome-wide organisation of structural alleles at human chromosome ends (catalogue "
        f"{VERSION}).** (a) The {len(arms)} catalogue arms (pooled groups: acroP = 13p/14p/15p/21p/22p, PAR1 = Xp/Yp, "
        f"PAR2 = Xq/Yq) in genome order, clockwise from 1p; ring colour = the arm's own CHM13 v2.0 end (pooled groups: "
        f"15p, Xp, Xq). Outer bars: structural alleles per arm, tier A/B and tier C (longest bar {amax}, "
        f"{int(maxn)}; {len(sa)} alleles in total). Chords: the {len(ch)} Release 2 ends in tier A/B alleles whose "
        f"outermost core block (>= 10 kb from the telomere boundary) is labelled neither home nor none, grouped by "
        f"recipient arm and by the best single CHM13 end of that block into {len(chords)} pairs; chord colour = that "
        f"best end, dot = recipient arm, width = fraction of the recipient arm's Release 2 ends with a structural-allele "
        f"call; dashed = mostly within-family ends ({nwf} pairs), i.e. the block is closer to another member of a "
        f"paralog family that contains the arm itself. Largest: {t4}. Blocks are labelled by paralog family in the "
        f"catalogue (e.g. PF04 = 3q/19p); the chord is drawn to the block's best single end. Labels are relative to "
        f"CHM13: the switch to home sequence lies inside shared paralogous sequence and is not necessarily a structural "
        f"breakpoint. (b) Nativity class of each arm's Release 2 ends (n = {min(comp[a]['n'] for a in order)}-"
        f"{max(comp[a]['n'] for a in order)} per arm; acroP pools five arms), from the end's own core structure: "
        f"native = only home and none blocks; within family = non-home blocks all in a family containing the arm; "
        f"foreign = at least one block from another family or arm. All Release 2 ends: native "
        f"{int(nat_tot.get('native', 0))}, within family {int(nat_tot.get('within_family', 0))}, foreign "
        f"{int(nat_tot.get('foreign', 0))}, no core structure {int(nat_tot.get('unassigned', 0))}. Arms with native "
        f"fraction < 0.5: {', '.join(mostly)}. Source: data/catalogue/alleles_structural.tsv, ends.tsv, "
        f"mosaic_blocks.tsv, paralog_families.tsv.")
    return over, {"n_arms": len(arms), "n_chords": len(chords), "n_chord_ends": len(ch), "within_family_pairs": nwf,
                  "mostly_nonnative": mostly, "top_chords": t4}


# ---------------------------------------------------------------- F2
F2_ARMS = ["20q", "3q", "6p", "8p", "16p", "11p"]


def spread(ys, mind, lo, hi):
    """Order-preserving label placement: forward pass enforces the minimum spacing, backward pass keeps labels within
    [lo, hi]; if the labels cannot fit, the spacing shrinks to fit."""
    n = len(ys)
    if n == 0:
        return []
    mind = min(mind, (hi - lo) / max(n - 1, 1))
    out = list(ys)
    out[0] = max(out[0], lo)
    for i in range(1, n):
        out[i] = max(out[i], out[i - 1] + mind)
    if out[-1] > hi:
        out[-1] = hi
        for i in range(n - 2, -1, -1):
            out[i] = min(out[i], out[i + 1] - mind)
    return out


def fig2(sa, ends, mb, captions, min_rows=10):
    ends_ok = ends[ends.set.isin(["HPRC_R2", "CHM13", "HG002_Q100"]) & ends.structural_allele_id.notna()]
    r2ids = set(ends.loc[ends.set == "HPRC_R2", "end_id"])
    raw = mb[mb.level == "raw"]
    rawg = {k: d for k, d in raw.groupby("end_id")}
    sai = sa.set_index("allele_id")
    BIN = 100
    NB = 100000 // BIN
    TEL = 60
    GAP = 6
    fig = plt.figure(figsize=(7.2, 8.6))
    W, H = 0.26, 0.235
    xs = [0.07, 0.565]
    ys = [0.685, 0.38, 0.075]
    axes = [fig.add_axes([xs[j], ys[i], W, H]) for i in range(3) for j in range(2)]
    used = set()
    stats = {}
    for ax, arm, letter in zip(axes, F2_ARMS, "abcdef"):
        d = ends_ok[ends_ok.catalogue_arm == arm]
        sizes = d.structural_allele_id.value_counts()
        shown = [a for a in sizes.index if sizes[a] >= min_rows]
        rest = [a for a in sizes.index if sizes[a] < min_rows]
        groups = [(a, d[d.structural_allele_id == a]) for a in shown]
        if rest:
            groups.append(("other", d[d.structural_allele_id.isin(rest)]))
        nrows = sum(len(g) for _, g in groups) + GAP * (len(groups) - 1)
        img = np.ones((nrows, TEL + NB, 3))
        img[:, :TEL - 8] = matplotlib.colors.to_rgb(TELO_COLOUR)
        row = 0
        labels = []

        def first_switch(e):
            b = rawg.get(e)
            return b.sub_end.iloc[0] if b is not None and len(b) > 1 else 1e9

        for gi, (aid, g) in enumerate(groups):
            if aid == "other":
                g = g.assign(_k=g.structural_allele_id.map(lambda x: sai.loc[x, "structure"]),
                             _s=g.end_id.map(first_switch)).sort_values(["_k", "_s"])
            else:
                g = g.assign(_s=g.end_id.map(first_switch)).sort_values("_s")
            r0 = row
            for e in g.end_id:
                b = rawg.get(e)
                if b is not None:
                    for s, en, lab, src in zip(b.sub_start, b.sub_end, b.label, b.best_end):
                        src = "none" if lab == "none" else (home_end(arm, src) if lab == "home" else src)
                        c0, c1 = int(max(0, s) // BIN), int(min(NB, math.ceil(min(en, 100000) / BIN)))
                        if c1 > c0:
                            img[row, TEL + c0:TEL + c1] = matplotlib.colors.to_rgb(arm_colour(src))
                            if isinstance(src, str) and src != "none":
                                used.add(src[3:])
                row += 1
            if aid == "other":
                txt = f"{len(rest)} other alleles ({len(g)} ends)"
                sw = []
            else:
                r = sai.loc[aid]
                txt = f"{aid.split('-')[-1]}  {r.structure}  {r.r2_freq:.2f}  {r.tier}"
                sw = [int(x) for x in str(r.switches_median).split(",")] if isinstance(r.switches_median, str) \
                    else ([int(r.switches_median)] if pd.notna(r.switches_median) else [])
            labels.append((r0, row, aid, txt, sw))
            row += GAP
        ax.imshow(img, aspect="auto", interpolation="nearest", extent=(-TEL * BIN / 1000, 100, nrows, 0),
                  rasterized=True)
        mind = nrows * 0.06
        lys = spread([(r0 + r1) / 2 for r0, r1, *_ in labels], mind, mind * 0.4, nrows - mind * 0.4)
        for (r0, r1, aid, txt, sw), ly in zip(labels, lys):
            ym = (r0 + r1) / 2
            # elbow leader: horizontal out of the band, then to the label; labels keep the band order, so no crossings
            ax.plot([100.5, 103.5, 108.0], [ym, ym, ly], color="#7a7a7a", lw=0.4, clip_on=False,
                    solid_joinstyle="miter")
            ax.text(108.6, ly, txt, fontsize=6, va="center", ha="left", clip_on=False)
            for s in sw:
                if 0 <= s <= 100000:
                    ax.plot(s / 1000, r0 - GAP / 2, marker="v", ms=2.6, color="#1c1c1c", mew=0, clip_on=False)
        ax.set_title(f"{arm}: {len(d)} ends in {len(sizes)} structural alleles", fontsize=8, loc="left", pad=4)
        panel_letter(ax, letter, x=-0.03, y=1.025)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlim(-TEL * BIN / 1000, 100)
        ax.set_ylim(nrows, -GAP)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xticklabels(["0", "25", "50", "75", "100"])
        ax.set_xlabel("kb from telomere boundary", labelpad=1)
        stats[arm] = {"ends": len(d), "alleles": len(sizes), "shown": len(shown), "other_alleles": len(rest),
                      "r2_ends": int(d.end_id.isin(r2ids).sum())}
    fig.text(0.07, 0.998, "Each row is one chromosome end, telomere at left (dark bar, not to scale); "
             "\u25bc = allele's median core switch.", fontsize=7, ha="left", va="top")
    fig.text(0.07, 0.982, "Labels: allele number, structure (telomere \u2192 centromere; home = own arm; "
             "PFnn = paralog family), Release 2 frequency, tier.", fontsize=7, ha="left", va="top")
    axk = fig.add_axes([0.05, -0.115, 0.80, 0.105])
    draw_palette_key(axk, used=used, title="Colour = CHM13 v2.0 end each block resembles most "
                                           "(faded: not present in this figure)")
    over = save(fig, "F2_allele_mosaics")
    captions["F2_allele_mosaics"] = (
        "**Figure 2 | Terminal mosaics of real chromosome ends, grouped by structural allele.** For six arms, every "
        "Release 2, CHM13 and HG002 Q100 end assigned to a structural allele of the arm is drawn as one row (pedigree "
        "validation ends are not drawn). Blocks are the raw mosaic blocks of mosaic_blocks.tsv (labels from per-base "
        "identity tracks on a 100 bp grid, base-level switch refinement), coloured by their best single CHM13 v2.0 end "
        "(column best_end; members of a paralog family share a hue and differ in lightness, see key); home blocks carry "
        "the arm's own colour even where a paralog is tied within the home margin (pooled groups: the member end); grey = "
        "none, no CHM13 end within 0.80 identity, in subtelomere coordinates (0 = telomere/subtelomere boundary; the telomere is "
        "drawn as a fixed-width dark bar, not to scale); white = beyond the end of a contig shorter than the 100 kb "
        "window. Ends are grouped by structural allele (defined on core blocks >= 10 kb; "
        f"alleles with >= {min_rows} ends are shown separately, the remainder pooled as 'other alleles' and sorted by "
        "structure) and ordered within a group by the first raw switch; leaders keep the order of the groups. Triangles above each group mark the allele's "
        "median core switch positions. Labels: allele number (SA-<arm>-NN), core structure (home, paralog family PFnn, "
        "CHM13 arm or none per core block), Release 2 frequency (r2_freq) and tier. Ends per arm: " + ", ".join(f"{a} {s['ends']} ({s['alleles']} alleles)" for a, s in
                                                         stats.items()) +
        ". Block labels are relative to CHM13: inside paralogous subtelomeric sequence a colour change is a change in "
        "which CHM13 end is most similar, not necessarily a structural breakpoint; short raw blocks inside a group "
        "(< 10 kb) do not change the allele call. Source: data/catalogue/ends.tsv, mosaic_blocks.tsv, "
        "alleles_structural.tsv.")
    return over, stats


# ---------------------------------------------------------------- F3
def common_forms(sa, min_freq=0.05):
    c = sa[(sa.nativity_class != "native") & (sa.tier == "A") & (sa.r2_freq >= min_freq)].copy()
    return c.sort_values("r2_freq", ascending=False)


def row_label(r):
    return f"{r.allele_id}  {r.structure}"


def fig3(sa, fr, captions):
    c = common_forms(sa)
    f = fr[(fr.allele_type == "structural_allele") & fr.allele_id.isin(c.allele_id)]
    piv = f.pivot_table(index="allele_id", columns="superpopulation", values="frequency").reindex(c.allele_id)
    nass = f.pivot_table(index="allele_id", columns="superpopulation", values="n_arm_ends_assessed")
    piv = piv.reindex(columns=SUPERPOPS + ["ALL"]).fillna(0.0)
    M = piv[SUPERPOPS].values
    n = len(c)
    fig = plt.figure(figsize=(5.2, 0.155 * n + 1.3))
    ax = fig.add_axes([0.42, 0.08, 0.38, 0.84])
    axm = fig.add_axes([0.805, 0.08, 0.012, 0.84], sharey=None)
    ax2 = fig.add_axes([0.835, 0.08, 0.13, 0.84], sharey=ax)
    vmax = max(0.8, float(np.nanmax(M)))
    im = ax.imshow(M, aspect="auto", cmap="viridis_r" if False else "YlGnBu", vmin=0, vmax=vmax,
                   interpolation="nearest")
    for i in range(n):
        for j in range(len(SUPERPOPS)):
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}".lstrip("0") if v > 0 else "0", ha="center", va="center", fontsize=6,
                    color="white" if v > 0.42 else "#1c1c1c")
    ax.set_xticks(range(len(SUPERPOPS)))
    ax.set_xticklabels(SUPERPOPS)
    ax.xaxis.tick_top()
    ax.set_yticks(range(n))
    ax.set_yticklabels([row_label(r) for r in c.itertuples()], fontsize=6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    # AFR minus non-AFR difference marker as a lollipop of ALL frequency
    for i, r in enumerate(c.itertuples()):
        axm.add_patch(Rectangle((0, i - 0.4), 1, 0.8, facecolor=NATIVITY_COLOURS[r.nativity_class], edgecolor="none"))
    axm.set_xlim(0, 1)
    axm.set_ylim(n - 0.5, -0.5)
    axm.set_axis_off()
    axm.text(0.5, -0.9, "class", fontsize=6, ha="center", va="bottom", rotation=90)
    nwf3 = int((c.nativity_class == "within_family").sum())
    fig.legend(handles=[Patch(facecolor=NATIVITY_COLOURS["foreign"], label="foreign"),
                        Patch(facecolor=NATIVITY_COLOURS["within_family"], label="within family")],
               loc="lower left", bbox_to_anchor=(0.02, 0.005), ncol=2, fontsize=6, handlelength=1.0)
    allf = piv["ALL"].values
    ax2.hlines(range(n), 0, allf, color="#9a9a9a", lw=0.6)
    ax2.plot(allf, range(n), "o", color="#1c1c1c", ms=2.5)
    ax2.set_xlim(0, max(0.52, float(np.nanmax(piv["ALL"].values)) * 1.08))
    ax2.set_xticks([0, 0.25, 0.5])
    ax2.set_xticklabels(["0", ".25", ".5"])
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position("top")
    ax2.set_xlabel("all Release 2", fontsize=6, labelpad=2)
    ax2.tick_params(axis="y", left=False, labelleft=False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["top"].set_visible(True)
    ax2.set_ylim(n - 0.5, -0.5)
    cax = fig.add_axes([0.42, 0.035, 0.38, 0.012])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Frequency among the arm's Release 2 ends", fontsize=6, labelpad=1)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8])
    afr_diff = (piv["AFR"] - piv[["AMR", "EAS", "EUR", "SAS"]].mean(axis=1))
    in_all5 = int((M > 0).all(axis=1).sum())
    fig.text(0.02, 0.975, f"{in_all5} of {n} common non-native alleles are seen in all five superpopulations",
             fontsize=8, ha="left", va="top")
    over = save(fig, "F3_population_frequencies")
    big = afr_diff.abs().sort_values(ascending=False).head(3)
    lo = int(nass[SUPERPOPS].min().min())
    hi = int(nass[SUPERPOPS].max().max())
    captions["F3_population_frequencies"] = (
        f"**Figure 3 | Frequencies of the common non-native structural alleles by 1000 Genomes superpopulation.** "
        f"Rows: the {n} tier A structural alleles of nativity class foreign or within family with Release 2 frequency "
        f">= 0.05 ({n - nwf3} foreign, {nwf3} within family; class strip left of the cells), with allele ID and core "
        f"structure (telomere -> centromere; 'home' = the arm's own sequence, PFnn = paralog family, none = no CHM13 "
        f"end within 0.80 identity), ordered by overall frequency. Cells: "
        f"carrier ends / Release 2 ends of the arm with a structural-allele call in that superpopulation "
        f"(denominators {lo}-{hi} ends per arm and superpopulation; allele_frequencies_by_superpop.tsv); right: "
        f"all Release 2 ends. {in_all5} of {n} alleles are observed in all five superpopulations. Largest AFR vs "
        f"mean non-AFR differences: " + "; ".join(f"{a} {afr_diff[a]:+.2f}" for a in big.index) +
        ". Frequencies are descriptive; AMR and part of AFR (ACB, ASW) are admixed; formal tests of population "
        "differentiation for terminal forms are in analyses/allele_ages.")
    return over, {"n_common": n, "in_all5": in_all5, "denoms": (lo, hi)}


# ---------------------------------------------------------------- F4
EVID = ["recurrent", "verkko", "reads", "pedigree", "t2t_assembly"]
EVID_LABEL = {"recurrent": "recurrent (>= 3 samples)", "verkko": "Verkko concordant", "reads": "reads span switch",
              "pedigree": "pedigree transmission", "t2t_assembly": "in T2T assembly"}


def fig4(sa, hg, ends, captions):
    tier_col = {"A": "#2b2b2b", "B": "#7a7a7a", "C": "#c4c4c4"}
    fig = plt.figure(figsize=(7.2, 5.8))
    ax_a = fig.add_axes([0.10, 0.70, 0.15, 0.22])
    ax_a2 = fig.add_axes([0.29, 0.70, 0.15, 0.22])
    ax_bt = fig.add_axes([0.20, 0.445, 0.33, 0.07])
    ax_bb = fig.add_axes([0.20, 0.33, 0.33, 0.105], sharex=ax_bt)
    ax_bm = fig.add_axes([0.20, 0.05, 0.33, 0.27], sharex=ax_bt)
    ax_c = fig.add_axes([0.775, 0.05, 0.215, 0.725])
    cats = ["structural\nalleles", "haplogroups"]
    for ax, val, xl in [(ax_a, "count", "Number"), (ax_a2, "n_r2_ends", "Release 2 ends")]:
        for i, tab in enumerate([sa, hg]):
            left = 0
            for tr in "ABC":
                v = (tab.tier == tr).sum() if val == "count" else tab.loc[tab.tier == tr, "n_r2_ends"].sum()
                ax.barh(i, v, left=left, color=tier_col[tr], height=0.6)
                left += v
        ax.set_yticks([0, 1])
        ax.set_yticklabels(cats if ax is ax_a else ["", ""])
        ax.set_ylim(1.5, -0.5)
        ax.set_xlabel(xl, labelpad=2)
        ax.tick_params(axis="y", length=0)
    ax_a2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: f"{x / 1000:g}k" if x else "0"))
    ax_a.set_title("Tier C is numerous but holds few ends", fontsize=8, loc="left")
    panel_letter(ax_a, "a", x=-0.62, y=1.04)
    ax_a2.legend(handles=[Patch(color=tier_col[t], label=f"tier {t}") for t in "ABC"], loc="upper left",
                 fontsize=6, handlelength=1.0, bbox_to_anchor=(1.03, 1.0))
    # (b) UpSet with a broken linear count axis
    ev = sa.evidence.fillna("").map(lambda s: tuple(e for e in EVID if e in s.split(",")))
    combos = Counter(zip(ev, sa.tier))
    tot = Counter(ev)
    keys = [k for k in sorted(tot, key=lambda k: -tot[k]) if tot[k] >= 3]
    x = np.arange(len(keys))
    for ax in (ax_bt, ax_bb):
        bottom = np.zeros(len(keys))
        for tr in "ABC":
            v = np.array([combos.get((k, tr), 0) for k in keys])
            ax.bar(x, v, bottom=bottom, color=tier_col[tr], width=0.7)
            bottom += v
    hi = max(tot.values())
    second = sorted(tot.values())[-2]
    ax_bb.set_ylim(0, 62)
    ax_bt.set_ylim(second - 20 if second > 100 else 180, hi * 1.12)
    for xi, k in zip(x, keys):
        axk_ = ax_bt if tot[k] > 62 else ax_bb
        axk_.text(xi, tot[k] + (hi * 0.02 if axk_ is ax_bt else 1.5), f"{tot[k]}", ha="center", va="bottom",
                  fontsize=6)
    ax_bt.spines["bottom"].set_visible(False)
    ax_bt.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bb.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bb.spines["bottom"].set_visible(False)
    ax_bt.set_yticks([int(round(second / 100.0) * 100), int(hi // 100 * 100)])
    ax_bb.set_yticks([0, 25, 50])
    kw = dict(color="#1c1c1c", clip_on=False, lw=0.6)
    for ax, yy in [(ax_bt, 0), (ax_bb, 1)]:
        ax.plot([-0.012, 0.012], [yy - 0.03, yy + 0.03], transform=ax.transAxes, **kw)
    ax_bb.set_ylabel("Alleles", labelpad=2)
    ax_bb.yaxis.set_label_coords(-0.13, 0.9)
    ax_bt.set_title("Most alleles lack an orthogonal line of evidence",
                    fontsize=8, loc="left", x=-0.38, pad=8)
    panel_letter(ax_bt, "b", x=-0.40, y=1.25)
    for j, e in enumerate(EVID):
        for xi, k in zip(x, keys):
            ax_bm.plot(xi, j, "o", ms=3.4, color="#1c1c1c" if e in k else "#e0e0e0")
    for xi, k in zip(x, keys):
        js = [j for j, e in enumerate(EVID) if e in k]
        if len(js) > 1:
            ax_bm.plot([xi, xi], [min(js), max(js)], color="#1c1c1c", lw=0.8)
    ax_bm.set_yticks(range(len(EVID)))
    ax_bm.set_yticklabels([EVID_LABEL[e] for e in EVID])
    ax_bm.set_ylim(len(EVID) - 0.5, -0.5)
    ax_bm.set_xlim(-0.6, len(keys) - 0.4)
    ax_bm.set_xticks([])
    for s in ["bottom", "left"]:
        ax_bm.spines[s].set_visible(False)
    ax_bm.tick_params(axis="y", length=0)
    ax_bm.text(0, len(EVID) - 0.35, "no line", fontsize=6, ha="center", va="top")
    hidden = sum(tot[k] for k in tot if tot[k] < 3)
    # (c)
    c = common_forms(sa)
    cols = [("n_samples", "samples"), ("n_vk_concordant", "Verkko concordant"),
            ("n_vk_discordant", "Verkko discordant"), ("n_layout_supported", "read layout"),
            ("n_hg002_reads_supported", "HG002 reads"), ("n_pedigree_intact", "pedigree intact"),
            ("_t2t", "T2T members")]
    c = c.assign(_t2t=c.n_chm13 + c.n_hg002_q100)
    n = len(c)
    for j, (col, lab) in enumerate(cols):
        v = c[col].values
        for i, vi in enumerate(v):
            if col == "n_samples":
                ax_c.text(j + 0.25, i, f"{int(vi)}", fontsize=6, va="center", ha="right")
            elif vi > 0:
                s = 1.8 + 2.2 * math.sqrt(vi / (c[col].max() or 1))
                colr = "#d62728" if col == "n_vk_discordant" else "#1c1c1c"
                ax_c.plot(j, i, "o", ms=s, color=colr, alpha=0.85)
                ax_c.text(j + 0.22, i, f"{int(vi)}", fontsize=6, va="center", ha="left")
            else:
                ax_c.plot(j, i, "_", ms=3, color="#c4c4c4")
    ax_c.set_xticks(range(len(cols)))
    ax_c.set_xticklabels([l for _, l in cols], rotation=90, fontsize=6)
    ax_c.xaxis.tick_top()
    ax_c.set_yticks(range(n))
    ax_c.set_yticklabels([row_label(r) for r in c.itertuples()], fontsize=6)
    ax_c.set_ylim(n - 0.5, -0.7)
    ax_c.set_xlim(-0.5, len(cols) - 0.1)
    for s in ax_c.spines.values():
        s.set_visible(False)
    ax_c.tick_params(length=0)
    fig.text(0.60, 0.985, f"Members with each evidence type,\n{n} common non-native alleles", fontsize=8,
             ha="left", va="top")
    fig.text(0.585, 0.985, "c", fontsize=10, fontweight="bold", ha="right", va="top")
    over = save(fig, "F4_validation")
    nA = int((sa.tier == "A").sum())
    nB = int((sa.tier == "B").sum())
    nC = int((sa.tier == "C").sum())
    endsA = int(sa.loc[sa.tier == "A", "n_r2_ends"].sum())
    endsAB = int(sa.loc[sa.tier.isin(["A", "B"]), "n_r2_ends"].sum())
    endsall = int(sa.n_r2_ends.sum())
    nsing = int(sa.singleton.sum())
    vkd = int(c.n_vk_discordant.sum())
    vkc = int(c.n_vk_concordant.sum())
    nnone = tot.get((), 0)
    nrec = tot.get(("recurrent",), 0)
    captions["F4_validation"] = (
        f"**Figure 4 | Evidence tiers and lines of evidence.** (a) Structural alleles ({len(sa)}: tier A {nA}, B {nB}, "
        f"C {nC}; {nsing} singletons) and sequence haplogroups ({len(hg)}; A {int((hg.tier == 'A').sum())}, "
        f"B {int((hg.tier == 'B').sum())}, C {int((hg.tier == 'C').sum())}), as counts and as the number of Release 2 "
        f"ends they contain: tier A alleles hold {endsA} and tier A/B {endsAB} of the {endsall} Release 2 ends with a "
        f"structural-allele call ({endsAB / endsall:.1%}). (b) Combinations of evidence lines per structural allele "
        f"(UpSet layout; bars stacked by tier on a broken linear axis; {nnone} alleles have no evidence line and "
        f"{nrec} only recurrence; combinations with < 3 alleles, {hidden} alleles in total, are omitted). Tiers: A = "
        f"recurrent (>= 3 unrelated samples) and >= 1 orthogonal line, or >= 2 orthogonal lines; B = recurrent without "
        f"an orthogonal line, or 2 samples with >= 1 orthogonal line; C = all others. Orthogonal lines: a member's "
        f"Verkko assembly agrees; HG002 Q100 ONT reads cross the switch or a hifiasm read layout supports the end; "
        f"intact parent-to-child transmission in the WashU pedigree; a member in HG002 Q100 (or CHM13 for non-native "
        f"alleles). (c) For the {n} common non-native alleles of Fig. 3, members with each evidence status (dot area "
        f"scaled within column; number = count; first column = samples; red = Verkko-discordant members, {vkd} in total against {vkc} "
        f"concordant). A missing orthogonal line mostly reflects that few ends were testable (Verkko assemblies and "
        f"read layouts exist for a subset of samples; HG002 and the pedigree are single families), not contrary "
        f"evidence. Source: data/catalogue/alleles_structural.tsv, haplogroups.tsv.")
    return over, {"tierA": nA, "tierB": nB, "tierC": nC, "endsAB_frac": endsAB / endsall}


# ---------------------------------------------------------------- palette figure
def fig_palette(captions):
    fig = plt.figure(figsize=(7.2, 1.45))
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.9])
    draw_palette_key(ax, used=None, title="CHM13 v2.0 end palette used in all atlas figures and on the website "
                                          "(paralog-family members share a hue)")
    over = save(fig, "palette_chm13_ends")
    captions["palette_chm13_ends"] = (
        "**Palette.** One colour per CHM13 v2.0 chromosome end (48), produced by pipeline/atlas_palette.py. Members of a "
        "paralog family (data/catalogue/paralog_families.tsv; outlined in the key) share one hue and differ in "
        "lightness, so a family-labelled block reads as its family whichever member is its best single end; every other "
        "end has its own colour. Colours were chosen once (seeded optimisation) to maximise the CIEDE2000 difference, in "
        "normal vision and in simulated deuteranopia (Machado et al. 2009, doi:10.1109/TVCG.2009.113), between CHM13 ends that are adjacent or "
        "co-occur in the catalogue mosaics; pairwise values in figures/atlas/palette_qc.tsv. " + PALETTE_QC_TEXT +
        " Grey = none (no CHM13 end within 0.80 identity). Hex codes: figures/atlas/PALETTE.tsv.")
    return over


def main():
    style()
    os.makedirs(OUT, exist_ok=True)
    sa, hg, fr, ends, mb = load()
    captions = {}
    report = {}
    report["palette"] = fig_palette(captions)
    o1, s1 = fig1(sa, ends, mb, captions)
    o2, s2 = fig2(sa, ends, mb, captions)
    o3, s3 = fig3(sa, fr, captions)
    o4, s4 = fig4(sa, hg, ends, captions)
    report.update({"F1": o1, "F2": o2, "F3": o3, "F4": o4})
    with open(os.path.join(OUT, "CAPTIONS.md"), "w") as fh:
        fh.write("# Atlas figures - captions\n\nGenerated by `pipeline/figures_atlas.py` from the catalogue v0.1 release "
                 "tables. Each figure is provided as PNG (300 dpi), PDF and SVG.\n\n")
        for k in ["F1_genome_overview", "F2_allele_mosaics", "F3_population_frequencies", "F4_validation",
                  "palette_chm13_ends"]:
            fh.write(f"## {k}\n\n![{k}]({k}.png)\n\n{captions[k]}\n\n")
    with open(os.path.join(OUT, "figure_stats.json"), "w") as fh:
        json.dump({"F1": s1, "F2": s2, "F3": s3, "F4": s4}, fh, indent=1, default=str)
    for k, v in report.items():
        print(k, "text overlaps:", len(v), v[:6])
    return captions


if __name__ == "__main__":
    main()
