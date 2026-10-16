#!/usr/bin/env python3
"""Family-aware colour palette for the 48 CHM13 v2.0 ends used in the atlas figures and on the website.

Usage (from the repository root):  python pipeline/atlas_palette.py      # writes figures/atlas/PALETTE.tsv, palette_qc.tsv

Design
  * Paralog families (data/catalogue/paralog_families.tsv): all members of a family share one hue and differ in
    lightness, so a family-labelled block reads as that family whichever member is its best single end.
  * Every other CHM13 end gets its own colour from an in-gamut CIE LCh candidate pool.
  * Hues of families and colours of single ends are chosen by seeded simulated annealing to maximise the CIEDE2000
    colour difference, in normal vision and in a simulated deuteranopia (Machado et al. 2009, severity 1), between
    CHM13 ends that are adjacent blocks or co-occur in the same end in the catalogue mosaics (mosaic_blocks.tsv, raw
    blocks, best_end; the end's home arm included), with a smaller weight on all other pairs.
  * palette_qc.tsv lists dE2000 (normal and deuteranopia) for every pair with its co-occurrence weight.
The optimisation is deterministic (seed below). The palette is fixed once PALETTE.tsv exists; re-run with --force.
"""
import math
import os
import random
import sys
from collections import Counter

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CAT = os.path.join(ROOT, "data", "catalogue")
OUT = os.path.join(ROOT, "figures", "atlas")
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
ENDS = [c + s for c in CHROMS for s in "pq"]
POOLED_HOME = {"acroP": "15p", "PAR1": "Xp", "PAR2": "Xq"}
SEED = 20261016
NONE_HEX = "#ececec"

# ------------------------------------------------------------------ colour science
M_RGB2XYZ = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750],
                      [0.0193339, 0.1191920, 0.9503041]])
WHITE = np.array([0.95047, 1.0, 1.08883])
# Machado, Oliveira & Fernandes 2009, IEEE TVCG 15:1291-1298, doi:10.1109/TVCG.2009.113;
# deuteranomaly severity 1.0 (linear RGB)
M_DEUT = np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                   [-0.011820, 0.042940, 0.968881]])


def lch_to_lab(L, C, h):
    return np.array([L, C * math.cos(math.radians(h)), C * math.sin(math.radians(h))])


def lab_to_linrgb(lab):
    L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]
    fy = (L + 16) / 116
    fx, fz = fy + a / 500, fy - b / 200

    def finv(t):
        return np.where(t ** 3 > 0.008856, t ** 3, (t - 16 / 116) / 7.787)

    xyz = np.stack([finv(fx) * WHITE[0], finv(fy) * WHITE[1], finv(fz) * WHITE[2]], axis=-1)
    return xyz @ np.linalg.inv(M_RGB2XYZ).T


def linrgb_to_lab(rgb):
    xyz = rgb @ M_RGB2XYZ.T / WHITE

    def f(t):
        return np.where(t > 0.008856, np.cbrt(t), 7.787 * t + 16 / 116)

    fx, fy, fz = f(xyz[..., 0]), f(xyz[..., 1]), f(xyz[..., 2])
    return np.stack([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)], axis=-1)


def srgb_gamma(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)


def srgb_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def hex_to_lab(h):
    rgb = np.array([int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return linrgb_to_lab(srgb_linear(rgb))


def lab_to_hex(lab):
    rgb = srgb_gamma(lab_to_linrgb(np.asarray(lab, float)))
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in rgb)


def in_gamut(lab):
    rgb = lab_to_linrgb(np.asarray(lab, float))
    return bool(np.all(rgb >= -1e-6) and np.all(rgb <= 1 + 1e-6))


def deuteranopia_lab(lab):
    """Simulated deuteranope appearance (on the displayed, gamut-clipped colour)."""
    rgb = np.clip(lab_to_linrgb(np.asarray(lab, float)), 0, 1)
    sim = np.clip(rgb @ M_DEUT.T, 0, 1)
    return linrgb_to_lab(sim)


def ciede2000(lab1, lab2):
    """Vectorised CIEDE2000 (Sharma et al. 2005). lab1, lab2: (..., 3)."""
    L1, a1, b1 = lab1[..., 0], lab1[..., 1], lab1[..., 2]
    L2, a2, b2 = lab2[..., 0], lab2[..., 1], lab2[..., 2]
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dLp = L2 - L1
    dCp = C2p - C1p
    dh = h2p - h1p
    dh = np.where(C1p * C2p == 0, 0, np.where(dh > 180, dh - 360, np.where(dh < -180, dh + 360, dh)))
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dh / 2))
    Lbp = (L1 + L2) / 2
    Cbp = (C1p + C2p) / 2
    hsum = h1p + h2p
    hbp = np.where(C1p * C2p == 0, hsum,
                   np.where(np.abs(h1p - h2p) <= 180, hsum / 2, np.where(hsum < 360, (hsum + 360) / 2, (hsum - 360) / 2)))
    T = (1 - 0.17 * np.cos(np.radians(hbp - 30)) + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6)) - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    dtheta = 30 * np.exp(-((hbp - 275) / 25) ** 2)
    Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    Sl = 1 + 0.015 * (Lbp - 50) ** 2 / np.sqrt(20 + (Lbp - 50) ** 2)
    Sc = 1 + 0.045 * Cbp
    Sh = 1 + 0.015 * Cbp * T
    Rt = -np.sin(np.radians(2 * dtheta)) * Rc
    return np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2 + Rt * (dCp / Sc) * (dHp / Sh))


# ------------------------------------------------------------------ inputs
def families():
    pf = pd.read_csv(os.path.join(CAT, "paralog_families.tsv"), sep="\t")
    return {r.family: r.members.split(",") for r in pf.itertuples()}


def weights():
    mb = pd.read_csv(os.path.join(CAT, "mosaic_blocks.tsv"), sep="\t", low_memory=False,
                     usecols=["end_id", "level", "block", "best_end"])
    ends = pd.read_csv(os.path.join(CAT, "ends.tsv"), sep="\t", low_memory=False, usecols=["end_id", "catalogue_arm"])
    arm = dict(zip(ends.end_id, ends.catalogue_arm))
    raw = mb[mb.level == "raw"].sort_values(["end_id", "block"])
    W = Counter()
    for e, d in raw.groupby("end_id", sort=False):
        s = [x[3:] for x in d.best_end if isinstance(x, str) and x.startswith("chr")]
        for a, b in zip(s, s[1:]):
            if a != b:
                W[tuple(sorted((a, b)))] += 1.0
        home = arm.get(e)
        home = POOLED_HOME.get(home, home)
        u = set(s) | ({home} if home in ENDS else set())
        u = sorted(u)
        for i in range(len(u)):
            for j in range(i + 1, len(u)):
                W[(u[i], u[j])] += 0.5
    return W


# ------------------------------------------------------------------ optimisation
def family_shades(n):
    return {1: [55], 2: [42, 74], 3: [36, 57, 78], 4: [34, 50, 66, 82], 5: [32, 45, 58, 71, 84]}[n]


def build_pool():
    pool = []
    for L in [32, 40, 48, 56, 64, 72, 80, 88]:
        for C in [28, 40, 52, 64]:
            for h in range(0, 360, 12):
                lab = lch_to_lab(L, C, h)
                if in_gamut(lab):
                    pool.append(lab)
    return np.array(pool)


def optimise(force=False):
    fam = families()
    W = weights()
    fam_of = {m: f for f, ms in fam.items() for m in ms}
    singles = [e for e in ENDS if e not in fam_of]
    idx = {e: i for i, e in enumerate(ENDS)}
    n = len(ENDS)
    Wm = np.zeros((n, n))
    for (a, b), w in W.items():
        if a in idx and b in idx:
            Wm[idx[a], idx[b]] = Wm[idx[b], idx[a]] = w
    base = Wm[Wm > 0].mean() * 0.02  # every pair matters a little (legend distinctness)
    Wt = Wm + base
    same_fam = np.zeros((n, n), bool)
    for f, ms in fam.items():
        for a in ms:
            for b in ms:
                same_fam[idx[a], idx[b]] = True
    np.fill_diagonal(Wt, 0)
    Wt[same_fam] = 0
    pool = build_pool()
    rng = random.Random(SEED)
    fam_names = sorted(fam)

    def fam_labs(fh):
        out = {}
        for f in fam_names:
            h, C = fh[f]
            for m, L in zip(fam[f], family_shades(len(fam[f]))):
                c = C
                lab = lch_to_lab(L, c, h)
                while not in_gamut(lab) and c > 8:
                    c -= 2
                    lab = lch_to_lab(L, c, h)
                out[m] = lab
        return out

    def assemble(fh, sing):
        labs = fam_labs(fh)
        for e, k in zip(singles, sing):
            labs[e] = pool[k]
        return np.array([labs[e] for e in ENDS])

    def cost(L):
        D = deuteranopia_lab(L)
        dn = ciede2000(L[:, None, :], L[None, :, :])
        dd = ciede2000(D[:, None, :], D[None, :, :])
        d = np.minimum(dn, dd)
        pen = np.exp(-d / 10.0)
        off = ~same_fam & ~np.eye(n, dtype=bool)
        # floors for every pair of different families/ends: distinct in a legend (normal) and not identical (deutan)
        floor = (20.0 * np.clip(10.0 - dn, 0, None) + 5.0 * np.clip(5.0 - dd, 0, None))[off].sum() / 2
        # within-family shades must stay distinguishable
        fam_pen = np.where(same_fam & ~np.eye(n, dtype=bool) & (dn < 12), 50.0, 0.0).sum() + floor
        # nothing may sit too close to the 'none' grey or to white/black text
        grey = hex_to_lab(NONE_HEX)
        dg = ciede2000(L, np.broadcast_to(grey, L.shape))
        return float((Wt * pen).sum() / 2 + fam_pen + 200.0 * np.sum(dg < 12))

    fh = {f: (rng.randrange(0, 360, 12), 52) for f in fam_names}
    sing = rng.sample(range(len(pool)), len(singles))
    L = assemble(fh, sing)
    c = cost(L)
    best = (c, dict(fh), list(sing))
    T = c * 0.05
    for it in range(45000):
        fh2, sing2 = dict(fh), list(sing)
        r = rng.random()
        if r < 0.15:
            f = rng.choice(fam_names)
            fh2[f] = ((fh2[f][0] + rng.choice([-24, -12, 12, 24])) % 360, rng.choice([40, 48, 56, 64]))
        elif r < 0.55:
            i = rng.randrange(len(singles))
            k = rng.randrange(len(pool))
            if k in sing2:
                continue
            sing2[i] = k
        else:
            i, j = rng.sample(range(len(singles)), 2)
            sing2[i], sing2[j] = sing2[j], sing2[i]
        c2 = cost(assemble(fh2, sing2))
        if c2 < c or rng.random() < math.exp((c - c2) / max(T, 1e-9)):
            fh, sing, c = fh2, sing2, c2
            if c < best[0]:
                best = (c, dict(fh2), list(sing2))
        T *= 0.99985
    L = assemble(best[1], best[2])
    return L, fam, fam_of, Wm, best[0]


def qc_table(L, fam_of, Wm):
    D = deuteranopia_lab(L)
    dn = ciede2000(L[:, None, :], L[None, :, :])
    dd = ciede2000(D[:, None, :], D[None, :, :])
    rows = []
    for i in range(len(ENDS)):
        for j in range(i + 1, len(ENDS)):
            a, b = ENDS[i], ENDS[j]
            rows.append({"end_a": a, "end_b": b, "cooccurrence_weight": round(Wm[i, j], 1),
                         "same_family": bool(fam_of.get(a) and fam_of.get(a) == fam_of.get(b)),
                         "dE2000": round(float(dn[i, j]), 2), "dE2000_deuteranopia": round(float(dd[i, j]), 2)})
    return pd.DataFrame(rows).sort_values("cooccurrence_weight", ascending=False)


def main(force=False):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, "PALETTE.tsv")
    if os.path.exists(p) and not force:
        old = pd.read_csv(p, sep="\t")
        if "family" in old.columns:
            print("PALETTE.tsv exists (family-aware); use --force to re-optimise")
            return
    L, fam, fam_of, Wm, c = optimise()
    hexes = [lab_to_hex(l) for l in L]
    rows = []
    for e, h, l in zip(ENDS, hexes, L):
        Lv, a, b = l
        rows.append({"chm13_end": "chr" + e, "arm": e, "hex": h, "family": fam_of.get(e, ""),
                     "L": round(float(Lv), 1), "C": round(float(math.hypot(a, b)), 1),
                     "h": round(math.degrees(math.atan2(b, a)) % 360, 1),
                     "rule": ("family hue, member shade" if e in fam_of else "own colour")})
    rows.append({"chm13_end": "none", "arm": "none", "hex": NONE_HEX, "family": "", "L": None, "C": None, "h": None,
                 "rule": "no CHM13 end within 0.80 identity"})
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    q = qc_table(np.array([hex_to_lab(h) for h in hexes]), fam_of, Wm)
    q.to_csv(os.path.join(OUT, "palette_qc.tsv"), sep="\t", index=False)
    top = q[(q.cooccurrence_weight >= 100) & ~q.same_family]
    print(f"cost {c:.1f}; pairs with weight >= 100 (between families/ends): {len(top)}; "
          f"min dE2000 {top.dE2000.min():.1f}, min deuteranopia {top.dE2000_deuteranopia.min():.1f}")


if __name__ == "__main__":
    main(force="--force" in sys.argv)
