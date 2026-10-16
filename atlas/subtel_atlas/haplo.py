"""Sequence haplogroups (catalogue step 2e; gc_ident is used by pipeline/catalogue_07_haplogroups.py and the query tool).

Compared sequence: the first HG_LEN (80 kb) of the subtelomere (window minus the telomere prefix s0). Similarity: Mash-style ANI of
FracMinHash sketches (sketch.py; k = 21, scale 1/10). Catalogue haplogroups are UPGMA clusters cut at ANI 0.995. The query tool
compares an end with the medoid of every represented haplogroup of its arm (representatives.fa.gz: haplogroups seen in >= 3
samples) over the same prefix length and checks the best one by alignment (gc_ident).
"""
import re
import numpy as np
from . import sketch as S

HG_LEN, K, SCALE, HG_ANI = 80_000, 21, 10, 0.995


def gc_ident(qs, ts):
    """gap-compressed identity (1 - de) over primary hits, weighted by aligned query length, and aligned fraction"""
    import mappy
    a = mappy.Aligner(seq=ts, preset='asm20')
    cov = np.zeros(len(qs), bool); w = d = 0.0
    for h in a.map(qs, cs=True):
        if not h.is_primary: continue
        cov[h.q_st:h.q_en] = True
        sub = len(re.findall(r'\*', h.cs)); ins = len(re.findall(r'\+', h.cs)); dl = len(re.findall(r'-', h.cs))
        mt = sum(int(x) for x in re.findall(r':(\d+)', h.cs)); L = mt + sub + ins + dl
        w += L; d += sub + ins + dl
    return (1 - d / w if w else np.nan), float(cov.mean())


def ani(ska, skb):
    inter = len(np.intersect1d(ska, skb, assume_unique=True))
    return float(S.ani_from_jaccard(inter / max(len(ska) + len(skb) - inter, 1), K))


def nearest(sub, medoids):
    """sub: query subtelomere (<= HG_LEN). medoids: list of (haplogroup_id, medoid subtelomere). Compared over the query's length.
    Returns list of (ani, haplogroup_id) sorted best first (ties by ID)."""
    Lp = len(sub); ske = S.sketch(sub, K, SCALE); out = []
    for hid, ms in medoids:
        out.append((ani(ske, S.sketch(ms[:Lp], K, SCALE)), hid))
    return sorted(out, key=lambda x: (-x[0], x[1]))
