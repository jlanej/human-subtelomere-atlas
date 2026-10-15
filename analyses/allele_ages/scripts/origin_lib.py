"""Shared helpers for the single-origin / allele-age analysis.

Coordinates: every sequence (HPRC tip window, CHM13 end, ape end) is oriented telomere-first, so a position is the
distance (bp, 0-based) from the chromosome terminus. HPRC L tips and CHM13/ape p-type ends are already telomere-first;
R tips / q ends are reverse-complemented.
"""
import os, re, glob
import numpy as np
import mappy as mp
import pysam

IN = os.environ.get('ATLAS_INPUT_DIR', '')
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
WORK = os.path.join(REPO, 'work', 'allele_ages')
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def rc(s):
    return s.translate(COMP)[::-1]


class TipStore:
    """Telomere-first 100 kb windows of named Release 2 tips, from the analysis-home fetch (foreign/tie tips) and this
    track's fetch of the remaining tips (work/allele_ages/r2_100k_native)."""
    def __init__(self):
        self.idx = {}
        for d in (os.path.join(IN, 'subtel', 'foreign', 'r2_100k'), os.path.join(WORK, 'r2_100k_native')):
            for f in glob.glob(os.path.join(d, '*.fa')):
                if not os.path.exists(f[:-3] + '.ok') or not os.path.exists(f + '.fai'):
                    if os.path.exists(f[:-3] + '.ok') and not os.path.exists(f + '.fai'):
                        pysam.faidx(f)
                    else:
                        continue
                asm = os.path.basename(f)[:-3]
                for line in open(f + '.fai'):
                    name = line.split('\t')[0]
                    c, rg = name.rsplit(':', 1); s, e = map(int, rg.split('-'))
                    self.idx[(asm, c, 'L' if s == 1 else 'R')] = (f, name)
        self._fh = {}

    def get(self, asm, contig, end):
        k = (asm, contig, end)
        if k not in self.idx:
            return None
        f, name = self.idx[k]
        if f not in self._fh:
            self._fh[f] = pysam.FastaFile(f)
        s = self._fh[f].fetch(name).upper()
        return s if end == 'L' else rc(s)


def chm13_end(name):
    """CHM13 v2.0 outer 500 kb of an end ('chr3q'), telomere-first."""
    fa = pysam.FastaFile(os.path.join(IN, 'subtel', 'foreign', 'ref', 'chm13_ends500k.fa'))
    s = fa.fetch(name).upper()
    return s if name.endswith('p') else rc(s)


TELO = re.compile(r'(?:CCCTAA){3,}')


def telo_len(s, maxscan=30000):
    """Length of the terminal (CCCTAA)n block of a telomere-first sequence: end of the last run of >=3 canonical repeats
    that lies within 1 kb of the previous run, starting within the first 1 kb."""
    end = 0
    for m in TELO.finditer(s[:maxscan]):
        if m.start() - end > 1000:
            break
        end = m.end()
    return end


def parse_cs(cs, r_st, q_st):
    """Yield (kind, ref_pos, ref_len, alt) from a minimap2 short cs string; positions on the reference."""
    r = r_st
    for op, val in re.findall(r'([:*+\-~])([0-9]+|[acgtn]+|[acgtn]{2})', cs):
        if op == ':':
            r += int(val)
        elif op == '*':
            yield ('S', r, 1, val[1].upper()); r += 1
        elif op == '+':
            yield ('I', r, 0, val.upper())
        elif op == '-':
            yield ('D', r, len(val), ''); r += len(val)


class RefAligner:
    def __init__(self, ref, preset='asm20'):
        self.ref = ref
        self.al = mp.Aligner(seq=ref, preset=preset)

    def align(self, q):
        """Greedy chain of + strand hits non-overlapping on reference and colinear on the query.
        Returns (covered intervals list, SNVs dict pos->alt, indel list [(pos, kind, len)], hits)."""
        hits = [h for h in self.al.map(q, cs=True) if h.strand == 1]
        hits.sort(key=lambda h: -h.mlen)
        chosen = []
        for h in hits:
            ok = True
            for c in chosen:
                if h.r_st < c.r_en and c.r_st < h.r_en:
                    ok = False; break
                if (h.r_st < c.r_st) != (h.q_st < c.q_st):
                    ok = False; break
            if ok:
                chosen.append(h)
        chosen.sort(key=lambda h: h.r_st)
        cov, snv, ind = [], {}, []
        for h in chosen:
            cov.append((h.r_st, h.r_en))
            for kind, p, ln, alt in parse_cs(h.cs, h.r_st, h.q_st):
                if kind == 'S':
                    snv[p] = alt
                else:
                    ind.append((p, kind, ln if kind == 'D' else len(alt)))
        return cov, snv, ind, chosen


def contiguous_extent(cov, start, max_gap=2000):
    """Furthest reference position reached from `start` through covered intervals separated by gaps <= max_gap."""
    reach = None
    for a, b in sorted(cov):
        if reach is None:
            if a <= start + max_gap and b > start:
                reach = b
        elif a <= reach + max_gap:
            reach = max(reach, b)
        else:
            break
    return reach or start
