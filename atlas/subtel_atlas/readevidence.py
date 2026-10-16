"""Read-level evidence for non-native terminal blocks (assembly-free; needs CHM13 v2.0 and the reads).

A claim says that the telomere-adjacent block of a chromosome end of recipient arm R is closer to the CHM13 copy of another end,
the donor D, than to CHM13 R, and that the sequence switches back to the proximal comparator P (normally CHM13 R itself, or the
next block's best CHM13 end) at a given position. A single long read supports the claim when it

  (a) starts in a terminal telomere (subtel_atlas.telomere.find_block on the read in tip frame),
  (b) carries D alleles at paralog-diagnostic units distal to the switch and P alleles proximal to it,
  (c) the maximum-likelihood change point along the read (with a profile-likelihood interval) separates the two, and
  (d) the read continues into an anchor that aligns genome-wide to CHM13 v2.0 with MAPQ 60 on the recipient chromosome,
  (e) with a colinear, single-strand structure (no foldback, internal telomere or adapter).

Diagnostic units are of two kinds, both defined from CHM13 only (independent of reads and assemblies):
  * SNV sites: single-base differences in a pairwise alignment of the D window to the P window (alignment blocks >= min_block bp
    at >= min_block_ident), excluding sites near other differences, homopolymer-adjacent sites, sites in low-complexity or
    locally repeated sequence and sites in the CHM13 telomere;
  * private bins: read bins (bin bp) that align to D (or P) sequence that has no counterpart in the other window at
    >= min_block_ident; the bin is called for the reference it matches better by >= bin_delta identity.

Thresholds live in a criteria JSON (evidence/criteria_v1.json; a copy ships in the package data). All computations are
deterministic for given inputs, criteria and minimap2/mappy versions.
"""
import os, re, json, math, gzip, subprocess, tempfile, hashlib
from importlib import resources
import numpy as np
from .common import rc, readfa, CHM13_ENDS
from .telomere import find_block, CAN

LN10 = math.log(10)
_CS = re.compile(r'(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})')
TELO_RUN = re.compile(r'(?:CCCTAA|TTAGGG)')
ADAPTERS = {   # public adapter/primer sequences (ONT ligation Y-adapter top strand, ONT rapid adapter, PacBio SMRTbell hairpin)
    'ont_lsk_y_top': 'AATGTACTTCGTTCAGTTACGTATTGCT',
    'ont_rap_y': 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA',
    'pacbio_smrtbell': 'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT',
}


# --------------------------------------------------------------------------------------------------------------------- criteria
def default_criteria_path():
    return str(resources.files('subtel_atlas').joinpath('data/criteria_v1.json'))


def load_criteria(path=None):
    with open(path or default_criteria_path()) as f:
        c = json.load(f)
    return c


def criteria_hash(c):
    return hashlib.sha256(json.dumps(c['thresholds'], sort_keys=True).encode()).hexdigest()[:16]


# --------------------------------------------------------------------------------------------------------------------- references
def tip(seq, end):
    """CHM13 end sequence in tip frame (position 0 = telomeric end); q ends of chm13_ends500k.fa finish at the telomere."""
    return seq if end.endswith('p') else rc(seq)


class Chm13Ends:
    """The 48 CHM13 v2.0 500 kb ends in tip frame, with the CHM13 telomere length of each end."""
    def __init__(self, fasta):
        self.seq = {n: tip(s.upper(), n) for n, s in readfa(fasta)}
        self.telo = {}
        for n, s in self.seq.items():
            b = find_block(s[:30000]); self.telo[n] = (b['junk_len'] + b['telo_len']) if b else 0

    def __getitem__(self, n): return self.seq[n]


def parse_cs(cs, r_st, q_st):
    """mappy/minimap2 short cs string -> (ref_pos, query_pos) arrays of aligned columns, arrays of mismatch columns
    (ref_pos, query_pos, ref_base, query_base) and the ref positions of every difference (mismatch or indel edge)."""
    rp, qp, mm, diff = [], [], [], []
    r, q = r_st, q_st
    for op in _CS.findall(cs):
        c = op[0]
        if c == ':':
            n = int(op[1:]); rp.append(np.arange(r, r + n)); qp.append(np.arange(q, q + n)); r += n; q += n
        elif c == '*':
            rp.append(np.array([r])); qp.append(np.array([q])); mm.append((r, q, op[1].upper(), op[2].upper())); diff.append(r)
            r += 1; q += 1
        elif c == '+':
            diff.append(r); q += len(op) - 1
        elif c == '-':
            diff.extend((r, r + len(op) - 2)); r += len(op) - 1
        else:
            raise ValueError('splice operators not expected')
    R = np.concatenate(rp) if rp else np.zeros(0, int); Q = np.concatenate(qp) if qp else np.zeros(0, int)
    return R, Q, mm, np.array(sorted(diff), int)


def _greedy_hits(hits, min_len, min_ident, max_ovl=0.1):
    """Forward-strand hits, best first by matching bases, kept when they overlap accepted hits by <= max_ovl in both coordinates."""
    out = []
    for h in sorted([h for h in hits if h.strand == 1 and (h.q_en - h.q_st) >= min_len and h.mlen / max(1, h.blen) >= min_ident],
                    key=lambda h: (-h.mlen, h.q_st, h.r_st)):
        ok = True
        for g in out:
            oq = min(h.q_en, g.q_en) - max(h.q_st, g.q_st); orr = min(h.r_en, g.r_en) - max(h.r_st, g.r_st)
            if oq > max_ovl * (h.q_en - h.q_st) or orr > max_ovl * (h.r_en - h.r_st): ok = False; break
        if ok: out.append(h)
    return out


def _hp_adjacent(s, i, k):
    """site i touches a homopolymer run of >= k bases (left or right neighbourhood, or the site extends a run)"""
    b = s[i]
    L = 0; j = i - 1
    while j >= 0 and s[j] == s[i - 1]: L += 1; j -= 1
    R = 0; j = i + 1
    while j < len(s) and s[j] == s[min(i + 1, len(s) - 1)]: R += 1; j += 1
    if L >= k or R >= k: return True
    # run including the site base
    run = 1; j = i - 1
    while j >= 0 and s[j] == b: run += 1; j -= 1
    j = i + 1
    while j < len(s) and s[j] == b: run += 1; j += 1
    return run >= k


def _low_complexity(s, i, half, min_tri):
    w = s[max(0, i - half): i + half + 1]
    return len({w[j:j + 3] for j in range(len(w) - 2)}) < min_tri


def _local_repeat(s, i, k, span):
    a = i - k // 2
    if a < 0 or a + k > len(s): return True
    kmer = s[a:a + k]; w = s[max(0, i - span): i + span]
    return w.count(kmer) > 1 or w.count(rc(kmer)) > 0


class PairMap:
    """Pairwise CHM13 D-window vs P-window map: homologous blocks, filtered SNV diagnostic sites and private masks.
    Coordinates are tip-frame positions on the CHM13 ends (distance from the telomeric end)."""
    def __init__(self, ends, D, P, d_win, p_win, th):
        self.D, self.P, self.d_win, self.p_win = D, P, tuple(d_win), tuple(p_win)
        Ds = ends[D][d_win[0]:d_win[1]]; Ps = ends[P][p_win[0]:p_win[1]]
        import mappy as mp
        al = mp.Aligner(seq=Ps, preset='asm20')
        hits = _greedy_hits(list(al.map(Ds, cs=True)), th['pair_min_block'], th['pair_min_block_ident'])
        self.blocks = []; sites = []
        self.th = th
        dcov = np.zeros(len(Ds), bool); pcov = np.zeros(len(Ps), bool)
        tD, tP = ends.telo[D], ends.telo[P]
        F = th['site_flank_clean']
        for h in sorted(hits, key=lambda h: h.q_st):
            Rr, Qq, mm, diff = parse_cs(h.cs, h.r_st, h.q_st)
            dcov[h.q_st:h.q_en] = True; pcov[h.r_st:h.r_en] = True
            self.blocks.append(dict(d_start=h.q_st + d_win[0], d_end=h.q_en + d_win[0], p_start=h.r_st + p_win[0],
                                    p_end=h.r_en + p_win[0], ident=round(h.mlen / h.blen, 5)))
            for (r, q, pb, db) in mm:
                k = np.searchsorted(diff, r)
                near = [x for x in diff[max(0, k - 2):k + 3] if x != r]
                if near and min(abs(x - r) for x in near) <= F: continue
                gd, gp = q + d_win[0], r + p_win[0]
                if gd < tD + th['site_min_from_chm13_telomere'] or gp < tP + th['site_min_from_chm13_telomere']: continue
                if pb not in 'ACGT' or db not in 'ACGT': continue
                bad = False
                for s, i in ((Ds, q), (Ps, r)):
                    if (_hp_adjacent(s, i, th['site_homopolymer'])
                            or _low_complexity(s, i, th['site_lowcomplex_half'], th['site_lowcomplex_min_trimers'])
                            or _local_repeat(s, i, th['site_repeat_k'], th['site_repeat_span'])):
                        bad = True; break
                if bad: continue
                sites.append((gd, gp, db, pb))
        self.sites = np.array([(a, b) for a, b, _, _ in sites], int).reshape(-1, 2)
        self.site_alleles = [(c, d) for _, _, c, d in sites]
        self.d_private = ~dcov; self.p_private = ~pcov      # window-relative masks
        self.n_snv_sites = len(sites)

    def summary(self):
        return dict(D=self.D, P=self.P, d_win=list(self.d_win), p_win=list(self.p_win), n_blocks=len(self.blocks),
                    homologous_bp=int(sum(b['d_end'] - b['d_start'] for b in self.blocks)), n_snv_sites=self.n_snv_sites,
                    d_private_bp=int(self.d_private.sum()), p_private_bp=int(self.p_private.sum()), blocks=self.blocks)


def shared_mask(ends, R, min_ident=0.95, min_len=1000):
    """Positions (tip frame) of CHM13 end R covered by an alignment of any other CHM13 end at >= min_ident over >= min_len bp."""
    import mappy as mp
    m = np.zeros(len(ends[R]), bool)
    al = mp.Aligner(seq=ends[R], preset='asm20')
    for o in CHM13_ENDS:
        if o == R or o not in ends.seq: continue
        for h in al.map(ends[o]):
            if h.blen >= min_len and h.mlen / h.blen >= min_ident: m[h.r_st:h.r_en] = True
    return m


# --------------------------------------------------------------------------------------------------------------------- claims
def claims_from_catalogue(ends_tsv, blocks_tsv, end_ids=None, th=None):
    """One claim per non-native (foreign or within_family) end: the outermost non-home core block.
    pattern DP  = the block is terminal (telomere -> D -> P);  PDP = the block is internal (telomere -> P -> D -> P).
    D = best CHM13 end covering most bases of the block (raw blocks); P = CHM13 end of the home arm when the next core block is
    home or absent, otherwise the next core block's dominant best end. d_lo/d_hi = span of the block on D (chm13_dist);
    p_switch = position on P where the next core block starts (assembly switch projected to CHM13 P), when defined."""
    import pandas as pd
    E = pd.read_csv(ends_tsv, sep='\t', low_memory=False); B = pd.read_csv(blocks_tsv, sep='\t', low_memory=False)
    E = E[E.nativity_class.isin(['foreign', 'within_family'])]
    if end_ids is not None: E = E[E.end_id.isin(set(end_ids))]
    out = []
    for e in E.itertuples():
        core = B[(B.end_id == e.end_id) & (B.level == 'core')].sort_values('block')
        raw = B[(B.end_id == e.end_id) & (B.level == 'raw')].sort_values('block')
        nh = core[~core.label.isin(['home', 'none'])]
        if not len(nh): continue
        k = nh.iloc[0]; idx = list(core.block).index(k.block)
        inside = raw[(raw.sub_start >= k.sub_start - 1) & (raw.sub_end <= k.sub_end + 1) & raw.best_end.notna()]
        if not len(inside): continue
        bp = (inside.sub_end - inside.sub_start).groupby(inside.best_end).sum().sort_values(ascending=False)
        D = bp.index[0]; dd = inside[inside.best_end == D]
        R = 'chr' + e.home_arm
        nxt = core.iloc[idx + 1] if idx + 1 < len(core) else None
        if nxt is None or nxt.label == 'home': P = R
        else:
            rn = raw[(raw.sub_start >= nxt.sub_start - 1) & (raw.sub_end <= nxt.sub_end + 1) & raw.best_end.notna()]
            P = (rn.sub_end - rn.sub_start).groupby(rn.best_end).sum().sort_values(ascending=False).index[0] if len(rn) else R
        p_sw = None
        if nxt is not None:
            first = raw[raw.sub_start >= nxt.sub_start - 1].head(1)
            if len(first) and first.best_end.iloc[0] == P: p_sw = int(first.chm13_dist_start.iloc[0])
        out.append(dict(claim_id=e.end_id, end_id=e.end_id, sample=e.sample, haplotype=e.haplotype, recipient_arm=e.home_arm, R=R,
                        D=D, P=P, pattern='DP' if idx == 0 else 'PDP', block_label=k.label, nativity_class=e.nativity_class,
                        core_structure=e.core_structure, switch_sub=int(k.sub_end) if nxt is not None else None,
                        switch_sub_distal=int(k.sub_start) if idx > 0 else None, window_sub_len=int(e.sub_len),
                        d_lo=int(dd.chm13_dist_start.min()), d_hi=int(dd.chm13_dist_end.max()), p_switch=p_sw,
                        contig=e.contig, side=e.side, contig_len=int(e.contig_len), s0=int(e.junk_len + e.telo_len)))
    return out


def claim_windows(c, th):
    """D and P windows (tip frame on CHM13) used for the diagnostic units of a claim."""
    pad = th['claim_window_pad']
    d_win = (max(0, c['d_lo'] - pad), c['d_hi'] + pad)
    reach = [x for x in (c.get('p_switch'), c.get('switch_sub'), c['d_hi']) if x is not None]
    p_hi = min(th['claim_p_window_max'], (max(reach) + th['claim_p_window_beyond']) if c.get('switch_sub') is not None
               else th['claim_p_window_max'])
    return d_win, (0, p_hi)


def shared_mask_cached(ends, R, th, cache_dir=None):
    """shared_mask() with an optional on-disk cache (packed bits, keyed by end, thresholds and the end's sequence hash)."""
    key = hashlib.sha256((ends[R] + f"{th['shared_mask_ident']}:{th['shared_mask_min_len']}").encode()).hexdigest()[:12]
    fn = os.path.join(cache_dir, f'{R}.{key}.shared.npy') if cache_dir else None
    if fn and os.path.exists(fn):
        return np.unpackbits(np.load(fn))[:len(ends[R])].astype(bool)
    m = shared_mask(ends, R, th['shared_mask_ident'], th['shared_mask_min_len'])
    if fn:
        os.makedirs(cache_dir, exist_ok=True); tmp = fn + '.tmp.npy'; np.save(tmp, np.packbits(m)); os.replace(tmp, fn)
    return m


# --------------------------------------------------------------------------------------------------------------------- reads
class Proj:
    """Projection of a tip-frame read onto one reference window (forward-strand mappy hits, greedy non-overlapping in the read)."""
    def __init__(self, aligner, read, win_len, th):
        n = len(read)
        self.x2q = np.full(win_len, -1, np.int64); self.q2x = np.full(n, -1, np.int64)
        self.match = np.zeros(n, bool); self.rev_bp = 0; self.hits = []
        hs = list(aligner.map(read, cs=True))
        self.rev_bp = int(sum(h.q_en - h.q_st for h in hs if h.strand == -1 and h.mapq >= th['foldback_min_mapq']))
        for h in _greedy_hits(hs, th['read_min_hit_bp'], 0.0):
            Rr, Qq, mm, _ = parse_cs(h.cs, h.r_st, h.q_st)
            self.x2q[Rr] = Qq; self.q2x[Qq] = Rr; self.match[Qq] = True
            if mm: self.match[np.array([q for _, q, _, _ in mm])] = False
            self.hits.append((h.q_st, h.q_en, h.r_st, h.r_en, round(h.mlen / h.blen, 4)))


def telomere_motifs(th, platform):
    """Canonical (CCCTAA)n rotations; for platforms listed in telomere_artefact_platforms also the published ONT R9 telomere
    basecalling artefacts in tip frame (Tan et al. 2022, doi:10.1186/s13059-022-02751-6)."""
    if th and platform in th.get('telomere_artefact_platforms', []):
        m = set(CAN)
        for u in th.get('telomere_artefact_motifs', []): m |= {u[i:] + u[:i] for i in range(len(u))}
        return m
    return CAN


def orient(seq, side=None, th=None, motifs=CAN):
    """Tip-frame orientation of a read: with side (p/L or q/R, for reads stored in reference orientation) the read is used as
    stored (p) or reverse-complemented (q); without side, the orientation with the longer terminal telomere block is chosen.
    Returns (tip_seq, orientation '+'/'-', telomere dict or None)."""
    s = seq.upper()
    cands = [('+', s), ('-', rc(s))]
    if side is not None:
        cands = [cands[0]] if side in ('p', 'L') else [cands[1]]
    best = None
    for o, x in cands:
        b = find_block(x[:70000], motifs=motifs, maxjunk=max(th['telomere_max_offset'] + 1, 1) if th else 8000)
        score = b['telo_len'] if b else -1
        if best is None or score > best[0]: best = (score, o, x, b)
    return best[2], best[1], best[3]


def _canon_mask(s):
    m = np.zeros(len(s), bool)
    for k in range(6):
        for mot in ('CCCTAA', 'TTAGGG'):
            u = mot[k:] + mot[:k]
            i = s.find(u)
            while i != -1:
                m[i:i + 6] = True; i = s.find(u, i + 1)
    return m


def internal_telomere(tipseq, start, th):
    """Largest block (bp) of >= 80% telomeric-hexamer coverage in 500 bp windows beyond `start`."""
    s = tipseq[start:]
    if len(s) < th['internal_telomere_min_bp']: return 0
    m = _canon_mask(s); W = th['internal_telomere_min_bp']
    cs = np.concatenate([[0], np.cumsum(m)]); fr = (cs[W:] - cs[:-W]) / W
    hit = np.where(fr >= 0.8)[0]
    if not len(hit): return 0
    # merge windows into blocks; report the largest block span
    runs = np.split(hit, np.where(np.diff(hit) > 1)[0] + 1)
    return int(max(r[-1] - r[0] + W for r in runs))


def adapter_hits(tipseq, th):
    try:
        import edlib
    except ImportError:
        return None
    out = []; n = len(tipseq); ex = th['adapter_edge_exclude']
    for name, a in ADAPTERS.items():
        for aa in (a, rc(a)):
            r = edlib.align(aa, tipseq, mode='HW', task='locations', k=int(th['adapter_max_edit_frac'] * len(aa)))
            for (st, en) in (r.get('locations') or []):
                if st is not None and st > ex and en < n - ex: out.append((name, int(st))); break
    return out


def _log10_binom_tail(k, n, p):
    """log10 P(X >= k), X ~ Binomial(n, p)"""
    if k <= 0: return 0.0
    if k > n: return -math.inf
    lp, lq = math.log(p), math.log1p(-p)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq for i in range(k, n + 1)]
    m = max(terms)
    return (m + math.log(sum(math.exp(t - m) for t in terms))) / LN10


def changepoint(y, pattern, e, drop):
    """Maximum-likelihood segmentation of unit labels y (1 = D, 0 = P; read order) for pattern 'DP' or 'PDP'.
    Returns dict with cuts (index of the first unit of each new segment), profile-likelihood index intervals, log-likelihoods."""
    y = np.asarray(y, int); n = len(y)
    l1, l0 = math.log(1 - e), math.log(e)
    cD = np.concatenate([[0], np.cumsum(y)]); idx = np.arange(n + 1)
    def fD(i, j): return (cD[j] - cD[i]) * l1 + ((j - i) - (cD[j] - cD[i])) * l0      # segment expected D
    def fP(i, j): return ((j - i) - (cD[j] - cD[i])) * l1 + (cD[j] - cD[i]) * l0
    allP, allD = fP(0, n), fD(0, n)
    if pattern == 'DP':
        if n < 2: return None
        c = idx[1:n]; LL = fD(0, c) + fP(c, n); k = int(np.argmax(LL)); best = float(LL[k])
        ok = c[LL >= best - drop]
        return dict(cuts=[int(c[k])], ci=[(int(ok.min()), int(ok.max()))], ll=best, ll_allP=float(allP), ll_allD=float(allD),
                    segs=[(0, int(c[k]), 1), (int(c[k]), n, 0)])
    if pattern == 'PDP':
        if n < 3: return None
        best = -math.inf; ba = bb = None; prof_a = np.full(n + 1, -math.inf); prof_b = np.full(n + 1, -math.inf)
        for a in range(1, n - 1):
            b = idx[a + 1:n]; LL = fP(0, a) + fD(a, b) + fP(b, n)
            j = int(np.argmax(LL)); prof_a[a] = LL[j]
            prof_b[a + 1:n] = np.maximum(prof_b[a + 1:n], LL)
            if LL[j] > best: best = float(LL[j]); ba, bb = a, int(b[j])
        oka = np.where(prof_a >= best - drop)[0]; okb = np.where(prof_b >= best - drop)[0]
        return dict(cuts=[ba, bb], ci=[(int(oka.min()), int(oka.max())), (int(okb.min()), int(okb.max()))], ll=best,
                    ll_allP=float(allP), ll_allD=float(allD), segs=[(0, ba, 0), (ba, bb, 1), (bb, n, 0)])
    raise ValueError(pattern)


class ClaimContext:
    """Everything reference-side that a claim needs: unit sets, per-window aligners and the recipient's shared mask."""
    def __init__(self, claim, ends, crit, preset='map-ont', mask_cache=None):
        import mappy as mp
        th = crit['thresholds']; self.th = th; self.crit = crit; self.c = claim; self.preset = preset
        self.d_win, self.p_win = claim_windows(claim, th)
        D, P, R = claim['D'], claim['P'], claim['R']
        self.pm = PairMap(ends, D, P, self.d_win, self.p_win, th)
        self.pmR = PairMap(ends, D, R, self.d_win, self.p_win, th) if R != P else None
        self.al = {D: mp.Aligner(seq=ends[D][self.d_win[0]:self.d_win[1]], preset=preset),
                   P: mp.Aligner(seq=ends[P][self.p_win[0]:self.p_win[1]], preset=preset)}
        if R != P: self.al[R] = mp.Aligner(seq=ends[R][self.p_win[0]:self.p_win[1]], preset=preset)
        self.mask = shared_mask_cached(ends, R, th, mask_cache)

    def summary(self):
        s = dict(claim=self.c, d_win=list(self.d_win), p_win=list(self.p_win), units_DP=self.pm.summary())
        if self.pmR is not None: s['units_DR'] = self.pmR.summary()
        return s


def _units(pm, read, prD, prX, th, start):
    """Diagnostic units of a read for unit set pm (D vs X): list of (read_pos, label, kind, d_pos, x_pos);
    label 1 = D allele, 0 = X allele, 2 = other (SNV only)."""
    out = []; d0, x0 = pm.d_win[0], pm.p_win[0]
    if len(pm.sites):
        qd = prD.x2q[pm.sites[:, 0] - d0]; qx = prX.x2q[pm.sites[:, 1] - x0]
        for i in range(len(pm.sites)):
            a, b = int(qd[i]), int(qx[i])
            if a < 0 and b < 0: continue
            if a >= 0 and b >= 0 and a != b and read[a] != read[b]: q, base = (a, None)
            else: q = a if a >= 0 else b; base = read[q]
            if q < start: continue
            dal, xal = pm.site_alleles[i]
            lab = 1 if base == dal else 0 if base == xal else 2
            out.append((q, lab, 'snv', int(pm.sites[i, 0]), int(pm.sites[i, 1])))
    B = th['bin_bp']; n = len(read)
    for s in range(start, n - B + 1, B):
        sl = slice(s, s + B)
        dx = prD.q2x[sl]; xx = prX.q2x[sl]
        dpriv = np.count_nonzero((dx >= 0) & pm.d_private[np.clip(dx, 0, None)]) / B
        xpriv = np.count_nonzero((xx >= 0) & pm.p_private[np.clip(xx, 0, None)]) / B
        iD = prD.match[sl].sum() / B; iX = prX.match[sl].sum() / B
        if dpriv >= th['bin_min_private_frac'] and iD >= th['bin_min_ident'] and iD - iX >= th['bin_delta']:
            dpos = int(np.median(dx[dx >= 0])) + d0
            out.append((s + B // 2, 1, 'bin', dpos, -1))
        elif xpriv >= th['bin_min_private_frac'] and iX >= th['bin_min_ident'] and iX - iD >= th['bin_delta']:
            xpos = int(np.median(xx[xx >= 0])) + x0
            out.append((s + B // 2, 0, 'bin', -1, xpos))
    out.sort(key=lambda u: (u[0], u[2]))
    return out


def evaluate_read(ctx, name, seq, side=None, meta=None):
    """Stage 1 (reference-window evidence) for one read. Returns (row dict, per-unit list, tip-frame read)."""
    th = ctx.th; c = ctx.c; D, P, R = c['D'], c['P'], c['R']
    t, o, tb = orient(seq, side, th, telomere_motifs(th, (meta or {}).get('platform')))
    row = dict(claim_id=c['claim_id'], read=name, read_len=len(t), orientation=o)
    if meta: row.update(meta)
    terminal = bool(tb and tb['junk_len'] <= th['telomere_max_offset'] and tb['telo_len'] >= th['telomere_min_bp'])
    row.update(telo_offset=tb['junk_len'] if tb else None, telo_len=tb['telo_len'] if tb else 0, terminal=terminal)
    s0 = (tb['junk_len'] + tb['telo_len']) if tb else 0; row['s0_read'] = s0
    pr = {k: Proj(a, t, (ctx.d_win[1] - ctx.d_win[0]) if k == D else (ctx.p_win[1] - ctx.p_win[0]), th) for k, a in ctx.al.items()}
    U = _units(ctx.pm, t, pr[D], pr[P], th, s0)
    snv = [u for u in U if u[2] == 'snv']; bins = [u for u in U if u[2] == 'bin']
    row.update(n_sites=int(ctx.pm.n_snv_sites), snv_D=sum(u[1] == 1 for u in snv), snv_P=sum(u[1] == 0 for u in snv),
               snv_other=sum(u[1] == 2 for u in snv), bin_D=sum(u[1] == 1 for u in bins), bin_P=sum(u[1] == 0 for u in bins))
    row['snv_uncovered'] = row['n_sites'] - len(snv)
    inf = [u for u in U if u[1] in (0, 1)]; y = [u[1] for u in inf]; pos = [u[0] for u in inf]
    cp = changepoint(y, c['pattern'], th['cp_error'], th['cp_ci_drop_ln']) if len(inf) >= 2 else None
    row['n_units'] = len(inf)
    sw_ok = False; ci_hi_read = None
    if cp:
        e = th['cp_error']
        row['log10_lr_vs_allP'] = round((cp['ll'] - cp['ll_allP']) / LN10, 3)
        row['log10_lr_vs_allD'] = round((cp['ll'] - cp['ll_allD']) / LN10, 3)
        segok = True; sg = []
        for (i, j, lab) in cp['segs']:
            nn_ = j - i; k = sum(1 for v in y[i:j] if v == lab); pur = k / nn_ if nn_ else 0.0
            sg.append(dict(n=nn_, purity=round(pur, 4), log10p=round(_log10_binom_tail(k, nn_, e), 2), lab='D' if lab else 'P'))
            if nn_ < th['min_units_segment'] or pur < th['min_segment_purity']: segok = False
        row['segments'] = ';'.join(f"{s['lab']}:{s['n']}:{s['purity']}" for s in sg)
        dseg = [s for s in sg if s['lab'] == 'D'][0]
        row.update(seg_D_units=dseg['n'], seg_D_purity=dseg['purity'], log10p_binom_D=dseg['log10p'])
        pseg = sg[-1]; row.update(seg_P_units=pseg['n'], seg_P_purity=pseg['purity'], log10p_binom_P=pseg['log10p'])
        # proximal (D -> P) switch = last cut
        cut = cp['cuts'][-1]; lo, hi = cp['ci'][-1]
        r_lo = pos[max(lo - 1, 0)]; r_hi = pos[min(hi, len(pos) - 1)]
        row.update(switch_read=int((pos[cut - 1] + pos[cut]) // 2), switch_ci_read_lo=int(r_lo), switch_ci_read_hi=int(r_hi),
                   switch_ci_sub_lo=int(r_lo - s0), switch_ci_sub_hi=int(r_hi - s0))
        pp = inf[min(hi, len(inf) - 1)]
        p_at_hi = pp[4] if pp[4] >= 0 else (int(pr[P].q2x[pp[0]]) + ctx.p_win[0] if pr[P].q2x[pp[0]] >= 0 else None)
        row['switch_ci_P_hi'] = p_at_hi
        row['switch_ci_P_lo'] = (p_at_hi - (r_hi - r_lo)) if p_at_hi is not None else None
        if c['pattern'] == 'PDP':
            a = cp['cuts'][0]; alo, ahi = cp['ci'][0]
            row.update(switch1_ci_sub_lo=int(pos[max(alo - 1, 0)] - s0), switch1_ci_sub_hi=int(pos[min(ahi, len(pos) - 1)] - s0))
        lr_ok = row['log10_lr_vs_allP'] >= th['min_log10_lr'] and row['log10_lr_vs_allD'] >= th['min_log10_lr']
        sw_ok = segok and lr_ok
        # distal segment vs R (when P is not R)
        if ctx.pmR is not None:
            dlo = pos[cp['segs'][[s[2] for s in cp['segs']].index(1)][0]]; dhi = r_lo
            UR = [u for u in _units(ctx.pmR, t, pr[D], pr[R], th, s0) if u[1] in (0, 1) and dlo <= u[0] <= dhi]
            nR = len(UR); kR = sum(u[1] for u in UR)
            row.update(dR_units=nR, dR_purity=round(kR / nR, 4) if nR else 0.0)
            if nR < th['min_units_segment'] or kR / max(nR, 1) < th['min_segment_purity']: sw_ok = False
        ci_hi_read = r_hi
    row['switch_pass'] = bool(sw_ok)
    # claimed donor interval (subtelomere coordinates of the read) for the contradiction test
    tol = th['concordance_tol_bp']
    a = (c.get('switch_sub_distal') or 0); b = (c['switch_sub'] if c.get('switch_sub') is not None else c['window_sub_len']) - tol
    di = [u for u in inf if a <= u[0] - s0 <= b]
    row['donor_interval_units'] = len(di)
    row['donor_interval_P_purity'] = round(sum(1 for u in di if u[1] == 0) / len(di), 4) if di else None
    # concordance with the assembly switch
    conc = None
    if sw_ok:
        conc = False
        if c.get('switch_sub') is not None:
            if row['switch_ci_sub_lo'] - tol <= c['switch_sub'] <= row['switch_ci_sub_hi'] + tol: conc = True
        if c.get('p_switch') is not None and row.get('switch_ci_P_hi') is not None:
            if row['switch_ci_P_lo'] - tol <= c['p_switch'] <= row['switch_ci_P_hi'] + tol: conc = True
    row['switch_concordant'] = conc
    # anchor segment start
    row['anchor_from'] = int(ci_hi_read) if sw_ok else int(s0 + (c['switch_sub'] if c.get('switch_sub') is not None else 0))
    row['internal_telomere_bp'] = internal_telomere(t, s0 + th['internal_telomere_gap'], th)
    ad = adapter_hits(t, th); row['adapter_hits'] = None if ad is None else ','.join(f'{a}@{p}' for a, p in ad)
    row['rev_hits_bp_D_P'] = pr[D].rev_bp + pr[P].rev_bp
    row['proj_hits'] = json.dumps({k: v.hits[:6] for k, v in pr.items()})
    row['_proj_iv'] = sorted({(h[0], h[1]) for v in pr.values() for h in v.hits})
    units = [dict(read_pos=int(u[0]), allele='D' if u[1] == 1 else 'P' if u[1] == 0 else 'other', kind=u[2], d_pos=u[3], p_pos=u[4])
             for u in U]
    return row, units, t


# --------------------------------------------------------------------------------------------------------------------- genome
def find_minimap2(exe=None):
    from shutil import which
    return exe or os.environ.get('SUBTEL_ATLAS_MINIMAP2') or which('minimap2')


def map_genome(records, genome, threads=1, exe=None, sam_out=None, preset='map-ont'):
    """records: list of (name, seq). minimap2 -a -Y (soft-clipped supplementary) --secondary=no against the CHM13 v2.0 genome
    (FASTA or .mmi). Returns {name: [alignment dict]} with query coordinates on the given sequence; optionally keeps the SAM."""
    import pysam
    mm = find_minimap2(exe)
    if not mm: raise SystemExit('minimap2 executable needed for the genome-wide anchor (conda install -c bioconda minimap2)')
    td = tempfile.mkdtemp(prefix='readev_'); fa = os.path.join(td, 'q.fa'); sam = sam_out or os.path.join(td, 'q.sam')
    with open(fa, 'w') as f:
        for n, s in records: f.write(f'>{n}\n{s}\n')
    with open(sam, 'w') as so:
        subprocess.run([mm, '-a', '-Y', '--secondary=no', '-x', preset, '-t', str(threads), genome, fa], stdout=so,
                       stderr=subprocess.DEVNULL, check=True)
    out = {}
    with pysam.AlignmentFile(sam, 'r') as f:
        L = dict(zip(f.references, f.lengths))
        for r in f:
            if r.is_unmapped: continue
            ql = r.infer_read_length(); qs, qe = r.query_alignment_start, r.query_alignment_end
            if r.is_reverse: qs, qe = ql - qe, ql - qs
            out.setdefault(r.query_name, []).append(dict(chrom=r.reference_name, chrom_len=L[r.reference_name], strand='-' if r.is_reverse else '+',
                                                         ref_start=r.reference_start, ref_end=r.reference_end, q_start=qs, q_end=qe,
                                                         mapq=r.mapping_quality, supp=r.is_supplementary, aln_len=qe - qs))
    os.remove(fa)
    if not sam_out: os.remove(sam)
    os.rmdir(td) if not os.listdir(td) else None
    return out


def finalize(ctx, row, whole_alns, anchor_alns):
    """Stage 2: anchor, integrity and verdict from genome-wide alignments of the whole read and of its proximal segment."""
    th = ctx.th; R = ctx.c['R']; chrom = R[:-1]; arm = R[-1]
    want = '+' if arm == 'p' else '-'
    best = None
    off = row['anchor_from']
    for a in anchor_alns or []:
        if a['chrom'] != chrom or a['strand'] != want or a['mapq'] < th['anchor_min_mapq'] or a['aln_len'] < th['anchor_min_bp']: continue
        dist = a['ref_start'] if arm == 'p' else a['chrom_len'] - a['ref_end']
        if dist > th['anchor_max_from_telomere']: continue
        if best is None or a['aln_len'] > best['aln_len']: best = a
    row['anchor_pass'] = best is not None
    if best:
        L = best['chrom_len']
        t0, t1 = (best['ref_start'], best['ref_end']) if arm == 'p' else (L - best['ref_end'], L - best['ref_start'])
        m = ctx.mask; inside = max(0, min(t1, len(m)) - t0) if t0 < len(m) else 0
        uniq = (inside - int(m[t0:min(t1, len(m))].sum()) if inside > 0 else 0) + max(0, t1 - max(t0, len(m)))
        row.update(anchor_chrom=best['chrom'], anchor_mapq=best['mapq'], anchor_len=best['aln_len'], anchor_tip_start=int(t0),
                   anchor_tip_end=int(t1), anchor_read_start=int(off + best['q_start']), anchor_unique_bp=int(uniq),
                   anchor_tier='unique' if uniq >= th['anchor_unique_tier_bp'] else 'mapq60_shared')
    else:
        lng = max([a for a in anchor_alns or []], key=lambda a: a['aln_len'], default=None)
        row.update(anchor_chrom=lng['chrom'] if lng else None, anchor_mapq=lng['mapq'] if lng else None,
                   anchor_len=lng['aln_len'] if lng else 0, anchor_unique_bp=0, anchor_tier=None)
    # integrity
    W = whole_alns or []
    fold = any(a['chrom'] == chrom and a['strand'] != want and a['mapq'] >= th['foldback_min_mapq'] and a['aln_len'] >= th['foldback_min_bp']
               for a in W)
    s0 = row['s0_read']; n = row['read_len']; cov = np.zeros(max(0, n - s0), bool)
    for a in W:
        cov[max(0, a['q_start'] - s0):max(0, a['q_end'] - s0)] = True
    for (qs, qe) in row.get('_proj_iv', []):        # reference-window alignments (donor blocks that minimap2 drops as secondary)
        cov[max(0, qs - s0):max(0, qe - s0)] = True
    frac = float(cov.mean()) if len(cov) else 0.0
    row.update(foldback=fold, subtel_aligned_frac=round(frac, 4),
               genome_alignments=';'.join(f"{a['chrom']}:{a['ref_start']}-{a['ref_end']}{a['strand']}q{a['q_start']}-{a['q_end']}"
                                          f"Q{a['mapq']}" for a in sorted(W, key=lambda a: a['q_start'])))
    it_ok = row['internal_telomere_bp'] < th['internal_telomere_min_bp']
    ad_ok = not row.get('adapter_hits')
    row['integrity_pass'] = bool(it_ok and ad_ok and not fold and frac >= th['min_subtel_aligned_frac'])
    # verdict
    base = row['terminal'] and row['anchor_pass'] and row['integrity_pass']
    contra = (base and row['donor_interval_units'] >= th['min_units_segment']
              and (row['donor_interval_P_purity'] or 0) >= th['min_segment_purity'] and not row['switch_pass'])
    if row['switch_pass'] and base: v = 'conclusive'
    elif contra: v = 'contradicts'
    elif row['switch_pass']: v = 'consistent'
    else: v = 'non-informative'
    row['verdict'] = v
    return row


def assembly_projection(ctx, window_seq, name='assembly_window'):
    """The claim's unit model applied to the person's assembly end window (tip frame, error-free sequence): the assembly switch in
    the same diagnostic-unit frame as the reads (secondary, assembly-concordance line)."""
    r, _, _ = evaluate_read(ctx, name, window_seq, side='L', meta=dict(platform='assembly'))
    return {k: r.get(k) for k in ('switch_pass', 'switch_ci_sub_lo', 'switch_ci_sub_hi', 'switch_ci_P_lo', 'switch_ci_P_hi', 'segments',
                                  'log10_lr_vs_allP', 'log10_lr_vs_allD', 'seg_D_purity')}


def _ovl(a0, a1, b0, b1, tol):
    return None not in (a0, a1, b0, b1) and not (isinstance(a0, float) and math.isnan(a0)) and a0 - tol <= b1 and b0 <= a1 + tol


def concordance(row, asm, th):
    """switch_concordant (final) and its reference: the assembly projection when it shows the switch, else the catalogue switch."""
    row['switch_concordant_catalogue'] = row.get('switch_concordant')
    if not row.get('switch_pass'):
        row['concordance_ref'] = None; return row
    tol = th['concordance_tol_bp']
    if asm and asm.get('switch_pass'):
        ok = _ovl(row['switch_ci_sub_lo'], row['switch_ci_sub_hi'], asm['switch_ci_sub_lo'], asm['switch_ci_sub_hi'], tol) or \
            _ovl(row.get('switch_ci_P_lo'), row.get('switch_ci_P_hi'), asm.get('switch_ci_P_lo'), asm.get('switch_ci_P_hi'), tol)
        row['switch_concordant'] = bool(ok); row['concordance_ref'] = 'assembly_projection'
    else:
        row['concordance_ref'] = 'catalogue'
    return row


def assembly_fields(rec, claim):
    """Own-assembly alignment line for a pysam record (reads aligned to the person's assembly, e.g. HG002 Q100)."""
    ql = rec.infer_read_length() or 1
    qcov = (rec.query_alignment_end - rec.query_alignment_start) / ql
    has_sa = rec.has_tag('SA')
    onend = rec.reference_name == claim.get('contig')
    return dict(input_contig=rec.reference_name, input_pos=rec.reference_start, input_mapq=rec.mapping_quality, input_has_sa=has_sa,
                input_qcov=round(qcov, 4), assembly_concordant=bool(onend and not has_sa and qcov >= 0.95))


def claim_verdict(rows, crit, same_hap_key=None):
    """Claim-level rule (CRITERIA.md section 4). rows: finalized read rows of one claim (one haplotype-assignable set or unphased)."""
    th = crit['thresholds']
    conc = [r for r in rows if r['verdict'] == 'conclusive' and r.get('switch_concordant')]
    disc = [r for r in rows if r['verdict'] == 'conclusive' and not r.get('switch_concordant')]
    cons = [r for r in rows if r['verdict'] == 'consistent']
    contra = [r for r in rows if r['verdict'] == 'contradicts']
    ds = sorted({r.get('dataset') for r in conc if r.get('dataset')})
    phased = same_hap_key is not None
    if len(conc) >= th['claim_min_conclusive'] and not (phased and contra): v = 'confirmed'
    elif len(conc) >= 1 and cons and not (phased and contra): v = 'supported'
    elif phased and len(contra) >= 2 and not conc: v = 'contradicted'
    else: v = 'not_confirmed'
    return dict(verdict=v, n_conclusive_concordant=len(conc), n_conclusive_discordant=len(disc), n_consistent=len(cons),
                n_contradicts=len(contra), n_reads=len(rows), datasets_conclusive=','.join(ds), multi_dataset=len(ds) >= 2)


def iter_reads(path, region=None, reference=None, index=None, min_len=0):
    """Yield (name, seq_as_stored, pysam record or None) from BAM/CRAM (optionally a region; local path or URL) or FASTA/FASTQ.
    BAM/CRAM: primary alignments only, each read once; seq is in reference orientation."""
    if path.endswith(('.bam', '.cram')) or '.bam?' in path or '.cram?' in path:
        import pysam
        kw = dict(reference_filename=reference) if reference else {}
        if index: kw['index_filename'] = index
        seen = set()
        with pysam.AlignmentFile(path, 'rc' if path.endswith('.cram') else 'rb', **kw) as f:
            it = f.fetch(region=region) if region else f.fetch(until_eof=True)
            for r in it:
                if r.is_secondary or r.is_supplementary or not r.query_sequence or r.query_name in seen: continue
                if len(r.query_sequence) < min_len: continue
                seen.add(r.query_name); yield r.query_name, r.query_sequence, r
    else:
        import mappy
        for name, seq, _q in mappy.fastx_read(path):
            if len(seq) >= min_len: yield name, seq, None


def stage1(ctx, reads, side=None, dataset=None, platform=None, tag=''):
    """Reference-window evidence for reads (iterable of (name, seq, pysam record or None)) -> (rows, units, genome queries)."""
    rows, units, queries = [], {}, []
    for name, seq, rec in reads:
        meta = dict(dataset=dataset, platform=platform)
        if rec is not None: meta.update(assembly_fields(rec, ctx.c))
        r, u, t = evaluate_read(ctx, name, seq, side=side, meta=meta)
        k = f"{tag}{len(rows)}"; r['_key'] = k
        rows.append(r); units[name] = u
        queries.append((f"{k}|w", t))
        if len(t) - r['anchor_from'] >= ctx.th['anchor_min_bp']: queries.append((f"{k}|a", t[r['anchor_from']:]))
    return rows, units, queries


def stage2(ctx, rows, G):
    for r in rows:
        finalize(ctx, r, G.get(f"{r['_key']}|w"), G.get(f"{r['_key']}|a"))
    return rows


def evidence_doc(ctx, rows, units):
    crit = ctx.crit
    return dict(criteria_version=crit['criteria_version'], criteria_status=crit['status'], criteria_sha=criteria_hash(crit),
                context=ctx.summary(), sites=[dict(d_pos=int(a), p_pos=int(b), D=x, P=y) for (a, b), (x, y) in
                                              zip(ctx.pm.sites.tolist(), ctx.pm.site_alleles)],
                reads={r['read']: dict(verdict=r['verdict'], dataset=r.get('dataset'), units=units[r['read']]) for r in rows})


def run_claim(claim, reads, ends, crit, genome, preset='map-ont', threads=1, side=None, dataset=None, platform=None,
              mask_cache=None, minimap2=None, sam_out=None, assembly_window=None, log=lambda *a: None):
    """Evaluate reads (iterable of (name, seq, record)) for one claim. Returns (rows, json_doc). assembly_window: the person's
    tip-frame end window (sequence) for the assembly-projection concordance (criteria change C2)."""
    ctx = ClaimContext(claim, ends, crit, preset=preset, mask_cache=mask_cache)
    asm = assembly_projection(ctx, assembly_window.upper()) if assembly_window else None
    rows, units, Q = stage1(ctx, reads, side, dataset, platform)
    log(f"  {claim['claim_id']}: {len(rows)} reads evaluated (stage 1)")
    G = map_genome(Q, genome, threads, minimap2, sam_out=sam_out) if Q else {}
    stage2(ctx, rows, G)
    for r in rows: concordance(r, asm, ctx.th)
    doc = evidence_doc(ctx, rows, units); doc['assembly_projection'] = asm
    return rows, doc


READ_COLS = ['claim_id', 'read', 'dataset', 'platform', 'read_len', 'orientation', 'telo_offset', 'telo_len', 'terminal', 's0_read',
             'n_sites', 'snv_D', 'snv_P', 'snv_other', 'snv_uncovered', 'bin_D', 'bin_P', 'n_units', 'segments', 'seg_D_units',
             'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'log10_lr_vs_allP', 'log10_lr_vs_allD', 'log10p_binom_D', 'log10p_binom_P',
             'dR_units', 'dR_purity', 'switch_pass', 'switch_read', 'switch_ci_read_lo', 'switch_ci_read_hi', 'switch_ci_sub_lo',
             'switch_ci_sub_hi', 'switch1_ci_sub_lo', 'switch1_ci_sub_hi', 'switch_ci_P_lo', 'switch_ci_P_hi', 'switch_concordant',
             'switch_concordant_catalogue', 'concordance_ref',
             'donor_interval_units', 'donor_interval_P_purity', 'anchor_from', 'anchor_pass', 'anchor_chrom', 'anchor_mapq', 'anchor_len',
             'anchor_tip_start', 'anchor_tip_end', 'anchor_unique_bp', 'anchor_tier', 'internal_telomere_bp', 'adapter_hits', 'foldback',
             'subtel_aligned_frac', 'integrity_pass', 'input_contig', 'input_pos', 'input_mapq', 'input_has_sa', 'input_qcov',
             'assembly_concordant', 'verdict', 'genome_alignments']


def write_tsv(rows, fh, cols=READ_COLS):
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        fh.write('\t'.join('' if r.get(k) is None else str(r.get(k)) for k in cols) + '\n')
