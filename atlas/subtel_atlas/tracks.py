"""Track-based mosaic labelling (catalogue v0.1 design; shared by pipeline/catalogue_04t/06t and the query tool).

Why: 5 kb bins tiled from the telomere boundary made block labels depend on the tiling phase, and the single best-matching CHM13
end flips among near-identical paralogous ends. Here labels come from per-base identity tracks and are expressed over paralog
families.

1. Alignment: the subtelomere (window[s0:]) is aligned once to the 48 CHM13 v2.0 500 kb ends, keeping all chains
   (minimap2 MM2_TRACK_ARGS: asm20, base-level alignment, cs tags, up to 200 secondary chains with score >= 1% of the primary;
   on 40 test ends this gave the same family-level core structures as -P, which retains all chains, at 1/10 of the run time).
2. Tracks: for every CHM13 end E, M_E[q] = 1 if query base q is aligned as a match in any alignment to E, else 0 (mismatches,
   inserted bases and the base before a deletion are 0). The position of the matched base on E is kept for CHM13 coordinates.
3. Identity: id_E(x) = mean of M_E over [x - WIN/2, x + WIN/2) (clipped to the subtelomere), evaluated on a GRID-bp grid anchored
   at the contig tip (tip-frame positions k * GRID >= s0), so that the grid does not depend on the boundary call.
4. Labels per grid point: 'none' if the best identity < MIN_ID; 'home' if a home CHM13 end (the arm or its pooled group) is within
   MARGIN of the best identity; otherwise the paralog family of the best end (families: data/catalogue/paralog_families.tsv;
   an end in no multi-member family is its own family, named by its arm).
5. Blocks: runs of equal labels; runs shorter than MIN_RAW are merged iteratively (shortest first) into the neighbour whose label
   track has the higher mean match over the run; each remaining boundary is refined at base level by the change point minimising
   errors(outer label before x) + errors(inner label after x) on the label tracks (max over member ends) within REFINE bp.
6. Fixed windows S5/S10/S20 (and legacy t20 = tip[s0:20 kb]): mean match per end over the window; best end; own (best end is
   home), tie (a home end within MARGIN), foreign (best end's family), none (best < MIN_ID).
"""
import re, collections, numpy as np
from .common import arm_of, group_ends, ARM_GROUPS

MM2_TRACK_ARGS = ['-x', 'asm20', '-c', '--cs=short', '--secondary=yes', '-N', '200', '-p', '0.01']
GRID, WIN, MARGIN, MIN_ID, MIN_RAW, REFINE, EL = 100, 5000, 0.005, 0.80, 2000, 5000, 500_000
FIXED = (('S5', 5000), ('S10', 10000), ('S20', 20000))
FIXED_FROM_LABELS = True
_CS = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[+\-][a-z]+)')


def parse_paf(line):
    """-> (qname, qlen, qs, qe, strand, tname, ts, te, cs)"""
    f = line.rstrip('\n').split('\t')
    cs = next((t[5:] for t in f[12:] if t.startswith('cs:Z:')), '')
    return f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), cs


def _tokens(cs):
    """cs (short) -> arrays: op code (0 match, 1 mismatch, 2 insertion, 3 deletion), length"""
    toks = _CS.findall(cs)
    op = np.fromiter((':*+-'.index(t[0]) for t in toks), np.int8, len(toks))
    ln = np.fromiter((int(t[1:]) if t[0] == ':' else (1 if t[0] == '*' else len(t) - 1) for t in toks), np.int64, len(toks))
    return op, ln


class Track:
    """per-base match track of one CHM13 end for one query, plus the alignments for target-coordinate look-up"""
    __slots__ = ('M', 'alns')
    def __init__(self, n): self.M = np.zeros(n, np.int8); self.alns = []

    def tpos(self, q):
        """CHM13 coordinate of query base q in an alignment where q is a match, else -1"""
        for qs, qe, st, ts, op, ln, qa, ta in self.alns:
            if not (qs <= q < qe): continue
            r = (q - qs) if st == '+' else (qe - 1 - q)          # offset along the alignment in alignment (target) direction
            k = int(np.searchsorted(qa, r, side='right')) - 1
            if k >= 0 and op[k] == 0 and r < qa[k] + ln[k]: return int(ta[k] + (r - qa[k]))
        return -1


def tracks(n, alns):
    """alns: list of (qs, qe, strand, tname, ts, te, cs) for one query of length n. Returns {tname: Track}."""
    out = {}
    for qs, qe, st, tn, ts, te, cs in alns:
        tr = out.setdefault(tn, Track(n))
        op, ln = _tokens(cs)
        qadv = np.where(op == 3, 0, ln); tadv = np.where(op == 2, 0, ln)
        qa = np.concatenate([[0], np.cumsum(qadv)[:-1]]); ta = ts + np.concatenate([[0], np.cumsum(tadv)[:-1]])
        m = np.zeros(qe - qs + 1, np.int64); mk = op == 0
        np.add.at(m, qa[mk], 1); np.add.at(m, qa[mk] + ln[mk], -1)
        mm = (np.cumsum(m)[:-1] > 0).astype(np.int8)                  # in alignment direction
        dl = qa[op == 3] - 1; dl = dl[(dl >= 0) & (dl < len(mm))]; mm[dl] = 0   # base before a deletion
        if st == '-': mm = mm[::-1]
        np.maximum(tr.M[qs:qe], mm, out=tr.M[qs:qe]); tr.alns.append((qs, qe, st, ts, op, ln, qa, ta))
    return out


def grid_points(s0, L):
    """tip-frame grid positions (multiples of GRID) inside the subtelomere [s0, L)"""
    first = -(-s0 // GRID) * GRID
    return np.arange(first, L, GRID)


def identity_grid(tr, s0, L, names):
    """id matrix (len(names) x grid) of sliding WIN identity; grid in tip frame."""
    n = L - s0; g = grid_points(s0, L) - s0; h = WIN // 2
    a = np.clip(g - h, 0, n); b = np.clip(g + h, 0, n); den = np.maximum(b - a, 1)
    I = np.zeros((len(names), len(g)), np.float32)
    for i, e in enumerate(names):
        if e in tr:
            c = np.concatenate([[0], np.cumsum(tr[e].M, dtype=np.int64)]); I[i] = (c[b] - c[a]) / den
    return g, I


def family_of(families, end):
    """families: {chm13 end -> family token}; ends not listed are their own family (token = arm, e.g. '3q')."""
    return families.get(end, arm_of(end))


def label_grid(I, names, homes, families):
    """per grid point: (label, best end, n_tied)"""
    hset = [i for i, e in enumerate(names) if e in set(homes)]
    best = I.argmax(0); bid = I.max(0); out = []
    for j in range(I.shape[1]):
        if bid[j] < MIN_ID: out.append(('none', '', 0)); continue
        tied = I[:, j] >= bid[j] - MARGIN
        b = names[best[j]]
        if any(tied[i] for i in hset): out.append(('home', b, int(tied.sum())))
        else: out.append((family_of(families, b), b, int(tied.sum())))
    return out


def label_track(tr, label, homes, families, n):
    """per-base match track of a label (max over its member ends)"""
    if label == 'none': return np.zeros(n, np.int8)
    mem = list(homes) if label == 'home' else [e for e in tr if family_of(families, e) == label]
    t = np.zeros(n, np.int8)
    for e in mem:
        if e in tr: np.maximum(t, tr[e].M, out=t)
    return t


def segment(tr, s0, L, homes, families, names):
    """Blocks of one end. Returns list of dicts (label, source, n_grid, sub_start, sub_end, tip_start, tip_end, ident, switch_res,
    chm13_dist_start, chm13_dist_end, mean_tied) in telomere -> centromere order, and the grid labels."""
    n = L - s0
    g, I = identity_grid(tr, s0, L, names)
    labs = label_grid(I, names, homes, families)
    if not len(g): return [], labs, g
    runs = []                                            # [label, first grid idx, last grid idx]
    for j, (lab, b, k) in enumerate(labs):
        if runs and runs[-1][0] == lab: runs[-1][2] = j
        else: runs.append([lab, j, j])
    cache = {}
    def tk(lab):
        if lab not in cache: cache[lab] = label_track(tr, lab, homes, families, n)
        return cache[lab]
    def span(r): return (0 if r[1] == 0 else int(g[r[1]] - GRID // 2)), (n if r[2] == len(g) - 1 else int(g[r[2]] + GRID // 2))
    while len(runs) > 1:
        lens = [span(r)[1] - span(r)[0] for r in runs]; i = int(np.argmin(lens))
        if lens[i] >= MIN_RAW: break
        a, b = span(runs[i]); cand = [k for k in (i - 1, i + 1) if 0 <= k < len(runs)]
        k = max(cand, key=lambda k: (tk(runs[k][0])[a:b].mean() if b > a else 0, -abs(k - i), -k))
        if k < i: runs[k][2] = runs[i][2]
        else: runs[k][1] = runs[i][1]
        del runs[i]
        j = 1
        while j < len(runs):
            if runs[j][0] == runs[j - 1][0]: runs[j - 1][2] = runs[j][2]; del runs[j]
            else: j += 1
    bnd = [span(r)[0] for r in runs[1:]]; res = []
    for i, x0 in enumerate(bnd):                          # base-level refinement of each boundary
        A, B = runs[i][0], runs[i + 1][0]
        lo = max(span(runs[i])[0], x0 - REFINE); hi = min(span(runs[i + 1])[1], x0 + REFINE)
        if 'none' in (A, B) or hi - lo < 200: res.append('grid'); continue
        ea = 1 - tk(A)[lo:hi]; eb = 1 - tk(B)[lo:hi]
        ca = np.concatenate([[0], np.cumsum(ea)]); cb = np.concatenate([[0], np.cumsum(eb)])
        x = int(np.argmin(ca + (cb[-1] - cb)))
        if 0 < x < hi - lo: bnd[i] = lo + x; res.append('base')
        else: res.append('grid')
    edges = [0] + bnd + [n]; blocks = []
    for i, r in enumerate(runs):
        a, b = edges[i], edges[i + 1]; js = range(r[1], r[2] + 1)
        bests = collections.Counter(labs[j][1] for j in js if labs[j][1])
        src = sorted(bests.items(), key=lambda x: (-x[1], x[0]))[0][0] if bests else ''
        idv = float(np.mean([I[:, j].max() for j in js])) if len(js) else float('nan')
        d0 = d1 = None
        if src and src in tr and b > a:
            Ms = tr[src].M[a:b]; ok = np.where(Ms == 1)[0]
            if len(ok):
                t0, t1 = tr[src].tpos(a + int(ok[0])), tr[src].tpos(a + int(ok[-1]))
                if t0 >= 0 and t1 >= 0: d0, d1 = (t0, t1) if src.endswith('p') else (EL - t0, EL - t1)
        blocks.append(dict(label=r[0], source=src, n_grid=r[2] - r[1] + 1, sub_start=a, sub_end=b, tip_start=a + s0, tip_end=b + s0,
                           ident=round(idv, 4), mean_tied=round(float(np.mean([labs[j][2] for j in js])), 2),
                           chm13_dist_start=d0, chm13_dist_end=d1, next_switch_res=res[i] if i < len(res) else ''))
    return blocks, labs, g


def fixed_state(tr, s0, L, a, b, homes, families, names):
    """state of the subtelomere span [a, b) (query coordinates)"""
    b = min(b, L - s0)
    if b - a < 1000: return dict(state='none', best_end='', family='', ident=float('nan'))
    ids = {e: float(tr[e].M[a:b].mean()) for e in names if e in tr}
    if not ids: return dict(state='none', best_end='', family='', ident=float('nan'))
    be = max(sorted(ids), key=ids.get); bi = ids[be]
    if bi < MIN_ID: return dict(state='none', best_end=be, family='', ident=round(bi, 4))
    hmax = max((ids.get(h, 0) for h in homes), default=0)
    st = 'own' if be in homes else ('tie' if hmax >= bi - MARGIN else 'foreign')
    return dict(state=st, best_end=be, family='' if st != 'foreign' else family_of(families, be), ident=round(bi, 4))


def fixed_state_labels(labs, g, I, names, homes, a, b):
    """state of the subtelomere span [a, b) from the grid labels (phase-independent): the majority label of the grid points in the
    span decides - home -> own (home end best at most of those points) or tie; a family -> foreign; none -> none."""
    js = [j for j, x in enumerate(g) if a <= x < b]
    if not js: return dict(state='none', best_end='', family='', ident=float('nan'))
    lab = collections.Counter(labs[j][0] for j in js); top = sorted(lab.items(), key=lambda x: (-x[1], x[0]))[0][0]
    be = collections.Counter(labs[j][1] for j in js if labs[j][1])
    best_end = sorted(be.items(), key=lambda x: (-x[1], x[0]))[0][0] if be else ''
    ident = round(float(np.mean([I[:, j].max() for j in js])), 4)
    if top == 'none': return dict(state='none', best_end=best_end, family='', ident=ident)
    if top == 'home':
        hb = sum(1 for j in js if labs[j][0] == 'home' and labs[j][1] in homes)
        return dict(state='own' if 2 * hb >= lab['home'] else 'tie', best_end=best_end, family='', ident=ident)
    return dict(state='foreign', best_end=best_end, family=top, ident=ident)


def mosaic_end(e, alns, s0, L, grp, families, names):
    """One end: alns for its subtelomere query (window[s0:]); L = window length. Returns (blocks, record)."""
    return mosaic_from_tracks(e, tracks(L - s0, alns), s0, L, grp, families, names)


def shift_tracks(tr, d):
    """Approximate the tracks of window[s0 + d:] from those of window[s0:] without re-alignment: d > 0 drops the first d bases,
    d < 0 prepends |d| unmatched bases. CHM13 coordinates are not carried over."""
    out = {}
    for e, t in tr.items():
        n = len(t.M); u = Track(max(0, n - d)); u.M[:] = t.M[d:] if d >= 0 else np.concatenate([np.zeros(-d, np.int8), t.M]); out[e] = u
    return out


def mosaic_from_tracks(e, tr, s0, L, grp, families, names):
    homes = group_ends(grp); n = L - s0
    blocks, labs, g = segment(tr, s0, L, homes, families, names)
    for i, b in enumerate(blocks): b.update(end_id=e, block=i + 1, n_blocks=len(blocks), is_home=b['label'] == 'home')
    prof = ''.join('O' if l[0] == 'home' else ('.' if l[0] == 'none' else 'F') for l in labs[::10])   # every 1 kb (10 grid points)
    rec = dict(end_id=e, catalogue_arm=grp, s0=s0, n_grid=len(g), n_blocks=len(blocks),
               n_foreign_blocks=sum(1 for b in blocks if b['label'] not in ('home', 'none')), profile=prof,
               structure='|'.join(b['label'] for b in blocks),
               switches_sub=','.join(str(b['sub_start']) for b in blocks[1:]),
               switch_res=','.join(b['next_switch_res'] for b in blocks[:-1]))
    gg, I = identity_grid(tr, s0, L, names)
    for tag, w in FIXED + (('t20', None),):
        bb = w if w else 20000 - s0
        st = fixed_state_labels(labs, gg, I, names, homes, 0, bb) if FIXED_FROM_LABELS else fixed_state(tr, s0, L, 0, bb, homes, families, names)
        rec.update({f'{tag}_state': st['state'], f'{tag}_best_end': st['best_end'], f'{tag}_family': st['family'], f'{tag}_best_ident': st['ident']})
    return blocks, rec


# ---------------------------------------------------------------- paralog families (built once by the catalogue pipeline)
FAM_STEP, FAM_MIN_POINTS, FAM_TAU = 10, 200, 0.75   # tau chosen by the boundary-shift stability test (QC.md)


def tie_counts(I, names, single=None, pair=None):
    """Accumulate co-tie statistics from one end's identity grid, sampled every FAM_STEP grid points (1 kb):
    single[E] = points where E is in the tie set (identity within MARGIN of the best, best >= MIN_ID); pair[(E, F)] likewise."""
    single = single if single is not None else collections.Counter(); pair = pair if pair is not None else collections.Counter()
    for j in range(0, I.shape[1], FAM_STEP):
        bi = I[:, j].max()
        if bi < MIN_ID: continue
        t = sorted(names[i] for i in np.where(I[:, j] >= bi - MARGIN)[0])
        for a in t: single[a] += 1
        for x in range(len(t)):
            for y in range(x + 1, len(t)): pair[(t[x], t[y])] += 1
    return single, pair


def build_families(single, pair, tau=FAM_TAU, min_points=FAM_MIN_POINTS):
    """Paralog families = connected components of the graph joining CHM13 ends E, F when they are co-tied at >= min_points sampled
    points and in >= tau of the points where the rarer of the two is tied (co-tie rate = pair / min(single)).
    Returns (families {end: token}, table rows). Multi-member families are named PF<nn> by size (ties: first member)."""
    adj = collections.defaultdict(set); edges = []
    for (a, b), c in pair.items():
        r = c / max(1, min(single[a], single[b]))
        if c >= min_points and r >= tau: adj[a].add(b); adj[b].add(a); edges.append((a, b, c, round(r, 4)))
    seen, comps = set(), []
    for v in sorted(adj):
        if v in seen: continue
        st, comp = [v], []
        while st:
            x = st.pop()
            if x in seen: continue
            seen.add(x); comp.append(x); st.extend(adj[x] - seen)
        comps.append(sorted(comp, key=lambda e: (len(e), e)))
    comps.sort(key=lambda c: (-len(c), c[0]))
    fam, rows = {}, []
    for k, c in enumerate(comps, 1):
        tok = f'PF{k:02d}'
        for e in c: fam[e] = tok
        es = [x for x in edges if x[0] in c]
        rows.append(dict(family=tok, n_members=len(c), members=','.join(arm_of(e) for e in sorted(c, key=_endkey)),
                         n_edges=len(es), min_cotie_rate=min(x[3] for x in es), median_cotie_rate=float(np.median([x[3] for x in es])),
                         cotie_points=sum(x[2] for x in es), tied_points=sum(single[e] for e in c)))
    return fam, rows, edges


def _endkey(e):
    c = e[3:-1]; return (int(c) if c.isdigit() else 100 + ord(c[0]), e[-1])


NATIVITY_CLASSES = ('native', 'within_family', 'foreign')


def nativity_class(structure, arm, families):
    """Nativity class of a core structure ('|'-separated tokens) of an end or allele of catalogue arm `arm` (an arm or pooled group).
    families: {chm13 end -> family token} (load_families). Tokens 'home' and 'none' (no CHM13 end >= MIN_ID) carry no foreign evidence.
      native        - only home (and none) blocks
      within_family - non-home blocks, all labelled with a paralog family that contains the own arm (or a member of its pooled group)
      foreign       - at least one block labelled with a family or single CHM13 arm that does not contain the own arm
    Returns '' for an empty structure."""
    if not isinstance(structure, str) or not structure: return ''
    own = set(ARM_GROUPS.get(arm, [arm])); members = collections.defaultdict(set)
    for e, f in families.items(): members[f].add(arm_of(e))
    cls = 'native'
    for t in structure.split('|'):
        if t in ('home', 'none', ''): continue
        if not (members.get(t, {t}) & own): return 'foreign'
        cls = 'within_family'
    return cls


def load_families(path):
    """paralog_families.tsv -> {chm13 end: family token}"""
    import csv
    fam = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        for a in r['members'].split(','): fam['chr' + a] = r['family']
    return fam


# ---------------------------------------------------------------- shared helpers for the pipeline and the query tool
def core_input(blocks):
    """mosaic blocks -> objects for alleles.core: home blocks carry their best home end, other blocks their label token (family/arm/none)"""
    from types import SimpleNamespace
    return [SimpleNamespace(is_home=b['label'] == 'home', source=b['source'] if b['label'] == 'home' else b['label'],
                            sub_start=b['sub_start'], sub_end=b['sub_end']) for b in sorted(blocks, key=lambda b: b['block'])]


def align_minimap2(ref_fa, queries, threads=1, exe='minimap2'):
    """queries: [(name, subtelomere)] -> {name: [(qs, qe, strand, tname, ts, te, cs)]} with MM2_TRACK_ARGS"""
    import os, subprocess, tempfile
    out = {}
    with tempfile.TemporaryDirectory(prefix='subtel_atlas_') as td:
        qfa = os.path.join(td, 'q.fa')
        with open(qfa, 'w') as f:
            for n, s in queries: f.write(f'>{n}\n{s}\n')
        p = subprocess.run([exe] + MM2_TRACK_ARGS + ['-t', str(threads), ref_fa, qfa], capture_output=True, text=True)
        if p.returncode != 0: raise RuntimeError('minimap2 failed: ' + p.stderr[-2000:])
        for line in p.stdout.splitlines():
            q, n, qs, qe, st, tn, ts, te, cs = parse_paf(line); out.setdefault(q, []).append((qs, qe, st, tn, ts, te, cs))
    return out


MM_F_ALL_CHAINS = 0x800000      # minimap.h; set by -P

def align_mappy(ref_fa, queries, threads=1):
    """mappy approximation of align_minimap2 (asm20, cs strings). The mappy API has no secondary-to-primary score ratio (-p), so
    with its default (0.8) most paralogous secondary chains are dropped; all chains are kept instead (MM_F_ALL_CHAINS, the -P mode
    of the command line), which yields a superset of the hits of MM2_TRACK_ARGS and the same grid labels (best per point)."""
    import mappy
    from concurrent.futures import ThreadPoolExecutor
    a = mappy.Aligner(ref_fa, preset='asm20', best_n=200, extra_flags=MM_F_ALL_CHAINS)
    def one(nq):
        n, s = nq; buf = mappy.ThreadBuffer()
        return n, [(h.q_st, h.q_en, '+' if h.strand == 1 else '-', h.ctg, h.r_st, h.r_en, h.cs) for h in a.map(s, buf=buf, cs=True)]
    out = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for n, res in ex.map(one, queries):
            if res: out[n] = res
    return out
