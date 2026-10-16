"""An end as an ordered mosaic of CHM13 end segments, telomere -> centromere (catalogue step 2c).

Used unchanged by pipeline/catalogue_06_mosaic.py and by the query tool (mosaic_end()).
Bin labels (5 kb bins from the telomere boundary): home (a home CHM13 end is best or within 2% of the best AS), <end> (another CHM13
end is best and every home end scores < 0.98 x best: a 'foreign' bin), none (no alignment covering >= 50% of the bin). Isolated bins
between two identical labels are relabelled when that label scores >= 0.90 x best (isolated 'none' bins are absorbed).
Blocks = maximal runs of one label; each boundary between two aligned blocks is refined at base level (mappy asm20 against both
source ends, per-base error tracks, changepoint minimising outer errors before + inner errors after), else kept at bin resolution.
Fixed-window states S5/S10/S20 (and t20): own / tie (home within 2% of the best foreign end) / foreign / none.
A 'switch' is a change in which CHM13 end the sequence most resembles inside shared, paralogous sequence, not necessarily a structural
breakpoint; labels are relative to CHM13.
"""
import re, collections
import numpy as np
from .common import arm_of, group_ends

NEAR, SMOOTH, EL, BIN = 0.98, 0.90, 500_000, 5000
FIX = ('S5', 'S10', 'S20', 't20')


def parse_cand(c):
    if not isinstance(c, str) or not c: return {}
    out = {}
    for x in c.split(';'):
        e, a, i = x.split(':'); out[e] = (int(a), float(i))
    return out


def label_bins(g, homes):
    labs = []
    for r in g.itertuples():
        cd = parse_cand(r.cand)
        if not cd: labs.append(['none', cd, r, None]); continue
        best = max(a for a, i in cd.values()); hs = [h for h in homes if h in cd and cd[h][0] >= NEAR * best]
        if hs:
            h = max(hs, key=lambda x: (cd[x][0], x)); cd = dict(cd); cd['home'] = cd[h]
            labs.append(['home', cd, r, h])
        else:
            hh = [h for h in homes if h in cd]
            if hh: cd = dict(cd); cd['home'] = max((cd[h] for h in hh), key=lambda v: v[0])
            labs.append([r.best_end, cd, r, None])
    nsm = 0
    for i in range(1, len(labs) - 1):
        y = labs[i - 1][0]
        if labs[i + 1][0] == y and labs[i][0] != y and y != 'none':
            cd = labs[i][1]; best = max((a for a, _ in cd.values()), default=0)
            if labs[i][0] == 'none' or (y in cd and cd[y][0] >= SMOOTH * best): labs[i][0] = y; nsm += 1
    return labs, nsm


def blocks_of(labs, s0):
    out = []
    for lab, cd, r, he in labs:
        if out and out[-1]['label'] == lab:
            b = out[-1]
        else:
            b = dict(label=lab, tip_start=int(r.qstart), idents=[], d0=[], d1=[], n_bins=0, hends=collections.Counter()); out.append(b)
        b['tip_end'] = int(r.qend); b['n_bins'] += 1
        if lab in cd: b['idents'].append(cd[lab][1])
        if lab == 'home' and he: b['hends'][he] += 1
        src = he if lab == 'home' else lab
        if src == r.best_end and r.best_dist0 == r.best_dist0 and str(r.best_dist0) != '':
            b['d0'].append(float(r.best_dist0)); b['d1'].append(float(r.best_dist1))
    for b in out:
        b['ident'] = round(float(np.mean(b['idents'])), 4) if b['idents'] else np.nan
        b['chm13_dist_start'] = int(min(b['d0'])) if b['d0'] else np.nan
        b['chm13_dist_end'] = int(max(b['d1'])) if b['d1'] else np.nan
        b['sub_start'] = b['tip_start'] - s0; b['sub_end'] = b['tip_end'] - s0
        b['source'] = sorted(b['hends'].items(), key=lambda x: (-x[1], x[0]))[0][0] if b['label'] == 'home' and b['hends'] else b['label']
        del b['idents'], b['d0'], b['d1'], b['hends']
    return out


class RefEnds:
    """The 48 CHM13 v2.0 500 kb ends (genome orientation) with lazily built per-end mappy asm20 aligners."""
    def __init__(self, seqs):
        self.seqs = dict(seqs); self._al = {}

    def aligner(self, end):
        import mappy
        if end not in self._al: self._al[end] = mappy.Aligner(seq=self.seqs[end], preset='asm20')
        return self._al[end]


def err_track(refs, ends, seq):
    n = len(seq); best_e = np.ones(n, np.int8); tmap = []
    for end, h in ((e, h) for e in ends for h in refs.aligner(e).map(seq, cs=True)):
        if h.q_en - h.q_st < 1000: continue
        e = np.ones(n, np.int8); tp = np.full(n, -1, np.int64); t = h.r_st
        q = h.q_st if h.strand == 1 else h.q_en - 1; st = 1 if h.strand == 1 else -1
        for op, val in re.findall(r'([:*+\-~])([0-9]+|[a-z]+)', h.cs):
            if op == ':':
                l = int(val); idx = np.arange(q, q + st * l, st); e[idx] = 0; tp[idx] = np.arange(t, t + l); q += st * l; t += l
            elif op == '*':
                if 0 <= q < n: e[q] = 1; tp[q] = t
                q += st; t += 1
            elif op == '+':
                l = len(val); idx = np.arange(q, q + st * l, st); e[idx] = 1; tp[idx] = t; q += st * l
            elif op == '-':
                if 0 <= q < n: e[q] = 1
                t += len(val)
        tmap.append((e, tp, end)); best_e = np.minimum(best_e, e)
    return best_e, tmap


def dist_from_tel(end, t): return int(t) if end.endswith('p') else EL - int(t)


def refine(refs, seq, A, B, homes):
    """base-level switch between outer block A and inner block B (tip-frame coordinates)"""
    bnd = B['tip_start']; lo = max(A['tip_start'], bnd - 10000); hi = min(B['tip_end'], bnd + 10000)
    sub = seq[lo:hi]
    if 'none' in (A['label'], B['label']) or len(sub) < 2000: return dict(switch_tip=bnd, switch_res='window')
    ends = lambda b: homes if b['label'] == 'home' else [b['label']]
    ea, ta = err_track(refs, ends(A), sub); eb, tb = err_track(refs, ends(B), sub)
    ca = np.concatenate([[0], np.cumsum(ea)]); cb = np.concatenate([[0], np.cumsum(eb)])
    cost = ca + (cb[-1] - cb); x = int(np.argmin(cost)); n = len(sub)
    ia = 1 - ea[max(0, x - 5000):x].mean() if x > 0 else 0; ib = 1 - eb[x:min(n, x + 5000)].mean() if x < n else 0
    res = dict(switch_tip=bnd, switch_res='window', sw_outer_ident_5kb=round(float(ia), 4), sw_inner_ident_5kb=round(float(ib), 4))
    if ia >= 0.95 and ib >= 0.95 and abs(lo + x - bnd) <= 7500 and 0 < x < n:
        res.update(switch_tip=lo + x, switch_res='base')
        ka = [(e[max(0, x - 2000):x].sum(), tp[x - 1], en) for e, tp, en in ta if tp[x - 1] >= 0]
        kb = [(e[x:x + 2000].sum(), tp[x], en) for e, tp, en in tb if tp[x] >= 0]
        if ka: k = min(ka, key=lambda v: (v[0], v[2])); res.update(sw_outer_chm13_end=k[2], sw_outer_chm13_dist=dist_from_tel(k[2], k[1]))
        if kb: k = min(kb, key=lambda v: (v[0], v[2])); res.update(sw_inner_chm13_end=k[2], sw_inner_chm13_dist=dist_from_tel(k[2], k[1]))
    return res


def t20_call(r, homes):
    cd = parse_cand(r.cand) if r is not None else {}
    if not cd: return dict(t20_state='none', t20_best_end='')
    b = r.best_end; best = cd[b][0]; hAS = max((cd[h][0] for h in homes if h in cd), default=0)
    if b in homes: st = 'own'
    elif hAS >= NEAR * best: st = 'tie'
    else: st = 'foreign'
    return dict(t20_state=st, t20_best_end=b, t20_best_ident=cd[b][1], t20_home_AS=hAS, t20_best_AS=best)


def mosaic_end(refs, e, g, fixed, s0, grp, seq):
    """Mosaic of one end. g = its tiled-bin rows (DataFrame sorted by bin number), fixed = {tag: bins row or None} for S5/S10/S20/t20,
    s0 = boundary, grp = catalogue arm, seq = tip-frame window. Returns (block dicts, per-end record) as in catalogue_06."""
    homes = group_ends(grp)
    labs, nsm = label_bins(g, homes); bl = blocks_of(labs, s0)
    sw = []
    for i in range(1, len(bl)):
        rf = refine(refs, seq, bl[i - 1], bl[i], homes); rf['switch_sub'] = rf['switch_tip'] - s0; sw.append(rf)
    BR = []
    for i, b in enumerate(bl):
        b.update(end_id=e, block=i + 1, n_blocks=len(bl), is_home=b['label'] == 'home')
        if i > 0: b['tip_start'] = sw[i - 1]['switch_tip']; b['sub_start'] = sw[i - 1]['switch_sub']
        if i < len(bl) - 1: b['tip_end'] = sw[i]['switch_tip']; b['sub_end'] = sw[i]['switch_sub']
        if i < len(bl) - 1: b.update({f'next_{k}': v for k, v in sw[i].items()})
        BR.append(b)
    rec = dict(end_id=e, catalogue_arm=grp, home_ends=','.join(homes), s0=s0, n_bins=len(g), n_smoothed=nsm, n_blocks=len(bl),
               n_foreign_blocks=sum(1 for b in bl if not b['is_home'] and b['source'] != 'none'),
               profile=''.join('O' if l[0] == 'home' else ('.' if l[0] == 'none' else 'F') for l in labs),
               structure='|'.join('home' if b['label'] == 'home' else arm_of(b['label']) for b in bl),
               switches_sub=','.join(str(x['switch_sub']) for x in sw), switch_res=','.join(x['switch_res'] for x in sw))
    for t in FIX:
        rec.update({k.replace('t20', t): v for k, v in t20_call(fixed.get(t), homes).items()})
    return BR, rec
