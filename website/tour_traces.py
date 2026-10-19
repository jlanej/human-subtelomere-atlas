"""Window-by-window alignment contrast for every tour molecule: is each stretch of the molecule closer to the CHM13 donor end
D or to the CHM13 end of the block after the switch P (= the recipient R except in 5 claims)?

Each molecule (tip frame, as in the release) is cut into WIN bp windows every STEP bp, from the end of its telomere to
MAX_SPAN bp. Each window is aligned on its own to D and to P (minimap2 via mappy, map-ont / map-hifi; best hit covering at
least MIN_COV of the window). Divergence d = 1000 x (1 - matching letters / alignment block length), in differences per
1,000 aligned letters. The contrast is delta = d_P - d_D: positive = closer to the donor end, negative = closer to P. Both
divergences come from the same read letters, so the read's own error rate largely cancels; windows in sequence shared by
D and P give delta near 0. A window that aligns to only one of the two ends gets delta = +/-CAP (drawn at the scale limit).

Outputs:
  docs/evidence/tour/data/claims/<slug>/trace_<n>.bedgraph   delta per 1 kb, on the person's own claimed end window
                                                            (tip frame; window centres projected through the molecule's
                                                            whole-read alignment to that window)
  report/pitch_data/trace_windows.tsv.gz                    every window: read coordinates, d_D, d_P, delta, side of switch
  report/pitch_data/trace_halves.tsv                        per molecule and side, claim-local: 'before' = telomere end to
                                                            the switch interval; 'after' = the first AFTER_MAIN bp past it
                                                            (also 10 and 40 kb, column after_span). Windows are classed
                                                            informative, shared (|delta| < SHARED) or poor (both d > POOR,
                                                            mostly tandem repeats); medians over informative windows.
                                                            Deeper blocks (beyond the claim) are in the traces, not here.
  tour_data.json                                            claims[].traces (file, read, n windows) for the tour pages

Usage (repository root, after website/tour_build.py; TOUR_CACHE as for the builder): python website/tour_traces.py
       [--summary-only]   (re-summarise trace_windows.tsv.gz without realigning)
       [--blocks-only]    (only the whole-block check: the assembled donor block as one piece against D and P)
"""
import os, sys, json, gzip
import numpy as np, pandas as pd, pysam, mappy as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tour_build as TB  # noqa: E402

WIN, STEP, MAX_SPAN, MIN_COV, CAP = 2000, 1000, 150_000, 0.7, 60.0
AFTER_SPANS, AFTER_MAIN = (10_000, 20_000, 40_000), 20_000   # claim-local stretch after the switch (bp)
SHARED, POOR = 2.0, 100.0   # |delta| < SHARED: shared by D and P; both divergences > POOR: poorly aligned to both (repeats)
OUT = TB.OUT; REPO = TB.REPO
PD = os.path.join(REPO, 'report', 'pitch_data')


def best_div(al, w):
    best = None
    for h in al.map(w):
        if (h.q_en - h.q_st) < MIN_COV * len(w): continue
        d = 1000.0 * (1 - h.mlen / h.blen)
        if best is None or d < best: best = d
    return best


def projector(ref_seq, read, preset):
    """read position -> claimed-window position, through the best whole-read alignment (cigar walk)"""
    al = mp.Aligner(seq=ref_seq, preset=preset)
    hs = [h for h in al.map(read) if h.is_primary and h.strand == 1]
    if not hs: return None
    h = max(hs, key=lambda x: x.mlen)
    q, r, qs, rs = h.q_st, h.r_st, [], []
    for n, op in h.cigar:
        if op in (0, 7, 8):
            qs.append(q); rs.append(r); q += n; r += n; qs.append(q); rs.append(r)
        elif op == 1: q += n
        elif op in (2, 3): r += n
    qs, rs = np.array(qs), np.array(rs)
    return lambda x: int(np.interp(x, qs, rs)) if qs[0] <= x <= qs[-1] else None


def main():
    J = json.load(open(os.path.join(OUT, 'tour_data.json')))
    rel = os.path.join(REPO, 'work', 'tour', 'release', 'read_evidence_v1')
    where = {f[:-len('.reads.tsv')]: r for r, _, fs in os.walk(rel) for f in fs if f.endswith('.reads.tsv')}
    lens = None; ends = {}; rows = []
    for c in J['claims']:
        d = where[c['slug']]
        R = pd.read_csv(os.path.join(d, f'{c["slug"]}.reads.tsv'), sep='\t').drop_duplicates('read_id').set_index('read_id')
        want = {m['read_id'] for m in c['mols']}; seqs = {}
        with pysam.AlignmentFile(os.path.join(d, f'{c["slug"]}.chm13.bam')) as b:
            if lens is None: lens = dict(zip(b.references, b.lengths))
            for r in b.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_name in seqs or r.query_name not in want: continue
                seqs[r.query_name] = TB.rc(r.query_sequence) if r.is_reverse else r.query_sequence
        for e in (c['D'], c['P']):
            if e not in ends: ends[e] = TB.fetch_end(e, lens, 260_000)
        own = {x.name: x.sequence for x in pysam.FastxFile(os.path.join(OUT, 'claims', c['slug'], 'combo.fa.gz'))}[c['claimed_window']]
        c['traces'] = []
        for i, m in enumerate(c['mols']):
            x = R.loc[m['read_id']]; t = TB.tip_frame(seqs[m['read_id']], x.platform, x.telo_len)
            preset = 'map-hifi' if str(m['platform']).startswith('HiFi') else 'map-ont'
            AD, AP = (mp.Aligner(seq=ends[e], preset=preset) for e in (c['D'], c['P']))
            proj = projector(own, t, preset)
            st0 = int(x.telo_offset) + int(x.telo_len) + 200; lo, hi = m['switch_lo'], m['switch_hi']
            bg = []
            for s in range(st0, min(len(t), st0 + MAX_SPAN) - WIN + 1, STEP):
                w = t[s:s + WIN]; dD, dP = best_div(AD, w), best_div(AP, w)
                if dD is None and dP is None: delta = np.nan
                elif dD is None: delta = -CAP
                elif dP is None: delta = CAP
                else: delta = max(-CAP, min(CAP, dP - dD))
                side = 'before' if s + WIN <= lo else ('after' if s >= hi else 'switch')
                rows.append(dict(slug=c['slug'], read_id=m['read_id'], platform=m['platform'], dataset=m['dataset'],
                                 win_start=s, win_end=s + WIN, d_D=dD, d_P=dP, delta=delta, side=side,
                                 switch_lo=lo, switch_hi=hi))
                if proj and not np.isnan(delta):
                    a, b2 = proj(s + WIN // 2 - STEP // 2), proj(s + WIN // 2 + STEP // 2)
                    if a is not None and b2 is not None and b2 > a: bg.append((a, b2, delta))
            fn = f'trace_{i + 1}.bedgraph'
            with open(os.path.join(OUT, 'claims', c['slug'], fn), 'w') as fh:
                fh.write(f'track type=bedGraph name="{m["read_id"][:8]} D-vs-P contrast"\n')
                prev = None
                for a, b2, v in bg:
                    if prev is not None and a < prev: a = prev
                    if b2 > a: fh.write(f'{c["claimed_window"]}\t{a}\t{b2}\t{v:.1f}\n'); prev = b2
            c['traces'].append(dict(file=fn, read_id=m['read_id'], dataset=m['dataset'], n_windows=len(bg)))
        print(c['slug'], len(c['mols']), 'molecules', flush=True)
    W = pd.DataFrame(rows)
    os.makedirs(PD, exist_ok=True)
    W.to_csv(os.path.join(PD, 'trace_windows.tsv.gz'), sep='\t', index=False, float_format='%.2f')
    json.dump(J, open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)
    print('windows', len(W)); summarise()


def window_class(W):
    poor = (W.d_D.fillna(1e9) > POOR) & (W.d_P.fillna(1e9) > POOR)
    shared = W.delta.abs() < SHARED
    return np.where(W.delta.isna() | poor, 'poor', np.where(shared, 'shared', 'informative'))


def summarise():
    W = pd.read_csv(os.path.join(PD, 'trace_windows.tsv.gz'), sep='\t'); W['cls'] = window_class(W)
    halves = []
    for span in AFTER_SPANS:
        sel = (W.side == 'before') | ((W.side == 'after') & (W.win_end <= W.switch_hi + span))
        for (slug, rid, side), g in W[sel].groupby(['slug', 'read_id', 'side']):
            inf = g[g.cls == 'informative']; both = inf.dropna(subset=['d_D', 'd_P'])
            md = inf.delta.median() if len(inf) else np.nan
            want = md > 0 if side == 'before' else md < 0
            halves.append(dict(slug=slug, read_id=rid, platform=g.platform.iloc[0], dataset=g.dataset.iloc[0], side=side,
                               after_span=span, n_windows=len(g), n_informative=len(inf), n_shared=int((g.cls == 'shared').sum()),
                               n_poor=int((g.cls == 'poor').sum()), n_only_D=int((inf.d_P.isna()).sum()),
                               n_only_P=int((inf.d_D.isna()).sum()), med_d_D=both.d_D.median() if len(both) else np.nan,
                               med_d_P=both.d_P.median() if len(both) else np.nan, med_delta=md,
                               verdict='not assessable' if not len(inf) else ('as claimed' if want else 'against claim')))
    H = pd.DataFrame(halves); H.to_csv(os.path.join(PD, 'trace_halves.tsv'), sep='\t', index=False, float_format='%.2f')
    for span, h in H.groupby('after_span'):
        print('after span', span, h.groupby(['side', 'verdict']).size().to_dict())
    return H


def block_check():
    """the person's assembled donor block (catalogue block next to the switch, on the claimed window), aligned as one piece
    (minimap2 asm20) to the CHM13 donor end D and to P. Per end: identity of the best hit (matching letters / alignment block
    length, gaps included), matching letters, and the cost of explaining the block with that end = differences per 1,000
    block letters, counting mismatched and gap letters of the hit plus every block letter left unaligned. Stored as
    claims[].block_check in tour_data.json."""
    J = json.load(open(os.path.join(OUT, 'tour_data.json'))); ends = {}
    fai = os.path.join(TB.CACHE, 'chm13', 'chm13v2.0.fa.gz.fai')
    lens = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)} if os.path.exists(fai) else None
    for c in J['claims']:
        d = os.path.join(OUT, 'claims', c['slug'])
        own = {x.name: x.sequence for x in pysam.FastxFile(os.path.join(d, 'combo.fa.gz'))}[c['claimed_window']]
        bl = [l.rstrip('\n').split('\t') for l in open(os.path.join(d, 'own_blocks.bed'))
              if l.startswith(c['claimed_window'] + '\t') and l.split('\t')[3].startswith(c['D'] + '-like')]
        if not bl: c['block_check'] = None; continue
        b = max(bl, key=lambda x: int(x[2]) - int(x[1])); s = own[int(b[1]):int(b[2])]; res = dict(start=int(b[1]), end=int(b[2]))
        for k, e in (('D', c['D']), ('P', c['P'])):
            if e not in ends:
                if lens is None:
                    with pysam.AlignmentFile(os.path.join(REPO, 'work', 'tour', 'release', 'read_evidence_v1', 'HG002_Q100', 'HG002_2_chr3_MATERNAL_R',
                                                          'HG002_2_chr3_MATERNAL_R.chm13.bam')) as bb: lens = dict(zip(bb.references, bb.lengths))
                ends[e] = TB.fetch_end(e, lens, 260_000)
            hs = list(mp.Aligner(seq=ends[e], preset='asm20').map(s))
            h = max(hs, key=lambda x: x.mlen) if hs else None
            res[k] = dict(ident=round(h.mlen / h.blen, 4), matches=h.mlen, r_start=h.r_st, r_end=h.r_en, q_start=h.q_st, q_end=h.q_en,
                          cost=round(1000.0 * ((h.blen - h.mlen) + (len(s) - (h.q_en - h.q_st))) / len(s), 1)) if h \
                else dict(ident=None, matches=0, cost=1000.0)
        c['block_check'] = res
        print('block', c['slug'], res['end'] - res['start'], {k: res[k] and (res[k]['ident'], res[k]['matches']) for k in 'DP'}, flush=True)
    json.dump(J, open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)


def claim_table():
    """report/pitch_data/alignment_goodness_by_claim.tsv: per claim, the window test per side (molecules as claimed /
    assessable; after = first AFTER_MAIN bp) and the whole-block check (matching letters of the assembled donor block to D
    and to P). Also stored as claims[].goodness in tour_data.json for the tour notes."""
    J = json.load(open(os.path.join(OUT, 'tour_data.json')))
    H = pd.read_csv(os.path.join(PD, 'trace_halves.tsv'), sep='\t'); H = H[H.after_span == AFTER_MAIN]
    rows = []
    for c in J['claims']:
        h = H[H.slug == c['slug']]; g = dict(slug=c['slug'], claim_id=c['claim_id'], D=c['D'], P=c['P'], n_molecules=len(c['mols']))
        for side in ('before', 'after'):
            s = h[h.side == side]
            g[f'{side}_as_claimed'] = int((s.verdict == 'as claimed').sum()); g[f'{side}_assessable'] = int((s.verdict != 'not assessable').sum())
        b = c.get('block_check') or {}
        g['block_bp'] = (b['end'] - b['start']) if b else None
        for k in ('D', 'P'):
            for f in ('cost', 'ident', 'matches'): g[f'block_{f}_{k}'] = (b.get(k) or {}).get(f) if b else None
        if b:
            cD, cP = g['block_cost_D'], g['block_cost_P']
            g['block_verdict'] = ('closer to D' if cD < cP - SHARED else ('closer to P' if cP < cD - SHARED else
                                  f'tie (within {SHARED:g} per 1,000)'))
        c['goodness'] = {k: v for k, v in g.items() if k not in ('slug', 'claim_id', 'D', 'P')}
        rows.append(g)
    T = pd.DataFrame(rows); T.to_csv(os.path.join(PD, 'alignment_goodness_by_claim.tsv'), sep='\t', index=False)
    json.dump(J, open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)
    return T


if __name__ == '__main__':
    if '--summary-only' in sys.argv: summarise()
    elif '--blocks-only' in sys.argv: block_check(); claim_table()
    else: main(); block_check(); claim_table()
