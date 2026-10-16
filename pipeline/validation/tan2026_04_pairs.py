"""Step 4: same-haplotype cross-check of HPRC Release 2 catalogue ends against the HPRCY2 records of Tan et al. 2026.
(a) Best match. For every Release 2 end (data/catalogue/ends.tsv, set HPRC_R2) of a sample present in the HPRCY2 records, the
    catalogue window (work/windows/ends.fa, tip frame) is mapped (mappy asm5) to all HPRCY2 windows of the same sample (both
    haplotypes). Best match = target with the largest summed matching bases (mlen) of primary forward-strand alignments.
(b) Haplotype labels. Haplotype numbering is not consistent between the catalogue assemblies and the Tan HPRCY2 (mostly Verkko)
    records outside trio-phased samples (labels can differ chromosome by chromosome), so labels are not used for pairing.
(c) Key pair = catalogue named end vs the Tan record of the same sample and arm (catalogue_arm vs the Tan arm label) with the
    most matching bases, either haplotype label (the best match of (a) when it is such a record); hap_label records whether
    the labels agree (same/opposite). hap_map.tsv lists the label relation per sample and chromosome.
    For each pair a single-target alignment (both sequences from 500 bp before their telomere boundary, so that the long
    telomere arrays do not seed spurious chains; Tan telomeres are trimmed to 2 kb) gives: identity = sum(mlen)/sum(blen) over
    primary forward alignments, i.e. essentially subtelomere identity; aligned
    fraction of the catalogue subtelomere (tip[s0:]); boundary projection: the catalogue position s0 + 1000 is projected through
    the alignment into the Tan window and 1000 subtracted (if s0 + 1000 is unaligned, the first aligned base within 5 kb of
    the boundary is used; boundary_anchor = its distance from the boundary); boundary_delta = projected - s0_tan (0 = same telomere/subtelomere
    boundary); diag_offset = Tan tip coordinate - catalogue tip coordinate at that point. Catalogue switches (mosaic and core,
    subtelomere coordinates) are projected the same way into Tan subtelomere coordinates.
Outputs: work/tan2026/best_match.tsv (cached; partner search only), pairs_all.tsv (cache of every aligned pair, keyed by its
         inputs; re-runs align only pairs whose catalogue s0/switches changed), pairs.tsv (selected pair per end), hap_map.tsv
Usage  : python tan2026_04_pairs.py
"""
import os, sys, collections, multiprocessing as mp, numpy as np, pandas as pd, pysam, mappy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
from catalogue_common import WIN

def union_len(iv, lo, hi):
    iv = sorted((max(a, lo), min(b, hi)) for a, b in iv if min(b, hi) > max(a, lo)); tot, cur = 0, None
    for a, b in iv:
        if cur is None or a > cur[1]:
            if cur: tot += cur[1] - cur[0]
            cur = [a, b]
        else: cur[1] = max(cur[1], b)
    return tot + (cur[1] - cur[0] if cur else 0)

def project(hits, q):
    for h in hits:
        if h.strand != 1 or not (h.q_st <= q < h.q_en): continue
        qp, tp = h.q_st, h.r_st
        for l, op in h.cigar:
            if op in (0, 7, 8):
                if qp + l > q: return tp + (q - qp)
                qp += l; tp += l
            elif op == 1:
                if qp + l > q: return tp
                qp += l
            elif op == 2: tp += l
    return None

def fwd_primary(hits): return [h for h in hits if h.is_primary and h.strand == 1]

_FA = {}
def fa(kind):
    if kind not in _FA: _FA[kind] = pysam.FastaFile(f'{WIN}/ends.fa' if kind == 'ours' else f'{T.TW}/tan_tips.fa')
    return _FA[kind]

def best_matches(job):
    sample, ours, tans = job
    tseq = {t: fa('tan').fetch(t) for t in tans}
    # multi-sequence index: a small per-worker FASTA, overwritten for each sample
    tmp = f'{T.TW}/pairs_tmp/worker{os.getpid() % 1000}.fa'
    with open(tmp, 'w') as f:
        for t, s in tseq.items(): f.write(f'>{t}\n{s}\n')
    al = mappy.Aligner(tmp, preset='asm5'); out = []
    for e, s0 in ours:
        q = fa('ours').fetch(e); hs = fwd_primary(al.map(q)); by = collections.defaultdict(list)
        for h in hs: by[h.ctg].append(h)
        if not by: out.append(dict(end_id=e, best_tan='', best_mlen=0)); continue
        t, H = max(by.items(), key=lambda kv: (sum(h.mlen for h in kv[1]), kv[0]))
        sec = sorted((sum(h.mlen for h in v) for k, v in by.items() if k != t), reverse=True)
        out.append(dict(end_id=e, best_tan=t, best_mlen=sum(h.mlen for h in H), best_ident=round(sum(h.mlen for h in H) / sum(h.blen for h in H), 5),
                        best_sub_frac=round(union_len([(h.q_st, h.q_en) for h in H], s0, len(q)) / max(1, len(q) - s0), 4),
                        second_mlen=sec[0] if sec else 0))
    return out

PAD = 500   # bp of telomere kept before each boundary in the pairwise alignment

def pair_metrics(job):
    e, t, s0, s0t, sws, csws = job
    q = fa('ours').fetch(e); r = fa('tan').fetch(t); qo, ro = max(0, s0 - PAD), max(0, s0t - PAD)
    hs = fwd_primary(mappy.Aligner(seq=r[ro:], preset='asm5').map(q[qo:]))
    L = len(q); rec = dict(end_id=e, tan_end_id=t)
    if not hs: return dict(rec, ident=np.nan, sub_aligned_frac=0.0, n_aln=0, mlen=0)
    qi = [(h.q_st + qo, h.q_en + qo) for h in hs]; ri = [(h.r_st + ro, h.r_en + ro) for h in hs]
    rec.update(mlen=sum(h.mlen for h in hs), ident=round(sum(h.mlen for h in hs) / sum(h.blen for h in hs), 6), n_aln=len(hs),
               nm_per_kb=round(1000 * sum(h.NM for h in hs) / max(1, sum(h.blen for h in hs)), 4),
               sub_aligned_frac=round(union_len(qi, s0, L) / max(1, L - s0), 4),
               tan_sub_aligned_frac=round(union_len(ri, s0t, min(len(r), s0t + (L - s0))) / max(1, min(len(r) - s0t, L - s0)), 4))
    proj = lambda x: (lambda y: None if y is None else y + ro)(project(hs, x - qo))
    anc = s0 + 1000; p = proj(anc)
    if p is None:   # fallback: first aligned catalogue base within 5 kb of the boundary (+50 bp inside the alignment)
        st = sorted(a + 50 for a, b in qi if b > max(a, s0) + 100 and a + 50 >= s0 and a + 50 < s0 + 5000)
        if st: anc = st[0]; p = proj(anc)
    if p is not None: rec.update(boundary_delta=p - (anc - s0) - s0t, diag_offset=p - anc, boundary_anchor=anc - s0)
    for tag, sw in (('switches', sws), ('core_switches', csws)):
        pr = []
        for x in sw:
            y = proj(s0 + x); pr.append('' if y is None else str(y - s0t))
        rec[f'{tag}_projected'] = ','.join(pr)
    return rec

def ints(s): return [int(float(x)) for x in str(s).split(',') if x not in ('', 'nan')] if isinstance(s, str) else []

def main():
    os.makedirs(f'{T.TW}/pairs_tmp', exist_ok=True)
    if not os.path.exists(f'{T.TW}/tan_tips.fa.fai'): pysam.faidx(f'{T.TW}/tan_tips.fa')
    E = pd.read_csv(f'{T.REPO}/data/catalogue/ends.tsv', sep='\t', low_memory=False)
    E = E[E.set == 'HPRC_R2'].copy(); E['s0'] = (E.junk_len + E.telo_len).astype(int)
    TW = pd.read_csv(f'{T.TW}/tan_windows.tsv', sep='\t'); H = TW[TW.cohort == 'HPRCY2']
    samples = sorted(set(E['sample']) & set(H['sample']))
    bm_f = f'{T.TW}/best_match.tsv'
    if not os.path.exists(bm_f):
        jobs = [(s, list(zip(E[E['sample'] == s].end_id, E[E['sample'] == s].s0)), list(H[H['sample'] == s].end_id)) for s in samples]
        with mp.Pool(T.THREADS) as p: BM = pd.DataFrame([r for rr in p.imap_unordered(best_matches, jobs) for r in rr])
        BM.to_csv(bm_f, sep='\t', index=False)
    BM = pd.read_csv(bm_f, sep='\t')
    # key pairs: each named catalogue end vs the Tan HPRCY2 record(s) of the same sample and arm, either haplotype label.
    # If the best match (a) is such a record it is used directly; otherwise every same-arm record of the sample is aligned and
    # the one with most matching bases is kept.
    N = E[E.arm_confidence == 'named'].copy()
    K = N.merge(H[['end_id', 'sample', 'hap', 'tan_arm', 'junk_len', 'telo_len']].rename(columns={'end_id': 'tan_end_id', 'hap': 'tan_hap',
                'tan_arm': 'catalogue_arm', 'junk_len': 'tjl', 'telo_len': 'ttl'}), on=['sample', 'catalogue_arm'])
    K['s0t'] = (K.tjl + K.ttl).astype(int)
    bt = dict(zip(BM.end_id, BM.best_tan)); K['is_best'] = [bt.get(e) == t for e, t in zip(K.end_id, K.tan_end_id)]
    hb = set(K.end_id[K.is_best]); K = K[K.is_best | ~K.end_id.isin(hb)]
    nd = len(set(zip(N['sample'], N.pansn_hap, N.catalogue_arm)) & set(zip(H['sample'], H.hap, H.tan_arm)))
    # incremental: a cached pair is reused when its inputs (sequences fixed; s0 values and catalogue switches) are unchanged
    K['sig'] = [f'{a}|{b}|{c}|{d}' for a, b, c, d in zip(K.s0, K.s0t, K.switches_sub.fillna(''), K.core_switches_sub.fillna(''))]
    cache_f = f'{T.TW}/pairs_all.tsv'
    cache = pd.read_csv(cache_f, sep='\t', dtype={'switches_projected': str, 'core_switches_projected': str, 'sig': str}) if os.path.exists(cache_f) else pd.DataFrame(columns=['end_id', 'tan_end_id', 'sig'])
    have = set(zip(cache.end_id, cache.tan_end_id, cache.sig)); todo = K[[(a, b, c) not in have for a, b, c in zip(K.end_id, K.tan_end_id, K.sig)]]
    jobs = [(r.end_id, r.tan_end_id, r.s0, r.s0t, ints(r.switches_sub), ints(r.core_switches_sub)) for r in todo.itertuples()]
    print('pairs to align', len(jobs), 'reused', len(K) - len(jobs), flush=True)
    new = []
    if jobs:
        with mp.Pool(T.THREADS) as p: new = list(p.imap_unordered(pair_metrics, jobs, chunksize=20))
    New = pd.DataFrame(new)
    if len(New): New = New.merge(todo[['end_id', 'tan_end_id', 'sig']], on=['end_id', 'tan_end_id'])
    Cn = pd.concat([cache, New], ignore_index=True).drop_duplicates(['end_id', 'tan_end_id', 'sig'], keep='last')
    Cn.to_csv(cache_f + '.tmp', sep='\t', index=False); os.replace(cache_f + '.tmp', cache_f)
    P = K[['end_id', 'tan_end_id', 'sig']].merge(Cn, on=['end_id', 'tan_end_id', 'sig']).drop(columns=['sig'])
    P = P.merge(K[['end_id', 'tan_end_id', 'pansn_hap', 'tan_hap', 'is_best']], on=['end_id', 'tan_end_id'])
    P['n_candidates'] = P.groupby('end_id').tan_end_id.transform('size')
    P = P.sort_values(['end_id', 'mlen', 'tan_end_id'], ascending=[True, False, True]).drop_duplicates('end_id')
    P['hap_label'] = np.where(P.pansn_hap == P.tan_hap, 'same', 'opposite')
    P = P.sort_values(['end_id', 'tan_end_id']); P.to_csv(f'{T.TW}/pairs.tsv', sep='\t', index=False)
    # haplotype-label relation per sample and chromosome (verified pairs only)
    v = P[(P.ident >= 0.99) & (P.sub_aligned_frac >= 0.5)].merge(E[['end_id', 'sample', 'catalogue_arm']], on='end_id')
    v['chrom'] = v.catalogue_arm.str[:-1]
    hm = v.groupby(['sample', 'chrom']).hap_label.agg(lambda x: ','.join(sorted(set(x)))).reset_index()
    hm.to_csv(f'{T.TW}/hap_map.tsv', sep='\t', index=False)
    print('samples', len(samples), 'direct-label keys', nd, 'pairs', len(P), P.hap_label.value_counts().to_dict(),
          'chromosomes with mixed p/q label relation', int(hm.hap_label.str.contains(',').sum()), 'of', len(hm))

if __name__ == '__main__':
    main()
