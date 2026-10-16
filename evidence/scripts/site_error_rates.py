#!/usr/bin/env python3
"""Per-read error rates at paralog-diagnostic SNV sites (calibration).

For each HG002 claim, every carrier read (reads aligned to the claimed contig end) is evaluated (stage 1 only) and its call at
every diagnostic SNV site is recorded. The haplotype consensus at a site is the majority call (D or P) among all carrier reads of
all platforms covering it (>= 5 reads, >= 80% agreement). A read's site error rate = calls that disagree with the consensus
(wrong allele) or match neither allele (other), over sites with a consensus. Strictly native controls are summarised as the
fraction of D calls at sites inside the claimed donor interval.
Output: <out>/site_error_rates.tsv (per read) and <out>/site_error_summary.tsv (per platform).
"""
import argparse, os, sys, collections, multiprocessing as mp_
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subtel_atlas import readevidence as RE  # noqa: E402
import run_calibration as RC  # noqa: E402


def work(j):
    ctx = RC._ctx(j['claim'], j['preset']); out = []
    for name, seq, rec in RE.iter_reads(j['bam'], min_len=RC.G['min_len']):
        r, u, t = RE.evaluate_read(ctx, name, seq, side=j['side'], meta=dict(platform=j['platform']))
        calls = [(x['d_pos'], x['allele']) for x in u if x['kind'] == 'snv']
        out.append(dict(read=name, platform=j['platform'], role=j['role'], claim=j['claim'], terminal=r['terminal'],
                        s0=r['s0_read'], calls=calls, sub=[x['read_pos'] - r['s0_read'] for x in u if x['kind'] == 'snv']))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--chm13-ends', required=True); ap.add_argument('--work', default=os.path.join(REPO, 'work', 'readevidence'))
    ap.add_argument('--out', required=True); ap.add_argument('-t', '--threads', type=int, default=8)
    a = ap.parse_args()
    import pandas as pd, csv
    crit = RE.load_criteria(a.criteria); ends = RE.Chm13Ends(a.chm13_ends)
    cat = os.path.join(REPO, 'data', 'catalogue'); hg = os.path.join(a.work, 'hg002_ends.tsv')
    claims = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'),
                                      end_ids=pd.read_csv(hg, sep='\t').end_id.tolist())
    k1g = [dict(r, role='ont1kg_native_strict') for r in csv.DictReader(open(os.path.join(a.work, 'ont1kg_strict_native_controls.tsv')), delimiter='\t')]
    jobs = [j for j in RC.build_jobs(claims, hg, k1g, a.work, 'all', RC.split(claims)) if j['role'] in ('carrier', 'hg002_otherhap_native', 'ont1kg_native_strict')]
    cmap = {c['claim_id']: c for c in claims}
    with mp_.get_context('fork').Pool(a.threads, initializer=RC.init, initargs=(cmap, ends, crit, os.path.join(a.work, 'masks'), 5000)) as pool:
        res = [x for out in pool.map(work, jobs, chunksize=1) for x in out]
    cons = {}
    for cid in cmap:
        votes = collections.defaultdict(collections.Counter)
        for r in res:
            if r['claim'] == cid and r['role'] == 'carrier':
                for d, al in r['calls']: votes[d][al] += 1
        for d, v in votes.items():
            n = v['D'] + v['P']; top, k = max((('D', v['D']), ('P', v['P'])), key=lambda x: x[1])
            if n >= 5 and k / n >= 0.8: cons[(cid, d)] = top
    rows = []
    for r in res:
        c = cmap[r['claim']]; sw = c.get('switch_sub')
        if r['role'] == 'carrier':
            ev = [(d, al) for d, al in r['calls'] if (r['claim'], d) in cons]
            wrong = sum(al in ('D', 'P') and al != cons[(r['claim'], d)] for d, al in ev); other = sum(al == 'other' for d, al in ev)
            rows.append(dict(claim=r['claim'], role=r['role'], platform=r['platform'], read=r['read'], n_sites=len(ev), wrong_allele=wrong,
                             other=other, error_rate=(wrong + other) / len(ev) if ev else None))
        else:
            ev = [al for (d, al), s in zip(r['calls'], r['sub']) if sw is None or s < sw]
            nd = sum(al == 'D' for al in ev); npp = sum(al == 'P' for al in ev); no = sum(al == 'other' for al in ev)
            rows.append(dict(claim=r['claim'], role=r['role'], platform=r['platform'], read=r['read'], n_sites=len(ev), D_calls=nd, P_calls=npp,
                             other=no, D_frac=nd / (nd + npp) if nd + npp else None))
    df = pd.DataFrame(rows); os.makedirs(a.out, exist_ok=True)
    df.to_csv(os.path.join(a.out, 'site_error_rates.tsv'), sep='\t', index=False)
    car = df[(df.role == 'carrier') & (df.n_sites >= 20)]
    s1 = car.groupby('platform').apply(lambda g: pd.Series(dict(reads=len(g), sites=int(g.n_sites.sum()),
                                       wrong_allele_rate=g.wrong_allele.sum() / g.n_sites.sum(), other_rate=g.other.sum() / g.n_sites.sum(),
                                       median_read_error=g.error_rate.median(), p95_read_error=g.error_rate.quantile(0.95))))
    ctl = df[(df.role != 'carrier') & (df.n_sites >= 20)]
    s2 = ctl.groupby(['role', 'platform']).apply(lambda g: pd.Series(dict(reads=len(g), sites=int(g.n_sites.sum()),
                                                 D_call_rate=g.D_calls.sum() / max(1, g.D_calls.sum() + g.P_calls.sum()),
                                                 other_rate=g.other.sum() / g.n_sites.sum(), max_read_D_frac=g.D_frac.max())))
    s1.reset_index().assign(kind='carrier_vs_consensus').to_csv(os.path.join(a.out, 'site_error_summary.tsv'), sep='\t', index=False)
    s2.reset_index().assign(kind='strict_native_control_donor_interval').to_csv(os.path.join(a.out, 'site_error_summary.tsv'), sep='\t',
                                                                                 index=False, mode='a')
    print(s1.to_string()); print(s2.to_string())


if __name__ == '__main__':
    main()
