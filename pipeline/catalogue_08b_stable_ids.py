"""Step 2g: keep structural-allele and haplogroup IDs stable across rebuilds.
Steps 07/08 number clusters by size (SA-<arm>-NN / HG-<arm>-NN). This step maps them onto the IDs of the previous release
(ATLAS_PREVIOUS_ENDS, default data/catalogue/ends.tsv as it exists before step 09 rewrites it; skipped if absent):
  * per arm, member ends (HPRC_R2/CHM13/HG002_Q100, any assignment) are cross-tabulated new cluster x previous ID;
  * a previous ID is inherited by a new cluster when they are reciprocal best matches (the cluster holds the largest share of
    the ID's previous members, and the ID is the largest previous source of the cluster's members; ties by ID); for structural
    alleles only between clusters of the same nativity (a native ID is never given to a non-native allele, or vice versa);
  * other clusters get new IDs numbered after the highest number ever used in the arm.
Status per ID: unchanged (same members), kept_changed (same ID, members changed), new; previous IDs not inherited are
retired; 'merged' marks clusters receiving >= 2 previous IDs' members (each >= 2 ends), 'split' previous IDs whose
members went to >= 2 clusters (each >= 2 ends).
Rewrites work/catalogue/{sa_members,sa_alleles,haplogroup_members,haplogroups_work}.tsv in place and writes
work/catalogue/id_history.tsv (copied to the release by step 09).
"""
import os, sys, re, collections, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
PREV = os.environ.get('ATLAS_PREVIOUS_ENDS', f'{C.CAT}/ends.tsv')
CORE = {'HPRC_R2', 'CHM13', 'HG002_Q100'}
num = lambda i: int(i.rsplit('-', 1)[1])

def stable(kind, members, idcol, prev_of, same_class=None):
    """members: end_id, arm, <idcol>, set. Returns mapping new->stable id and history rows."""
    M = members[members.set.isin(CORE) & members[idcol].fillna('').ne('')]
    hist, mp = [], {}
    for arm, g in M.groupby('arm', sort=True):
        g = g.assign(prev=g.end_id.map(prev_of).fillna(''))
        ct = collections.Counter(zip(g[idcol], g.prev))
        ctm = ct if same_class is None else collections.Counter({(n, p): c for (n, p), c in ct.items() if not p or same_class(n, p)})
        # ctm (same-class pairs) decides inheritance; ct (all member moves) is reported in the history
        new_ids = sorted(set(g[idcol]), key=num)
        prev_ids = sorted({p for p in g.prev if p}, key=num)
        prev_all = {p for p in prev_of.values() if isinstance(p, str) and p.startswith(f'{kind}-{arm}-')}
        best_prev = {n: max(((c, -num(p), p) for (n2, p), c in ctm.items() if n2 == n and p), default=(0, 0, ''))[2] for n in new_ids}
        best_new = {p: max(((c, -num(n), n) for (n, p2), c in ctm.items() if p2 == p), default=(0, 0, ''))[2] for p in prev_ids}
        used = set(); top = max([num(p) for p in prev_all] or [0])
        for n in new_ids:
            p = best_prev[n]
            if p and best_new.get(p) == n and p not in used: mp[n] = p; used.add(p)
        for n in new_ids:
            if n not in mp: top += 1; mp[n] = f'{kind}-{arm}-{top:02d}'
        prev_members = collections.defaultdict(set); new_members = collections.defaultdict(set)
        for r in g.itertuples():
            new_members[mp[getattr(r, idcol)]].add(r.end_id)
            if r.prev: prev_members[r.prev].add(r.end_id)
        for n in new_ids:
            sid = mp[n]; src = collections.Counter(g.prev[g[idcol] == n]); src.pop('', None)
            status = 'new' if sid not in prev_all else ('unchanged' if new_members[sid] == prev_members.get(sid, set()) else 'kept_changed')
            merged = sum(1 for c in src.values() if c >= 2) >= 2
            hist.append(dict(type=kind, arm=arm, id=sid, status=status, merged=merged, n_members=len(new_members[sid]),
                             previous_sources=';'.join(f'{p}:{c}' for p, c in sorted(src.items(), key=lambda x: (-x[1], x[0])))))
        for p in sorted(prev_all, key=num):
            dest = collections.Counter(mp[n] for (n, p2), c in ct.items() if p2 == p for _ in range(c))
            split = sum(1 for c in dest.values() if c >= 2) >= 2
            if p not in used:
                hist.append(dict(type=kind, arm=arm, id=p, status='retired', merged=False, n_members=0, split=split,
                                 previous_sources='', now_in=';'.join(f'{k}:{v}' for k, v in sorted(dest.items()))))
            else:
                for h in hist:
                    if h['id'] == p and h['type'] == kind: h['split'] = split; h['now_in'] = ';'.join(f'{k}:{v}' for k, v in sorted(dest.items()))
    return mp, hist

def main():
    if not os.path.exists(PREV):
        print('no previous release table; IDs left as numbered by size'); return
    P = pd.read_csv(PREV, sep='\t', low_memory=False, usecols=['end_id', 'structural_allele_id', 'haplogroup_id'])
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t', usecols=['end_id', 'set'])
    SM = pd.read_csv(f'{C.CW}/sa_members.tsv', sep='\t'); SA = pd.read_csv(f'{C.CW}/sa_alleles.tsv', sep='\t')
    HM = pd.read_csv(f'{C.CW}/haplogroup_members.tsv', sep='\t'); HW = pd.read_csv(f'{C.CW}/haplogroups_work.tsv', sep='\t')
    # always map both types against the previous release: re-mapping IDs that are already stable is an identity (status unchanged)
    prevA = os.path.join(os.path.dirname(PREV), 'alleles_structural.tsv')
    pnat = dict(pd.read_csv(prevA, sep='\t', usecols=['allele_id', 'is_native']).values) if os.path.exists(prevA) else {}
    nnat = dict(zip(SA.allele_id, SA.is_native))
    mp_sa, h1 = stable('SA', SM, 'allele_id', dict(zip(P.end_id, P.structural_allele_id)),
                       same_class=lambda n, p: bool(nnat.get(n)) == bool(pnat.get(p, nnat.get(n))))
    hm = HM.merge(W, on='end_id', how='left')
    mp_hg, h2 = stable('HG', hm, 'haplogroup_id', dict(zip(P.end_id, P.haplogroup_id)))
    f = lambda m: (lambda x: m.get(x, x) if isinstance(x, str) and x else x)
    SM['allele_id'] = SM.allele_id.map(f(mp_sa))
    SM['candidates'] = SM.candidates.map(lambda c: ','.join(mp_sa.get(i, i) for i in c.split(',')) if isinstance(c, str) and c else c)
    SA['allele_id'] = SA.allele_id.map(f(mp_sa)); SA['id_stabilised'] = True
    HM['haplogroup_id'] = HM.haplogroup_id.map(f(mp_hg)); HW['haplogroup_id'] = HW.haplogroup_id.map(f(mp_hg)); HW['id_stabilised'] = True
    assert SA.allele_id.is_unique and HW.haplogroup_id.is_unique
    for d, fn in ((SM, 'sa_members'), (SA, 'sa_alleles'), (HM, 'haplogroup_members'), (HW, 'haplogroups_work')):
        C.write_tsv_atomic(d, f'{C.CW}/{fn}.tsv')
    H = pd.DataFrame(h1 + h2); C.write_tsv_atomic(H, f'{C.CW}/id_history.tsv')
    print(H.groupby(['type', 'status']).size().to_string()); print('merged', H.groupby('type').merged.sum().to_dict(), 'split', H.groupby('type').split.sum().to_dict())

if __name__ == '__main__':
    main()
