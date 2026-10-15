"""Step 2e: sequence haplogroups per arm from all-vs-all comparison of end windows.
Sequence compared: the first HG_LEN (80 kb) of the subtelomere (window minus masked telomere/junk prefix s0), so that all full
  ends are compared over the same telomere-anchored span; ends with < HG_LEN of subtelomere are 'partial'.
Similarity: Mash-style ANI from FracMinHash sketches (catalogue_sketch.py; k = 21, scale 1/10, ~8,000 hashes per end) computed for
  all pairs of full ends of an arm by a sparse matrix product (exact for the sketches). k-mer ANI integrates base differences
  and SV/indel content. It is calibrated against alignment identity on CALIB_PAIRS random within-arm pairs (mappy asm20:
  gap-compressed identity of the aligned bases and aligned fraction).
Clustering: average linkage (UPGMA) on 1 - ANI, cut at 1 - HG_ANI (threshold chosen from the ANI distribution, see QC.md).
  Haplogroup ID HG-<arm>-NN (NN by size, ties by medoid end_id). Medoid = member with the highest mean ANI to the other members.
Clustered (defining) ends: full-length HPRC_R2/CHM13/HG002_Q100 ends of the catalogue arm with a reliable arm call
  (catalogue_common.DEFINING_CONFIDENCE). Other ends - partial windows, pedigree_aux ends, ends with a low-confidence or
  conflicting arm call - are compared over their available subtelomere (up to HG_LEN) with the same-length prefix of every
  clustered end of the arm and assigned to the haplogroup of the nearest one if that ANI >= HG_ANI (assignment 'partial',
  'aux' or 'arm_lowconf'), else 'unassigned_<kind>'.
Validation: every member aligned to its medoid (mappy asm20): identity and aligned fraction (columns medoid_aln_*).
Outputs: work/catalogue/haplogroup_members.tsv, haplogroups_work.tsv, hg_pair_ani_hist.tsv, hg_nn_ani.tsv, hg_calibration.tsv
Usage  : python catalogue_07_haplogroups.py [HG_ANI]
"""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, multiprocessing as mp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C, catalogue_sketch as S
HG_LEN, K, SCALE = 80_000, 21, 10
HG_ANI = float(sys.argv[1]) if len(sys.argv) > 1 else float(os.environ.get('ATLAS_HG_ANI', '0.995'))
CALIB_PAIRS, SEED = 2000, 20260715

def jacc_matrix(sk):
    allh = np.unique(np.concatenate(sk)) if sk else np.zeros(0, np.uint64)
    rows = np.concatenate([np.full(len(s), i) for i, s in enumerate(sk)]); cols = np.concatenate([np.searchsorted(allh, s) for s in sk])
    X = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(sk), len(allh)))
    I = (X @ X.T).toarray().astype(np.float64); sz = np.array([len(s) for s in sk], float)
    return I / np.maximum(sz[:, None] + sz[None, :] - I, 1)

_SEQ = {}
def gc_ident(args):
    """gap-compressed identity (1 - de) over primary hits, weighted by aligned query length, and aligned fraction"""
    import mappy, re
    q, t = args; a = mappy.Aligner(seq=_SEQ[t], preset='asm20'); qs = _SEQ[q]
    cov = np.zeros(len(qs), bool); w = d = 0.0
    for h in a.map(qs, cs=True):
        if not h.is_primary: continue
        cov[h.q_st:h.q_en] = True
        sub = len(re.findall(r'\*', h.cs)); ins = len(re.findall(r'\+', h.cs)); dl = len(re.findall(r'-', h.cs))
        mt = sum(int(x) for x in re.findall(r':(\d+)', h.cs)); L = mt + sub + ins + dl
        w += L; d += sub + ins + dl
    return q, t, (1 - d / w if w else np.nan), float(cov.mean())

def main():
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t').set_index('end_id')
    A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
    sub = {}
    for n, s in C.readfa(f'{C.WIN}/ends.fa'):
        s0 = int(W.at[n, 'junk_len'] + W.at[n, 'telo_len']); sub[n] = s[s0:s0 + HG_LEN]
    _SEQ.update(sub)
    full = {e for e, s in sub.items() if len(s) >= HG_LEN and W.at[e, 'set'] in ('HPRC_R2', 'CHM13', 'HG002_Q100')
            and A.at[e, 'arm_confidence'] in C.DEFINING_CONFIDENCE}
    kind_of = lambda e: 'aux' if W.at[e, 'set'] == 'pedigree_aux' else ('partial' if len(sub[e]) < HG_LEN else 'arm_lowconf')
    rng = np.random.default_rng(SEED)
    mem, hgs, hist, nn, calib_pool, scan = [], [], [], [], [], []
    GRID = [0.98, 0.985, 0.99, 0.9925, 0.995, 0.996, 0.997, 0.998, 0.999]
    bins = np.round(np.arange(0.90, 1.00001, 0.0005), 4)
    for arm, g in sorted(A.groupby('catalogue_arm')):
        ids = sorted(e for e in g.index if e in full); part = sorted(e for e in g.index if e not in full)
        if not ids: continue
        sk = [S.sketch(sub[e], K, SCALE) for e in ids]
        ani = S.ani_from_jaccard(jacc_matrix(sk), K); np.fill_diagonal(ani, 1.0)
        n = len(ids)
        if n >= 2:
            iu = np.triu_indices(n, 1); v = ani[iu]
            h, _ = np.histogram(np.clip(v, 0.90, 1.0), bins=np.append(bins, 1.00051)); hist.append(pd.Series(h, index=bins, name=arm))
            Z = linkage(squareform(1 - ani, checks=False), 'average'); lab = fcluster(Z, t=1 - HG_ANI, criterion='distance')
            for th in GRID:
                lb = fcluster(Z, t=1 - th, criterion='distance'); cnt = np.bincount(lb)[1:]
                scan.append(dict(arm=arm, threshold=th, n_ends=n, n_groups=len(cnt), n_singletons=int((cnt == 1).sum()), largest=int(cnt.max())))
            for i, e in enumerate(ids):
                o = np.delete(ani[i], i); nn.append(dict(end_id=e, arm=arm, nn_ani=float(o.max())))
            for _ in range(min(CALIB_PAIRS // 46 + 1, len(v))):
                a, b = rng.choice(n, 2, replace=False); calib_pool.append((ids[a], ids[b], float(ani[a, b]), arm))
        else:
            lab = np.array([1])
        cl = {}
        for i, l in enumerate(lab): cl.setdefault(l, []).append(i)
        groups = []
        for l, ix in cl.items():
            sub_ani = ani[np.ix_(ix, ix)]; meanani = (sub_ani.sum(1) - 1) / max(len(ix) - 1, 1)
            order = sorted(range(len(ix)), key=lambda j: (-meanani[j], ids[ix[j]])); med = ix[order[0]]
            groups.append((len(ix), ids[med], ix, med, meanani))
        groups.sort(key=lambda x: (-x[0], x[1]))
        hid_of = {}
        for r, (sz, medid, ix, med, meanani) in enumerate(groups, 1):
            hid = f'HG-{arm}-{r:02d}'
            off = [j for j in range(n) if j not in set(ix)]
            nearest_other = float(ani[med, off].max()) if off else np.nan
            within = ani[np.ix_(ix, ix)][np.triu_indices(len(ix), 1)] if len(ix) > 1 else np.array([])
            hgs.append(dict(haplogroup_id=hid, arm=arm, n_full=sz, medoid_end_id=medid,
                            min_within_ani=round(float(within.min()), 5) if len(within) else np.nan,
                            mean_within_ani=round(float(within.mean()), 5) if len(within) else np.nan,
                            medoid_nearest_other_member_ani=round(nearest_other, 5)))
            for jj, j in enumerate(ix):
                hid_of[ids[j]] = hid
                mem.append(dict(end_id=ids[j], arm=arm, haplogroup_id=hid, assignment='full',
                                ani_to_medoid=round(float(ani[j, med]), 5), mean_ani_to_members=round(float(meanani[jj]), 5),
                                nn_ani=round(float(np.delete(ani[j], j).max()), 5) if n > 1 else np.nan, subtel_len=len(sub[ids[j]])))
        # partial ends: nearest full end over the same prefix length
        for e in part:
            Lp = len(sub[e])
            if Lp < 5000 or not ids:
                mem.append(dict(end_id=e, arm=arm, haplogroup_id='', assignment=f'unassigned_{kind_of(e)}', subtel_len=Lp)); continue
            ske = S.sketch(sub[e], K, SCALE); best = (-1, None)
            for f_ in ids:
                skf = S.sketch(sub[f_][:Lp], K, SCALE); inter = len(np.intersect1d(ske, skf, assume_unique=True))
                a_ = float(S.ani_from_jaccard(inter / max(len(ske) + len(skf) - inter, 1), K))
                if a_ > best[0] or (a_ == best[0] and f_ < best[1]): best = (a_, f_)
            ok = best[0] >= HG_ANI
            mem.append(dict(end_id=e, arm=arm, haplogroup_id=hid_of[best[1]] if ok else '', assignment=kind_of(e) if ok else f'unassigned_{kind_of(e)}',
                            nn_ani=round(best[0], 5), nn_end_id=best[1], subtel_len=Lp))
        print(arm, 'full', n, 'partial', len(part), 'haplogroups', len(groups), flush=True)
    Mm = pd.DataFrame(mem); H = pd.DataFrame(hgs)
    med = dict(zip(H.haplogroup_id, H.medoid_end_id))
    jobs = [(r.end_id, med[r.haplogroup_id]) for r in Mm.itertuples() if r.haplogroup_id and r.end_id != med[r.haplogroup_id]]
    cj = [(a, b) for a, b, _, _ in calib_pool]
    with mp.get_context('fork').Pool(C.THREADS) as pool:
        va = pool.map(gc_ident, jobs, chunksize=20); ca = pool.map(gc_ident, cj, chunksize=20)
    V = pd.DataFrame(va, columns=['end_id', 'medoid', 'medoid_aln_ident', 'medoid_aln_frac'])
    Mm = Mm.merge(V[['end_id', 'medoid_aln_ident', 'medoid_aln_frac']], on='end_id', how='left')
    Mm.loc[Mm.end_id.isin(set(H.medoid_end_id)), ['medoid_aln_ident', 'medoid_aln_frac']] = [1.0, 1.0]
    CB = pd.DataFrame([dict(q=a, t=b, ani=ani_, arm=arm, aln_ident=r[2], aln_frac=r[3]) for (a, b, ani_, arm), r in zip(calib_pool, ca)])
    st = Mm[Mm.assignment != 'aux'].groupby('haplogroup_id').agg(n_members=('end_id', 'size'), n_partial=('assignment', lambda x: (x == 'partial').sum()),
                                         min_medoid_aln_ident=('medoid_aln_ident', 'min'), median_medoid_aln_ident=('medoid_aln_ident', 'median'),
                                         min_medoid_aln_frac=('medoid_aln_frac', 'min')).reset_index()
    H = H.merge(st, on='haplogroup_id', how='left'); H['threshold_ani'] = HG_ANI; H['compared_len'] = HG_LEN
    C.write_tsv_atomic(Mm, f'{C.CW}/haplogroup_members.tsv'); C.write_tsv_atomic(H, f'{C.CW}/haplogroups_work.tsv')
    C.write_tsv_atomic(pd.concat(hist, axis=1).fillna(0).astype(int).rename_axis('ani_bin').reset_index(), f'{C.CW}/hg_pair_ani_hist.tsv')
    C.write_tsv_atomic(pd.DataFrame(nn), f'{C.CW}/hg_nn_ani.tsv'); C.write_tsv_atomic(CB, f'{C.CW}/hg_calibration.tsv')
    C.write_tsv_atomic(pd.DataFrame(scan), f'{C.CW}/hg_threshold_scan.tsv')
    print('haplogroups', len(H), 'members', len(Mm), Mm.assignment.value_counts().to_dict())

if __name__ == '__main__':
    main()
