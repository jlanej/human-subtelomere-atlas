#!/usr/bin/env python3
"""Single origin and relative age of common foreign subtelomere forms (A>B: own arm A carries a terminal block that
best matches CHM13 arm B).

For each donor CHM13 end B used by a selected form, every Release 2 haplotype segment that should carry B-like terminal
sequence is aligned (minimap2/mappy asm20, telomere-first orientation) to the CHM13 v2.0 B end:
  carriers      : tips with the focal form A>B (t20_state foreign)
  native B      : tips on arm B placed on their own arm (t20_state own)
  other recipients: tips of other forms C>B
A common interval on B, [telomere(B)+1 kb, d1), is defined per form with d1 = 20th percentile of the carriers' contiguous
B-alignment extent (gaps <= 2 kb) minus 500 bp (cap 95 kb). Haplotypes covering < 90% of the interval are dropped.
SNVs vs CHM13 B inside the interval form a genotype matrix; positions within 10 bp of an indel and SNV clusters
(>= 4 SNVs within 50 bp in one haplotype) are masked for that haplotype. Pairwise distance = differing SNV sites /
jointly covered bp.
Statistics per form: pi (mean pairwise SNVs/kb) within carriers, natives, other recipients; dxy; nearest-neighbour label
of each carrier (permutation p, label shuffling, N_PERM); average-linkage (UPGMA) tree: smallest clade containing all /
95% of carriers and its non-carrier content; diagnostic SNVs (alt freq >= 0.9 in carriers and <= 0.1 in natives);
nearest native / other-recipient haplotypes to the carrier consensus; age: rho = mean SNV distance of carriers to the
carrier consensus per bp, T = rho / mu generations (mu = 1.25e-8 /bp/generation, 29 years/generation), and the
mu-independent ratio pi_carrier / pi_native; 95% CIs by bootstrap over 2 kb blocks (N_BOOT).
Inputs : $ATLAS_INPUT_DIR/results/subtel_foreign_tips.tsv, results/subtel_foreign_forms.tsv,
         subtel/foreign/ref/chm13_ends500k.fa, 100 kb tip windows (origin_lib.TipStore); data/metadata/samples.tsv
Outputs: results/origin_forms.tsv (one row per form), results/origin_members.tsv (one row per haplotype x form),
         results/origin_nearest_to_consensus.tsv, work/allele_ages/origin/<form>.npz (distance matrices)
Foreign tips are truncated at their own switch (trunc_of) before alignment so that own-arm sequence beyond the junction,
which can resemble B when A and B share duplicons, is never compared.
Usage  : python origin_age.py [--min-carriers 30] [--forms 20q>3q,3q>11p]
"""
import os, sys, json, argparse, pickle
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(__file__))
import origin_lib as OL

MU, GEN, N_PERM, N_BOOT, SEED = 1.25e-8, 29, 2000, 1000, 20261015
REF_LEN = 140_000
OUT = os.path.join(OL.REPO, 'analyses', 'allele_ages', 'results')
ALN = os.path.join(OL.WORK, 'aln'); NPZ = os.path.join(OL.WORK, 'origin')
for d in (OUT, ALN, NPZ):
    os.makedirs(d, exist_ok=True)

_RA = None; _TS = None


def _init(bend):
    global _RA, _TS
    _RA = OL.RefAligner(OL.chm13_end(bend)[:REF_LEN]); _TS = OL.TipStore()


def _align(rec):
    key, asm, contig, end, trunc = rec
    s = _TS.get(asm, contig, end)
    if s is None:
        return key, None
    if trunc is not None:
        s = s[:trunc]
    if len(s) < 2000:
        return key, ([], {}, [], len(s))
    cov, snv, ind, _ = _RA.align(s)
    return key, (cov, snv, ind, len(s))


def align_end(bend, recs, threads):
    p = os.path.join(ALN, f'{bend}.pkl')
    res = pickle.load(open(p, 'rb')) if os.path.exists(p) else {}
    todo = [r for r in recs if res.get(r[0]) is None]  # retry tips whose window was not yet available
    if todo:
        with Pool(threads, initializer=_init, initargs=(bend,)) as pool:
            for q, v in pool.imap_unordered(_align, todo, chunksize=4):
                res[q] = v
        pickle.dump(res, open(p + '.tmp', 'wb')); os.replace(p + '.tmp', p)
    return res


def trunc_of(r):
    """Query truncation (bp from the contig end) that keeps only the foreign block of a foreign tip: switch_final - 500
    (base-resolved switch) or - 5000 (5 kb window-resolved switch); -1 = no truncation (native tips, 'through' tips)."""
    if r.t20_state != 'foreign' or not np.isfinite(r.switch_final) or r.switch_res == 'through':
        return -1
    return int(r.switch_final - (500 if r.switch_res == 'base' else 5000))


def masked_cov(cov, snv, ind, d0, d1):
    L = d1 - d0
    c = np.zeros(L, bool)
    for a, b in cov:
        a, b = max(a, d0), min(b, d1)
        if b > a:
            c[a - d0:b - d0] = True
    for p, kind, ln in ind:
        a, b = max(p - 10, d0), min(p + max(ln, 1) + 10, d1)
        if b > a:
            c[a - d0:b - d0] = False
    ps = np.array(sorted(x for x in snv if d0 <= x < d1))
    if len(ps) >= 4:
        # clusters: >= 4 SNVs within 50 bp
        j = np.searchsorted(ps, ps + 50)
        for i in np.where(j - np.arange(len(ps)) >= 4)[0]:
            a, b = ps[i] - d0, min(ps[i] + 50, d1) - d0
            c[a:b] = False
    return c


def dist_matrices(G, C, blocks=None):
    """G: n x S int8 allele codes (0 ref, 1-4 alt, -1 uncovered); C: n x L bool coverage.
    Returns diff counts (n x n) and jointly covered bp (n x n)."""
    Cf = C.astype(np.float32)
    joint = Cf @ Cf.T
    cov_s = (G >= 0).astype(np.float32)
    both = cov_s @ cov_s.T
    match = np.zeros_like(both)
    for a in range(5):
        M = (G == a).astype(np.float32); match += M @ M.T
    return both - match, joint


def mean_offdiag(D, idx, jdx=None):
    if jdx is None:
        if len(idx) < 2:
            return np.nan
        sub = D[np.ix_(idx, idx)]; return sub[np.triu_indices(len(idx), 1)].mean()
    if len(idx) == 0 or len(jdx) == 0:
        return np.nan
    return D[np.ix_(idx, jdx)].mean()


def smallest_clade(Z, n, members, frac=1.0):
    """Smallest average-linkage cluster containing >= frac of `members`; returns (size, members_in, nonmembers_in)."""
    root = to_tree(Z)
    target = int(np.ceil(frac * len(members))); mset = set(members); best = None
    stack = [root]
    while stack:
        nd = stack.pop()
        leaves = nd.pre_order()
        k = sum(1 for x in leaves if x in mset)
        if k >= target:
            if best is None or nd.count < best[0]:
                best = (nd.count, k, nd.count - k)
            if not nd.is_leaf():
                stack.extend([nd.left, nd.right])
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-carriers', type=int, default=30)
    ap.add_argument('--min-natives', type=int, default=8)
    ap.add_argument('--forms', default='')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('ATLAS_THREADS', 3)))
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    T = pd.read_csv(os.path.join(OL.IN, 'results', 'subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    S = pd.read_csv(os.path.join(OL.REPO, 'data', 'metadata', 'samples.tsv'), sep='\t')
    T = T.merge(S[['sample', 'superpopulation', 'population']], on='sample', how='left')
    FO = pd.read_csv(os.path.join(OL.IN, 'results', 'subtel_foreign_forms.tsv'), sep='\t')
    forms = a.forms.split(',') if a.forms else FO[FO.n_tips >= a.min_carriers].form.tolist()
    FO = FO.set_index('form')
    rows, mem_rows, near_rows = [], [], []
    for B in sorted({FO.loc[f, 'foreign_arm'] for f in forms}):
        bend = 'chr' + B
        sel = T[((T.t20_state == 'foreign') & (T.foreign_arm == B)) | ((T.own == B) & (T.t20_state == 'own'))]
        sel = sel.assign(trunc=[trunc_of(r) for r in sel.itertuples()])
        sel = sel.assign(key=[f'{q}|{t}' for q, t in zip(sel.q, sel.trunc)])
        recs = [(r.key, r.asm, r.contig, r.end, None if r.trunc < 0 else int(r.trunc)) for r in sel.itertuples()]
        res0 = align_end(bend, recs, a.threads)
        res = {q: res0.get(k) for q, k in zip(sel.q, sel.key)}
        d0 = OL.telo_len(OL.chm13_end(bend)) + 1000
        for f in [x for x in forms if FO.loc[x, 'foreign_arm'] == B]:
            Aarm = FO.loc[f, 'own_arm']
            grp = sel.assign(group=np.where((sel.form == f) & (sel.t20_state == 'foreign'), 'carrier',
                                            np.where(sel.own == B, 'native', 'other_recipient')))
            car = grp[grp.group == 'carrier']
            ext = [OL.contiguous_extent(res[q][0], d0) for q in car.q if res.get(q)]
            row = dict(form=f, own_arm=Aarm, donor_arm=B, carriers_total=len(car), carriers_aligned=len(ext),
                       natives_total=int((grp.group == 'native').sum()), others_total=int((grp.group == 'other_recipient').sum()),
                       interval_start=d0)
            if len(ext) < 10:
                row['status'] = 'too few aligned carriers'; rows.append(row); continue
            d1 = int(min(95_000, np.percentile(ext, 20) - 500)); row['interval_end'] = d1; row['interval_bp'] = d1 - d0
            if d1 - d0 < 3000:
                row['status'] = 'interval < 3 kb'; rows.append(row); continue
            keep, covs, snvs = [], [], []
            for r in grp.itertuples():
                v = res.get(r.q)
                if not v:
                    continue
                c = masked_cov(v[0], v[1], v[2], d0, d1)
                rawcov = masked_cov(v[0], {}, [], d0, d1).mean()
                if rawcov < 0.9:
                    continue
                keep.append(r); covs.append(c); snvs.append({p: b for p, b in v[1].items() if d0 <= p < d1 and c[p - d0]})
            K = pd.DataFrame(keep)
            gi = {g: np.where(K.group.values == g)[0] for g in ('carrier', 'native', 'other_recipient')}
            row.update(carriers_used=len(gi['carrier']), natives_used=len(gi['native']), others_used=len(gi['other_recipient']))
            if len(gi['carrier']) < 10 or len(gi['native']) < a.min_natives:
                row['status'] = 'too few carriers/natives covering interval'; rows.append(row); continue
            sites = sorted({p for s in snvs for p in s}); si = {p: i for i, p in enumerate(sites)}
            code = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
            C = np.vstack(covs); n = len(K)
            G = np.full((n, len(sites)), -1, np.int8)
            sp = np.array(sites) - d0
            for i in range(n):
                G[i, C[i, sp]] = 0
                for p, b in snvs[i].items():
                    G[i, si[p]] = code.get(b, 1)
            diff, joint = dist_matrices(G, C)
            D = diff / np.maximum(joint, 1) * 1000  # SNVs per kb
            np.savez_compressed(os.path.join(NPZ, f.replace('>', '_to_') + '.npz'), D=D, q=K.q.values, group=K.group.values)
            c_, n_, o_ = gi['carrier'], gi['native'], gi['other_recipient']
            pc, pn, po = mean_offdiag(D, c_), mean_offdiag(D, n_), mean_offdiag(D, o_)
            dxy = mean_offdiag(D, c_, n_)
            # nearest neighbours
            Dn = D.copy(); np.fill_diagonal(Dn, np.inf)
            lab = K.group.values
            nn = Dn.argmin(1)
            nn_car = (lab[nn[c_]] == 'carrier').mean()
            perm = []
            for _ in range(N_PERM):
                pl = rng.permutation(lab); ci = np.where(pl == 'carrier')[0]
                perm.append((pl[nn[ci]] == 'carrier').mean())
            p_nn = (1 + np.sum(np.array(perm) >= nn_car)) / (N_PERM + 1)
            # nearest native per carrier
            near_nat = D[np.ix_(c_, n_)].min(1)
            near_car = Dn[np.ix_(c_, c_)].min(1)
            # UPGMA tree
            Z = linkage(squareform(D, checks=False), method='average')
            cl_all = smallest_clade(Z, n, set(c_), 1.0); cl95 = smallest_clade(Z, n, set(c_), 0.95)
            # diagnostic sites
            def altfreq(idx):
                g = G[idx]; cv = (g >= 0).sum(0); return np.where(cv > 0, (g > 0).sum(0) / np.maximum(cv, 1), np.nan)
            fc, fn = altfreq(c_), altfreq(n_)
            diag_c = int(np.sum((fc >= 0.9) & (fn <= 0.1))); diag_n = int(np.sum((fn >= 0.9) & (fc <= 0.1)))
            # carrier consensus and rho
            cons = np.zeros(len(sites), np.int8)
            for j in range(len(sites)):
                g = G[c_, j]; g = g[g >= 0]
                cons[j] = np.bincount(g, minlength=5).argmax() if len(g) else 0
            dcons = np.array([((G[i] >= 0) & (G[i] != cons)).sum() / max(C[i].sum(), 1) for i in range(n)])
            rho = dcons[c_].mean()
            # bootstrap over 2 kb blocks: recompute pi ratio and rho
            L = d1 - d0; nb = max(L // 2000, 1); blk = np.minimum(sp // 2000, nb - 1)
            bc = np.vstack([C[:, b * 2000:min((b + 1) * 2000, L)].sum(1) for b in range(nb)]).T  # n x nb covered bp
            mism = (G >= 0) & (G != cons)
            bdc = np.vstack([mism[:, blk == b].sum(1) for b in range(nb)]).T
            # per block pairwise diffs for carriers & natives (sum over pairs)
            def pair_block(idx):
                out = np.zeros(nb); cov = np.zeros(nb)
                gsub = G[idx]
                for b in range(nb):
                    gb = gsub[:, blk == b]
                    cvb = (gb >= 0).astype(np.float32)
                    both = cvb @ cvb.T; match = sum(((gb == k).astype(np.float32) @ (gb == k).astype(np.float32).T) for k in range(5))
                    cb = C[idx][:, b * 2000:min((b + 1) * 2000, L)].astype(np.float32); jb = cb @ cb.T
                    iu = np.triu_indices(len(idx), 1)
                    out[b] = (both - match)[iu].sum(); cov[b] = jb[iu].sum()
                return out, cov
            pcb, pcc = pair_block(c_); pnb, pnc = pair_block(n_)
            boots_ratio, boots_rho = [], []
            for _ in range(N_BOOT):
                ix = rng.integers(0, nb, nb)
                pcx = pcb[ix].sum() / max(pcc[ix].sum(), 1); pnx = pnb[ix].sum() / max(pnc[ix].sum(), 1)
                boots_ratio.append(pcx / pnx if pnx > 0 else np.nan)
                boots_rho.append((bdc[c_][:, ix].sum(1) / np.maximum(bc[c_][:, ix].sum(1), 1)).mean())
            ratio = pc / pn if pn > 0 else np.nan
            T_gen = rho / MU
            lo_r, hi_r = np.nanpercentile(boots_ratio, [2.5, 97.5]); lo_h, hi_h = np.nanpercentile(boots_rho, [2.5, 97.5])
            # carrier diversity by AFR / non-AFR and superpopulation composition of carriers
            afr = c_[K.superpopulation.values[c_] == 'AFR']; non = c_[K.superpopulation.values[c_] != 'AFR']
            row.update(status='ok', snv_sites=len(sites),
                       pi_carrier=round(pc, 3), pi_native=round(pn, 3), pi_other=round(po, 3) if not np.isnan(po) else np.nan,
                       dxy_carrier_native=round(dxy, 3), da_net=round(dxy - (pc + pn) / 2, 3),
                       dxy_carrier_other=round(mean_offdiag(D, c_, o_), 3) if len(o_) else np.nan,
                       pi_carrier_AFR=round(mean_offdiag(D, afr), 3) if len(afr) > 1 else np.nan,
                       pi_carrier_nonAFR=round(mean_offdiag(D, non), 3) if len(non) > 1 else np.nan,
                       carriers_AFR_used=len(afr),
                       nn_carrier_frac=round(nn_car, 4), nn_perm_expected=round(float(np.mean(perm)), 4), p_nn=p_nn,
                       median_nearest_native=round(float(np.median(near_nat)), 3), median_nearest_carrier=round(float(np.median(near_car)), 3),
                       frac_carriers_nearer_native_than_carrier=round(float(np.mean(near_nat < near_car)), 4),
                       clade_all_size=cl_all[0], clade_all_noncarriers=cl_all[2],
                       clade95_size=cl95[0], clade95_noncarriers=cl95[2],
                       clade95_natives=None, diag_snv_carrier=diag_c, diag_snv_native=diag_n,
                       rho_per_kb=round(rho * 1000, 4), rho_ci95=f'{lo_h * 1000:.4f}-{hi_h * 1000:.4f}',
                       age_gen_mu=round(T_gen), age_kyr_mu=round(T_gen * GEN / 1000, 1),
                       age_kyr_ci95=f'{lo_h / MU * GEN / 1000:.1f}-{hi_h / MU * GEN / 1000:.1f}',
                       pi_ratio_carrier_native=round(ratio, 4), pi_ratio_ci95=f'{lo_r:.4f}-{hi_r:.4f}',
                       native_pairwise_tmrca_kyr_mu=round(pn / 1000 / (2 * MU) * GEN / 1000, 1))
            rows.append(row)
            for i, r in enumerate(K.itertuples()):
                mem_rows.append(dict(form=f, q=r.q, sample=r.sample, superpopulation=r.superpopulation, own_arm=r.own,
                                     group=r.group, member_form=r.form if isinstance(r.form, str) else '', covered_bp=int(C[i].sum()),
                                     snvs_vs_chm13=int((G[i] > 0).sum()), dist_to_carrier_consensus_per_kb=round(dcons[i] * 1000, 4),
                                     nearest_label=lab[nn[i]], nearest_dist=round(Dn[i, nn[i]], 4)))
            order = np.argsort(dcons)
            for i in [j for j in order if lab[j] != 'carrier'][:10]:
                r = K.iloc[i]
                near_rows.append(dict(form=f, q=r.q, sample=r['sample'], superpopulation=r.superpopulation, group=r.group,
                                      member_form=r.form if isinstance(r.form, str) else '', own_arm=r.own,
                                      dist_to_carrier_consensus_per_kb=round(dcons[i] * 1000, 4),
                                      carrier_median_dist_to_consensus=round(float(np.median(dcons[c_])) * 1000, 4)))
            print(f, {k: row[k] for k in ('carriers_used', 'natives_used', 'others_used', 'interval_bp', 'pi_carrier', 'pi_native', 'nn_carrier_frac', 'clade95_noncarriers', 'age_kyr_mu', 'pi_ratio_carrier_native')}, flush=True)
    R = pd.DataFrame(rows).drop(columns=['clade95_natives'], errors='ignore')
    R.to_csv(os.path.join(OUT, 'origin_forms.tsv'), sep='\t', index=False)
    pd.DataFrame(mem_rows).to_csv(os.path.join(OUT, 'origin_members.tsv'), sep='\t', index=False)
    pd.DataFrame(near_rows).to_csv(os.path.join(OUT, 'origin_nearest_to_consensus.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
