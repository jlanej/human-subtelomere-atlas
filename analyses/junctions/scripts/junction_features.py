"""H2 step 2 -- sequence features at base-resolved foreign-block junctions, against matched random positions.

Junction set: Release 2 tips with switch_res == 'base' (results/subtel_foreign_tips.tsv); ONE representative per
junction_cluster (same form, junction coordinates within 1 kb): the member closest to the cluster median junction
(ties: higher min(sw_for_ident_5kb, sw_own_ident_5kb), then tip name).
For each representative the subtelomere (tip frame, telomere masked) is re-aligned (mappy asm20, hits >= 2 kb) to the donor
CHM13 end used for its base-level switch (form_foreign_end) and to the recipient (own) end; per-base error tracks give the
changepoint x (recomputed; must equal the published switch) and the parental coordinates: donor telomere-frame position a-1 of
the last carrier base before x, recipient position b of the carrier base at x.  The carrier is modelled as D_tf[0:a] + R_tf[b:].
Features
  mh_parental      bases of identical sequence shared by the two parental sequences at the junction: forward
                   (D_tf[a+i] == R_tf[b+i]) + backward (D_tf[a-1-i] == R_tf[b-1-i]), capped at 5,000; = breakpoint ambiguity.
                   0 bp = blunt, 1-25 bp = microhomology, > 25 bp = homology-length identity
  bp_interval      carrier-track breakpoint interval: bases between the last donor-specific site before x (carrier = donor !=
                   recipient) and the first recipient-specific site at/after x
  dr_ident_pw      donor-vs-recipient parental identity over +-500 bp (edlib global alignment of D_tf[a-500:a+500] vs
                   R_tf[b-500:b+500]; matches / alignment columns)
  dr_score_trk     in the carrier +-500 bp window: fraction of bases aligned to BOTH parents with both matching the carrier
                   (local donor-recipient sequence identity seen through the carrier; 0 if the parents are not homologous there)
  dr_span_both     True if donor alignment continues >= 500 bp past x AND recipient alignment starts >= 500 bp before x
  its_*            telomere-like repeats within +-1 kb: exact TTAGGG/CCCTAA hexamers (its_exact), degenerate hexamers
                   (<= 1 mismatch, its_deg), starts of in-phase tandem runs >= 3 degenerate copies (its_run) -- in the
                   carrier window and in each parental window
Null models (matched random positions; NPERM draws per junction, seed fixed)
  parental pairs   random (a', b') uniform in telomere-frame [5 kb, 200 kb) of the SAME donor and recipient CHM13 ends:
                   mh_parental, dr_ident_pw (first 200 draws), parental ITS
  carrier          random x' uniform in the SAME carrier subtelomere [5 kb, n-1 kb) excluding x +- 2 kb: dr_score_trk, carrier ITS
                   (5 kb lower bound keeps null windows clear of the terminal telomere array, as for the parental null)
Junction clusters are defined per form, so one physical junction can appear under two forms when its donor block is equally
close to paralogous CHM13 ends (e.g. 4q>15p#1 and 4q>4p#0 share recipient coordinate 80,572); 'dedup' summaries keep one
representative per recipient arm and recipient coordinate (single linkage <= 1 kb on junction_R_tf; largest cluster kept).
ITS tests are reported for all junctions and for the subset >= 5 kb from the telomere in carrier and both parents
('far' subset), because subtelomere-proximal (TTAGGG)n-like tracts make junctions < 5 kb from the telomere trivially ITS-rich.
junction_class: homology_spanning (donor alignment continues >= 500 bp past x and recipient alignment starts >= 500 bp before
x, i.e. the two parental ends are collinear-homologous across the switch); otherwise by mh_parental: blunt (0 bp),
microhomology (1-25 bp), short_homology (> 25 bp).
Statistics: per junction empirical p (fraction of draws >= observed); aggregate permutation test of the mean over junctions
(observed mean vs the means of draw k across junctions, k = 1..NPERM; p = (1 + #null >= obs)/(NPERM + 1)).
Own terminal sequence replaced: distribution of sw_own_dist over all base-resolved tips and over cluster representatives,
window-resolved tips separately; 'through' tips give lower bounds (subtelomere length) only.
Inputs : results/subtel_foreign_tips.tsv, subtel/foreign/r2/tips_input.tsv, subtel/foreign/ref/chm13_ends500k.fa
Outputs: results/junction_features.tsv (per representative), results/junction_null_summary.tsv, results/junction_windows.fa
         (carrier/donor/recipient +-500 bp), results/junction_lost_sequence.tsv, results/junction_summary.json
Usage  : python scripts/junction_features.py
"""
import os, sys, json, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J
import edlib

NPERM = 2000; NPW = 200; MHCAP = 5000; HW = 500; IW = 1000; PLO, PHI = 5000, 200000

def mh(D, R, a, b, cap=MHCAP):
    f = 0
    while f < cap and a + f < len(D) and b + f < len(R) and D[a + f] == R[b + f] and D[a + f] != 'N': f += 1
    k = 0
    while k < cap and a - 1 - k >= 0 and b - 1 - k >= 0 and D[a - 1 - k] == R[b - 1 - k] and D[a - 1 - k] != 'N': k += 1
    return f + k, f, k

def pw_ident(s, t):
    if not s or not t: return np.nan
    r = edlib.align(s, t, mode='NW', task='path')
    nice = edlib.getNiceAlignment(r, s, t)
    m = nice['matched_aligned']; return round(m.count('|') / len(m), 4) if m else np.nan

def its_counts(tr, c, lo, hi):
    ex, dg, run = tr
    return J.wsum(c[0], lo, hi), J.wsum(c[1], lo, hi), J.wsum(c[2], lo, hi)

def cum(tr): return [np.concatenate([[0], np.cumsum(t)]) for t in tr]

def tf_pos(end, gpos): return int(gpos) if end.endswith('p') else J.EL - 1 - int(gpos)

G = {}
def init():
    ends = J.chm13_ends()
    G['tf'] = {e: J.chm13_teloframe(e) for e in ends}
    G['its'] = {}
def chm_its(e):
    if e not in G['its']:
        tr = J.its_tracks(G['tf'][e][:PHI + IW + 10]); G['its'][e] = (tr, cum(tr))
    return G['its'][e]

def one(rec):
    rng = np.random.default_rng(J.SEED + int(rec['idx']))
    sub = rec['seq'][rec['s0']:]; n = len(sub)
    De, Re = rec['dend'], 'chr' + rec['own']
    ef, alF, HF = J.tracks(De, sub); eo, alO, HO = J.tracks(Re, sub)
    out = dict(q=rec['q'])
    if not HF or not HO: out['status'] = 'no_alignment'; return out, None
    x, c = J.changepoint(ef, eo)
    out.update(x_recomputed=x, x_published=rec['sw_q'], x_match=bool(x == rec['sw_q']))
    lo, hi = max(0, x - 2000), min(n, x + 2000)
    fi = int(np.argmin([E[0][lo:x].sum() if x > lo else 1e9 for H, E in HF])); oi = int(np.argmin([E[0][x:hi].sum() if hi > x else 1e9 for H, E in HO]))
    tD = HF[fi][1][1]; tR = HO[oi][1][1]
    if x < 1 or tD[x - 1] < 0 or tR[x] < 0: out['status'] = 'unmapped_junction_base'; return out, None
    a = tf_pos(De, tD[x - 1]) + 1; b = tf_pos(Re, tR[x])
    # orientation check: carrier must be co-linear with telomere-frame parents
    out['strand_D'] = HF[fi][0][2]; out['strand_R'] = HO[oi][0][2]
    out['colinear'] = bool((out['strand_D'] == 1) == De.endswith('p') and (out['strand_R'] == 1) == Re.endswith('p'))
    Dt, Rt = G['tf'][De], G['tf'][Re]
    out.update(status='ok', junction_D_tf=a, junction_R_tf=b, dist_D=rec['sw_for_dist'], dist_R=rec['sw_own_dist'])
    m, mf, mb = mh(Dt, Rt, a, b); out.update(mh_parental=m, mh_fwd=mf, mh_bwd=mb)
    # carrier-track breakpoint interval
    dspec = np.where((ef == 0) & (eo == 1))[0]; rspec = np.where((eo == 0) & (ef == 1))[0]
    Ld = dspec[dspec < x].max() if (dspec < x).any() else -1; Fr = rspec[rspec >= x].min() if (rspec >= x).any() else n
    out.update(bp_interval=int(Fr - Ld - 1), bp_last_D_site=int(Ld), bp_first_R_site=int(Fr))
    w0, w1 = max(0, x - 1000), min(n, x + 1000)
    out['mosaic_D_sites_after_x_1kb'] = int(((dspec >= x) & (dspec < w1)).sum()); out['mosaic_R_sites_before_x_1kb'] = int(((rspec < x) & (rspec >= w0)).sum())
    both = alF & alO
    def drs(xx):
        l, h = max(0, xx - HW), min(n, xx + HW + 1)
        return float(((ef[l:h] == 0) & (eo[l:h] == 0) & both[l:h]).mean()) if h > l else np.nan
    out['dr_score_trk'] = drs(x)
    out['dr_span_both'] = bool(alF[x:min(n, x + HW)].all() and alO[max(0, x - HW):x].all())
    out['D_aln_beyond_x'] = int(np.argmin(alF[x:]) if not alF[x:].all() else n - x)
    out['R_aln_before_x'] = int(np.argmin(alO[:x][::-1]) if not alO[:x].all() else x)
    out['dr_ident_pw'] = pw_ident(Dt[max(0, a - HW):a + HW], Rt[max(0, b - HW):b + HW])
    # ITS: carrier and parents
    ctr = J.its_tracks(sub); cc = cum(ctr)
    out['its_exact_c'], out['its_deg_c'], out['its_run_c'] = its_counts(ctr, cc, x - IW, x + IW)
    (dtr, dc), (rtr, rcc) = chm_its(De), chm_its(Re)
    out['its_exact_D'], out['its_deg_D'], out['its_run_D'] = its_counts(dtr, dc, a - IW, a + IW)
    out['its_exact_R'], out['its_deg_R'], out['its_run_R'] = its_counts(rtr, rcc, b - IW, b + IW)
    # nulls
    pa = rng.integers(PLO, PHI, NPERM); pb = rng.integers(PLO, PHI, NPERM)
    cand = np.arange(5000, max(5001, n - 1000)); cand = cand[np.abs(cand - x) >= 2000]
    px = rng.choice(cand, NPERM) if len(cand) else np.full(NPERM, -1)
    N = dict(mh_parental=np.array([mh(Dt, Rt, int(u), int(v))[0] for u, v in zip(pa, pb)]),
             dr_ident_pw=np.array([pw_ident(Dt[u - HW:u + HW], Rt[v - HW:v + HW]) if k < NPW else np.nan for k, (u, v) in enumerate(zip(pa, pb))]),
             dr_score_trk=np.array([drs(int(u)) if u >= 0 else np.nan for u in px]),
             its_exact_c=np.array([its_counts(ctr, cc, u - IW, u + IW)[0] if u >= 0 else np.nan for u in px]),
             its_run_c=np.array([its_counts(ctr, cc, u - IW, u + IW)[2] if u >= 0 else np.nan for u in px]),
             its_exact_P=np.array([its_counts(dtr, dc, u - IW, u + IW)[0] + its_counts(rtr, rcc, v - IW, v + IW)[0] for u, v in zip(pa, pb)]),
             its_run_P=np.array([its_counts(dtr, dc, u - IW, u + IW)[2] + its_counts(rtr, rcc, v - IW, v + IW)[2] for u, v in zip(pa, pb)]))
    out['its_exact_P'] = out['its_exact_D'] + out['its_exact_R']; out['its_run_P'] = out['its_run_D'] + out['its_run_R']
    for k, v in N.items():
        vv = v[~np.isnan(v)] if v.dtype.kind == 'f' else v
        out[f'p_{k}'] = round(float((1 + (vv >= out[k]).sum()) / (1 + len(vv))), 5) if len(vv) else np.nan
        out[f'null_median_{k}'] = float(np.median(vv)) if len(vv) else np.nan
    win = dict(c=sub[max(0, x - HW):x + HW], D=Dt[max(0, a - HW):a + HW], R=Rt[max(0, b - HW):b + HW])
    return out, (N, win)

def main():
    T = pd.read_csv(J.inp('results/subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    B = T[(T.t20_state == 'foreign') & (T.switch_res == 'base') & T.junction_cluster.notna()].copy()
    B['minid'] = B[['sw_for_ident_5kb', 'sw_own_ident_5kb']].min(axis=1)
    reps = []
    for cid, g in B.groupby('junction_cluster'):
        med = g.junction.median(); g = g.assign(dev=(g.junction - med).abs()).sort_values(['dev', 'minid', 'q'], ascending=[True, False, True])
        r = g.iloc[0]; reps.append(dict(junction_cluster=cid, q=r.q, cluster_tips=len(g), cluster_samples=int(r.junction_cluster_samples), junction_median=med))
    Rp = pd.DataFrame(reps).merge(T[['q', 'asm', 'sample', 'own', 'form', 'foreign_arm', 'form_foreign_end', 's0', 'switch_final', 'sw_own_dist', 'sw_for_dist', 'sw_for_ident_5kb', 'sw_own_ident_5kb']], on='q')
    Rp = Rp.sort_values('junction_cluster').reset_index(drop=True)
    if os.environ.get('JUNC_DEBUG_MAX_REP'): Rp = Rp.iloc[:int(os.environ['JUNC_DEBUG_MAX_REP'])]  # debugging only; never set by run.sh
    seqs = J.load_tip_seqs(Rp.q)
    jobs = [dict(idx=i, q=r.q, seq=seqs[r.q]['seq'], s0=seqs[r.q]['s0'], dend=r.form_foreign_end, own=r.own, sw_q=int(r.switch_final - r.s0),
                 sw_own_dist=r.sw_own_dist, sw_for_dist=r.sw_for_dist) for i, r in Rp.iterrows()]
    with Pool(J.THREADS, initializer=init) as pool: res = pool.map(one, jobs, chunksize=4)
    O = pd.DataFrame([r[0] for r in res]); X = Rp.merge(O, on='q', how='left')
    X.to_csv(os.path.join(J.RES, 'junction_features.tsv'), sep='\t', index=False)
    with open(os.path.join(J.RES, 'junction_windows.fa'), 'w') as f:
        for (o, extra), r in zip(res, Rp.itertuples()):
            if extra is None: continue
            for k, s in extra[1].items():
                lab = dict(c='carrier', D='donor_' + r.form_foreign_end, R='recipient_chr' + r.own)[k]
                f.write(f'>{r.junction_cluster}|{lab}|{r.q}\n{s}\n')
    def jclass(r):
        if r.dr_span_both: return 'homology_spanning'
        return 'blunt' if r.mh_parental == 0 else ('microhomology' if r.mh_parental <= 25 else 'short_homology')
    X['junction_class'] = [jclass(r) if r.status == 'ok' else None for r in X.itertuples()]
    X['far_from_telomere'] = (X.x_recomputed >= 5000) & (X.junction_D_tf >= 5000) & (X.junction_R_tf >= 5000)
    X.to_csv(os.path.join(J.RES, 'junction_features.tsv'), sep='\t', index=False)
    ok = [i for i, (o, e) in enumerate(res) if e is not None]
    rows = []
    feats = [('all', k) for k in ['mh_parental', 'dr_ident_pw', 'dr_score_trk', 'its_exact_c', 'its_run_c', 'its_run_c_any', 'its_exact_P', 'its_run_P', 'its_run_P_any']] + \
            [('far', k) for k in ['its_exact_c', 'its_run_c', 'its_run_c_any', 'its_exact_P', 'its_run_P', 'its_run_P_any']]
    for i in ok:
        for k in ('its_run_c', 'its_run_P'): res[i][1][0][k + '_any'] = (res[i][1][0][k] > 0).astype(float)
    X['its_run_c_any'] = (X.its_run_c > 0).astype(float); X['its_run_P_any'] = (X.its_run_P > 0).astype(float)
    for sub_, k in feats:
        sel = [i for i in ok if sub_ == 'all' or bool(X.far_from_telomere.iloc[i])]
        Xo = X.iloc[sel]
        obs = Xo[k].astype(float).values
        M = np.vstack([res[i][1][0][k][:NPERM].astype(float) for i in sel])  # junctions x draws
        nd = NPW if k == 'dr_ident_pw' else NPERM; M = M[:, :nd]
        keep = ~np.isnan(obs) & ~np.isnan(M).all(axis=1); obs_m = np.nanmean(obs[keep]); null_m = np.nanmean(M[keep], axis=0)
        p = (1 + (null_m >= obs_m).sum()) / (1 + len(null_m))
        pk = Xo[f'p_{k}'].values[keep] if f'p_{k}' in Xo else np.full(keep.sum(), np.nan)
        rows.append(dict(subset=sub_, feature=k, junctions=int(keep.sum()), obs_mean=round(float(obs_m), 4), obs_median=float(np.nanmedian(obs[keep])),
                         null_mean=round(float(np.mean(null_m)), 4), null_median=float(np.nanmedian(M[keep])), null_draws=nd,
                         perm_p_greater=round(float(p), 5), perm_p_less=round(float((1 + (null_m <= obs_m).sum()) / (1 + len(null_m))), 5),
                         frac_junctions_p_lt_0_05=round(float((pk < 0.05).mean()), 4) if not np.isnan(pk).all() else np.nan,
                         obs_frac_gt0=round(float((obs[keep] > 0).mean()), 4), null_frac_gt0=round(float((M[keep] > 0).mean()), 4)))
    NS = pd.DataFrame(rows); NS.to_csv(os.path.join(J.RES, 'junction_null_summary.tsv'), sep='\t', index=False)
    Xo = X.iloc[ok]
    # de-duplicate physical junctions shared between paralogous-donor forms
    D_ = Xo.sort_values(['own', 'junction_R_tf']).copy(); gid = []; prev = None; k = -1
    for r in D_.itertuples():
        if prev is None or r.own != prev[0] or r.junction_R_tf - prev[1] > 1000: k += 1
        gid.append(k); prev = (r.own, r.junction_R_tf)
    D_['phys_junction'] = gid
    Dd = D_.sort_values(['phys_junction', 'cluster_samples', 'junction_cluster'], ascending=[True, False, True]).drop_duplicates('phys_junction')
    X['phys_junction'] = X.q.map(dict(zip(D_.q, D_.phys_junction))); X['phys_junction_rep'] = X.q.isin(set(Dd.q))
    X.to_csv(os.path.join(J.RES, 'junction_features.tsv'), sep='\t', index=False)
    # mh classes
    def mclass(v): return 'blunt_0' if v == 0 else ('micro_1_4' if v <= 4 else ('micro_5_25' if v <= 25 else ('homology_26_200' if v <= 200 else 'homology_gt200')))
    obs_cls = Xo.mh_parental.map(mclass).value_counts().to_dict()
    nullmh = np.concatenate([res[i][1][0]['mh_parental'] for i in ok]); null_cls = pd.Series(nullmh).map(mclass).value_counts(normalize=True).round(5).to_dict()
    # lost own terminal sequence
    F = T[T.t20_state == 'foreign']
    L = []
    for nm, g, col in [('base_all_tips', F[F.switch_res == 'base'], 'sw_own_dist'), ('base_cluster_reps', X, 'sw_own_dist'),
                       ('window_tips', F[F.switch_res == 'window'], 'sw_own_dist'), ('base_all_tips_donor_block', F[F.switch_res == 'base'], 'sw_for_dist'),
                       ('base_all_tips_net_gain', F[F.switch_res == 'base'].assign(net=lambda d: d.sw_for_dist - d.sw_own_dist), 'net'),
                       ('through_tips_lower_bound', F[F.switch_res == 'through'].assign(lb=lambda d: d.tip_len - d.s0), 'lb')]:
        v = g[col].dropna().astype(float)
        L.append(dict(set=nm, n=len(v), min=v.min(), q10=v.quantile(.1), q25=v.quantile(.25), median=v.median(), q75=v.quantile(.75), q90=v.quantile(.9), max=v.max(),
                      frac_lt_10kb=round((v < 10000).mean(), 4), frac_10_50kb=round(((v >= 10000) & (v < 50000)).mean(), 4), frac_ge_50kb=round((v >= 50000).mean(), 4)))
    LS = pd.DataFrame(L); LS.to_csv(os.path.join(J.RES, 'junction_lost_sequence.tsv'), sep='\t', index=False)
    summ = dict(clusters=len(Rp), analysed=len(ok), status=X.status.value_counts().to_dict(), x_match=int(X.x_match.sum()),
                mh_classes_obs=obs_cls, mh_classes_null_frac=null_cls,
                dr_span_both=int(Xo.dr_span_both.sum()), junction_class=Xo.junction_class.value_counts().to_dict(),
                junction_class_clusters_ge5_samples=Xo[Xo.cluster_samples >= 5].junction_class.value_counts().to_dict(),
                far_from_telomere=int(Xo.far_from_telomere.sum()),
                dedup_physical_junctions=len(Dd), dedup_junction_class=Dd.junction_class.value_counts().to_dict(),
                dedup_mh_classes_obs=Dd.mh_parental.map(mclass).value_counts().to_dict(),
                dedup_mh_gt25_frac=round(float((Dd.mh_parental > 25).mean()), 4),
                its_run_c_far=int(((Xo.its_run_c > 0) & Xo.far_from_telomere).sum()), its_run_P_far=int(((Xo.its_run_P > 0) & Xo.far_from_telomere).sum()), bp_interval_median=float(Xo.bp_interval.median()),
                its_run_c_junctions_with_run=int((Xo.its_run_c > 0).sum()), its_run_P_junctions_with_run=int((Xo.its_run_P > 0).sum()),
                params=dict(NPERM=NPERM, NPW=NPW, MHCAP=MHCAP, half_window_homology=HW, half_window_its=IW, parental_null_range=[PLO, PHI], seed=J.SEED))
    json.dump(summ, open(os.path.join(J.RES, 'junction_summary.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(summ, default=str)); print(NS.to_string()); print(LS.to_string())

if __name__ == '__main__':
    main()
