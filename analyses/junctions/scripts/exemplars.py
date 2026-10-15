"""H2 step 3 -- junction exemplars: base-level carrier vs donor (CHM13) vs recipient (CHM13) alignments, identity tracks, read support.

Exemplars (fixed list, chosen for population recurrence and independent read evidence):
  E1  20q>3q  the most-shared junction cluster of the release (20q>3q#3); carrier = among members whose junction equals the
              cluster's modal coordinate, the one with Verkko junction concordance and most hifiasm layout reads across the
              switch (ties: tip name); read support = hifiasm layout reads across the switch and Verkko concordance
              (columns lay_*, vk_* of results/subtel_foreign_tips.tsv)
  E2  3q>11p  HG002 Q100 v1.1 chr3_MATERNAL q (200 kb tip, subtel/foreign/hg002/hg002_tips.tsv); read support = ONT UL /
              PacBio Revio primary alignments over switch +-5 kb (results/subtel_hg002_junction_reads.tsv)
  E3  16p>9q  HG002 Q100 v1.1 chr16_PATERNAL p (600 kb extension, subtel/foreign/hg002/chr16PAT_600k.input.tsv); same read data
For each: the carrier subtelomere is aligned (mappy asm20, hits >= 2 kb) to the donor end used for its switch and to the
recipient CHM13 end; the changepoint x is recomputed (must equal the published switch); per carrier base the parental base of
the hit that best explains x +- 2 kb on each side (genome-orientation hit, minus-strand bases complemented) is recorded.
Outputs: results/exemplar_<id>_track.tsv (carrier pos relative to x, carrier/donor/recipient base and telomere-frame coordinates,
         match flags, site class, 100 bp rolling identities, over x +- 5 kb), results/exemplar_<id>_alignment.txt (x +- 150 bp,
         three-row alignment), results/exemplars.tsv (coordinates, identities, microhomology, read support)
Usage  : python scripts/exemplars.py
"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J

TW = 5000; AW = 150; RW = 100
COMP = dict(A='T', C='G', G='C', T='A', N='N')

def parental_bases(H, E, n, end):
    """per carrier base: parental base (carrier orientation) and telomere-frame coordinate for one hit"""
    (qs, qe, strand, ts, te, cs), (e, tpos) = H, E
    seq = J.chm13_ends()[end]; b = np.array(['-'] * n, dtype='<U1'); tf = np.full(n, -1)
    ok = tpos >= 0
    for i in np.where(ok)[0]:
        g = tpos[i]; base = seq[g] if g < len(seq) else 'N'; b[i] = base if strand == 1 else COMP.get(base, 'N')
        tf[i] = g if end.endswith('p') else J.EL - 1 - g
    return b, tf, e

def run_one(eid, label, q, seq, s0, dend, own, sw_q):
    sub = seq[s0:]; n = len(sub); Re = 'chr' + own
    ef, alF, HF = J.tracks(dend, sub); eo, alO, HO = J.tracks(Re, sub)
    x, c = J.changepoint(ef, eo)
    lo, hi = max(0, x - 2000), min(n, x + 2000)
    fi = int(np.argmin([E[0][lo:x].sum() if x > lo else 1e9 for H, E in HF])); oi = int(np.argmin([E[0][x:hi].sum() if hi > x else 1e9 for H, E in HO]))
    bD, tfD, eD = parental_bases(HF[fi][0], HF[fi][1], n, dend); bR, tfR, eR = parental_bases(HO[oi][0], HO[oi][1], n, Re)
    a0, a1 = max(0, x - TW), min(n, x + TW)
    rows = []
    for i in range(a0, a1):
        cb = sub[i]; mD = int(bD[i] == cb); mR = int(bR[i] == cb)
        cls = 'shared' if mD and mR else ('donor_specific' if mD else ('recipient_specific' if mR else ('both_unaligned' if bD[i] == '-' and bR[i] == '-' else 'neither')))
        rows.append(dict(rel=i - x, subtel_pos=i, tip_pos=i + s0, carrier=cb, donor=bD[i], donor_tf=int(tfD[i]), recipient=bR[i], recipient_tf=int(tfR[i]),
                         match_donor=mD, match_recipient=mR, site=cls))
    Tk = pd.DataFrame(rows)
    for k in ('donor', 'recipient'):
        Tk[f'ident_{k}_{RW}bp'] = Tk[f'match_{k}'].rolling(RW, center=True, min_periods=RW // 2).mean().round(3)
    Tk.to_csv(os.path.join(J.RES, f'exemplar_{eid}_track.tsv'), sep='\t', index=False)
    # text alignment x +- AW
    w = Tk[(Tk.rel >= -AW) & (Tk.rel < AW)]
    lines = [f'# {eid} {label}', f'# carrier {q} subtelomere positions {int(w.subtel_pos.min())}-{int(w.subtel_pos.max())} (tip frame = +{s0}); '
             f'junction x = {x} (carrier base x is the first recipient-derived base)',
             f'# donor {dend} telomere-frame {int(w[w.donor_tf >= 0].donor_tf.min())}-{int(w[w.donor_tf >= 0].donor_tf.max())}; '
             f'recipient {Re} telomere-frame {int(w[w.recipient_tf >= 0].recipient_tf.min())}-{int(w[w.recipient_tf >= 0].recipient_tf.max())}',
             '# rows: D = donor CHM13 base aligned to the carrier base, C = carrier, R = recipient CHM13 base; "-" = no aligned base',
             '# marks: "|" shared by all three, "d" carrier = donor only, "r" carrier = recipient only, "." neither; "^" = junction', '']
    for k in range(0, len(w), 60):
        b = w.iloc[k:k + 60]
        mk = ''.join('|' if s == 'shared' else 'd' if s == 'donor_specific' else 'r' if s == 'recipient_specific' else '.' for s in b.site)
        jm = ''.join('^' if r == 0 else ' ' for r in b.rel)
        lines += [f'{int(b.rel.iloc[0]):>6} D ' + ''.join(b.donor), f'{"":>6} C ' + ''.join(b.carrier), f'{"":>6} R ' + ''.join(b.recipient),
                  f'{"":>6}   ' + mk, f'{"":>6}   ' + jm, '']
    open(os.path.join(J.RES, f'exemplar_{eid}_alignment.txt'), 'w').write('\n'.join(lines))
    dspec = Tk[(Tk.site == 'donor_specific') & (Tk.rel < 0)].rel; rspec = Tk[(Tk.site == 'recipient_specific') & (Tk.rel >= 0)].rel
    ld = int(dspec.max()) if len(dspec) else None; fr = int(rspec.min()) if len(rspec) else None
    win = lambda col, l, h: round(float(Tk[(Tk.rel >= l) & (Tk.rel < h)][col].mean()), 4)
    from junction_features import mh
    a = int(tfD[x - 1]) + 1; b = int(tfR[x]); m, mf, mb = mh(J.chm13_teloframe(dend), J.chm13_teloframe(Re), a, b)
    span = bool(alF[x:min(n, x + 500)].all() and alO[max(0, x - 500):x].all())
    return dict(exemplar=eid, label=label, mh_parental=m, mh_fwd=mf, mh_bwd=mb, parents_homologous_across_500bp=span, q=q, donor_end=dend, recipient_end=Re, s0=s0, x_recomputed=x, x_published=sw_q, x_match=bool(x == sw_q),
                junction_tip_pos=x + s0, junction_donor_tf=int(tfD[x - 1]) + 1 if tfD[x - 1] >= 0 else None, junction_recipient_tf=int(tfR[x]) if tfR[x] >= 0 else None,
                last_donor_specific_rel=ld, first_recipient_specific_rel=fr, breakpoint_interval=(fr - ld - 1) if ld is not None and fr is not None else None,
                ident_donor_left_5kb=win('match_donor', -TW, 0), ident_recipient_left_5kb=win('match_recipient', -TW, 0),
                ident_donor_right_5kb=win('match_donor', 0, TW), ident_recipient_right_5kb=win('match_recipient', 0, TW),
                donor_specific_sites_left_5kb=int(((Tk.site == 'donor_specific') & (Tk.rel < 0)).sum()),
                recipient_specific_sites_right_5kb=int(((Tk.site == 'recipient_specific') & (Tk.rel >= 0)).sum()),
                discordant_sites=int(((Tk.site == 'donor_specific') & (Tk.rel >= 0)).sum() + ((Tk.site == 'recipient_specific') & (Tk.rel < 0)).sum()))

def main():
    T = pd.read_csv(J.inp('results/subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    JF = pd.read_csv(os.path.join(J.RES, 'junction_features.tsv'), sep='\t')
    out = []
    # E1
    r = JF[JF.junction_cluster == '20q>3q#3'].iloc[0]
    M0 = M = T[T.junction_cluster == '20q>3q#3'].copy(); mode = M.junction.mode().min()
    M = M[M.junction == mode].assign(vc=lambda d: (d.vk_junction_concordant == True).astype(int), lr=lambda d: d.lay_reads_span_1kb.fillna(-1))
    M = M.sort_values(['vc', 'lr', 'q'], ascending=[False, False, True]); qsel = M.q.iloc[0]
    t = T.set_index('q').loc[qsel]; r = r.copy(); r['q'] = qsel
    s = J.load_tip_seqs([r.q])[r.q]
    d = run_one('E1_20q-3q', '20q>3q cluster #3 representative (HPRC Release 2)', r.q, s['seq'], s['s0'], t.form_foreign_end, t.own, int(t.switch_final - t.s0))
    d.update(cluster_tips=int(r.cluster_tips), cluster_samples=int(r.cluster_samples), junction_own_dist=t.sw_own_dist, junction_donor_dist=t.sw_for_dist,
             cluster_modal_junction=float(mode), cluster_members_at_mode=int(len(M)), cluster_members_at_mode_samples=int(M['sample'].nunique()),
             junction_contig_pos_1based=(int(t.switch_final) + 1) if t.end == 'L' else int(t.contig_len - t.switch_final),
             read_support=f"selected tip: hifiasm layout reads spanning switch +-1 kb: {t.lay_reads_span_1kb}; Verkko junction status: {t.vk_junction_status}; "
                          f"cluster: Verkko junction concordant {int((M0.vk_junction_concordant == True).sum())}/{int(M0.vk_junction_status.isin(['concordant', 'discordant', 'outer_differs']).sum())} tested, "
                          f"hifiasm layout tested {int(M0.lay_reads_span_1kb.notna().sum())} (min reads across switch +-1 kb {M0.lay_reads_span_1kb.min()})")
    out.append(d)
    # E2, E3 (HG002)
    HT = pd.read_csv(J.inp('subtel/foreign/hg002/hg002_tips.tsv'), sep='\t').set_index('q')
    HP = pd.read_csv(J.inp('subtel/foreign/hg002/place.tips.tsv'), sep='\t').set_index('q')
    E16 = pd.read_csv(J.inp('subtel/foreign/hg002/chr16PAT_600k.input.tsv'), sep='\t').set_index('q')
    P16 = pd.read_csv(J.inp('subtel/foreign/hg002/chr16PAT_600k.tips.tsv'), sep='\t').set_index('q')
    RD = pd.read_csv(J.inp('results/subtel_hg002_junction_reads.tsv'), sep='\t')
    for eid, q, seqsrc, plc, lab in [('E2_3q-11p', 'HG002|chr3_MATERNAL|q', HT, HP, '3q>11p HG002 Q100 v1.1 chr3_MATERNAL q'),
                                     ('E3_16p-9q', 'HG002|chr16_PATERNAL|p', E16, P16, '16p>9q HG002 Q100 v1.1 chr16_PATERNAL p (600 kb window)')]:
        p = plc.loc[q]
        d = run_one(eid, lab, q, seqsrc.at[q, 'seq'], int(seqsrc.at[q, 's0']), p.form_foreign_end, p.own, int(p.sw_q))
        d.update(junction_own_dist=p.sw_own_dist, junction_donor_dist=p.sw_for_dist)
        rr = RD[RD.q == q]; rs = []
        for ds, g in rr.groupby('dataset'):
            cr = g[g.crosses_10kb == True]
            d[f'{ds}_reads'] = len(g); d[f'{ds}_crossing_10kb'] = len(cr); d[f'{ds}_conflict'] = int((g.conflict == True).sum())
            d[f'{ds}_crossing_ident_median'] = round(float(cr.ident_10kb.median()), 4) if len(cr) else None
            rs.append(f"{ds}: {len(cr)} primary alignments cross switch +-5 kb (median identity {d[f'{ds}_crossing_ident_median']}), {d[f'{ds}_conflict']} clipped at the switch")
        d['read_support'] = '; '.join(rs)
        # contig coordinate of the junction: p end -> tip position; q end -> contig_len - tip position + 1 (1-based)
        L = int(HT.at[q, 'contig_len']) if 'contig_len' in HT.columns else None
        tp = d['junction_tip_pos']; d['junction_contig_pos_1based'] = (tp + 1) if q.endswith('|p') else (L - tp if L else None)
        out.append(d)
    X = pd.DataFrame(out); X.to_csv(os.path.join(J.RES, 'exemplars.tsv'), sep='\t', index=False)
    print(X.T.to_string())

if __name__ == '__main__':
    main()
