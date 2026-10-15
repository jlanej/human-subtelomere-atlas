#!/usr/bin/env python3
"""Diagnose every ape hit in which the human carrier junction probe aligned across the switch point (ape_forms.py).

For each (junction cluster, species) with carrier_present == True, re-align the 6 kb probes and report:
  * ape locus of the carrier-probe hit (end, start/end on the telomere-first ape end) and the ape-side position of the
    switch point (distance from the ape chromosome terminus) vs the human switch distance in the representative carrier;
  * per-flank identity of the carrier probe at that ape locus (foreign/donor-like half = probe[0:3000],
    own/recipient half = probe[3000:6000]) and, for comparison, the identity of each half to CHM13 B and CHM13 A;
  * whether the native own-arm probe (CHM13 A) and native donor probe (CHM13 B) align to the SAME ape locus
    (reciprocal overlap of ape coordinates >= 50%), and with which identity;
  * whether the carrier probe aligns contiguously across the switch point to CHM13 A or CHM13 B themselves
    (human_internal_span): if so the 'switch' is a change in paralog similarity inside sequence shared by A and B, not a
    structural breakpoint, and ape presence of contiguous sequence is not evidence for the human junction.
Classification per hit:
  shared_duplicon   : a native probe aligns to the same ape locus spanning the switch with identity within 0.005 of the carrier
                      probe, or the carrier probe spans CHM13 A or B contiguously (human_internal_span);
  paralog           : carrier-specific by identity but the ape end is not orthologous to human A, or identity < 0.97;
  candidate_shared_junction : carrier-specific, orthologous end (ape end ortholog == A) and identity >= 0.97 on both flanks.
Output: results/ape_junction_hits_diagnosis.tsv
"""
import os, sys
import numpy as np, pandas as pd, mappy as mp
sys.path.insert(0, os.path.dirname(__file__))
import origin_lib as OL
import ape_forms as AF

OUT = os.path.join(OL.REPO, 'analyses', 'allele_ages', 'results')


def hits(al, probe):
    return [h for h in al.map(probe)]


def span_hit(hs, center=3000, flank=1000):
    sp = [h for h in hs if h.q_st <= center - flank and h.q_en >= center + flank]
    return max(sp, key=lambda h: h.mlen) if sp else None


def ident_on(seq_ref, q):
    a = mp.Aligner(seq=seq_ref, preset='asm20')
    hs = list(a.map(q))
    if not hs:
        return np.nan, 0
    h = max(hs, key=lambda h: h.mlen)
    return round(h.mlen / h.blen, 4), round((h.q_en - h.q_st) / len(q), 3)


def main():
    P = pd.read_csv(os.path.join(OUT, 'ape_junction_probes.tsv'), sep='\t')
    O = pd.read_csv(os.path.join(OUT, 'ape_end_orthology.tsv'), sep='\t')
    omap = {f'{r.species}|{r.ape_seq}|{r.ape_end}': r.ortholog for r in O.itertuples()}
    T = pd.read_csv(os.path.join(OL.IN, 'results', 'subtel_foreign_tips.tsv'), sep='\t', low_memory=False).set_index('q')
    E = AF.ape_ends(); TS = OL.TipStore()
    H = {}
    rows = []
    for r in P[P.carrier_present].itertuples():
        t = T.loc[r.rep_tip]
        s = TS.get(t.asm, t.contig, t.end); sf = int(t.switch_final)
        probe = s[sf - 3000: sf + 3000]
        A, B = 'chr' + t.own, 'chr' + t.foreign_arm
        for e in (A, B):
            if e not in H:
                H[e] = OL.chm13_end(e)
        jA = int(t.junction)
        natA = H[A][max(0, jA - 3000): jA + 3000]
        jB = int(r.chm13_B_switch_coord) if np.isfinite(r.chm13_B_switch_coord) else None
        natB = H[B][max(0, jB - 3000): jB + 3000] if jB else ''
        ape = E[r.carrier_end]
        al = mp.Aligner(seq=ape, preset='asm20')
        hc = span_hit(hits(al, probe))
        if hc is None:
            continue
        ape_switch = hc.r_st + (3000 - hc.q_st) if hc.strand == 1 else hc.r_en - (3000 - hc.q_st)
        loc = ape[hc.r_st:hc.r_en]
        idF, cvF = ident_on(loc, probe[:3000]); idO, cvO = ident_on(loc, probe[3000:])
        idF_B, _ = ident_on(H[B][:300_000], probe[:3000]); idO_A, _ = ident_on(H[A][:300_000], probe[3000:])
        def same_locus(q):
            if not len(q):
                return False, np.nan
            h = span_hit(hits(al, q))
            if h is None:
                return False, np.nan
            ov = min(h.r_en, hc.r_en) - max(h.r_st, hc.r_st)
            return ov >= 0.5 * min(h.r_en - h.r_st, hc.r_en - hc.r_st), round(h.mlen / h.blen, 4)
        sA, iA = same_locus(natA); sB, iB = same_locus(natB)
        # human-internal: does the carrier probe span the switch contiguously on CHM13 A or B?
        hiA = span_hit(hits(mp.Aligner(seq=H[A][:300_000], preset='asm20'), probe))
        hiB = span_hit(hits(mp.Aligner(seq=H[B][:300_000], preset='asm20'), probe))
        human_internal = (hiA is not None and hiA.mlen / hiA.blen >= 0.9) or (hiB is not None and hiB.mlen / hiB.blen >= 0.9)
        cid = round(hc.mlen / hc.blen, 4)
        orth = omap.get(r.carrier_end, 'none')
        if human_internal or (sA and cid - iA < 0.005) or (sB and cid - iB < 0.005):
            cls = 'shared_duplicon'
        elif orth != A or min(idF, idO) < 0.97:
            cls = 'paralog'
        else:
            cls = 'candidate_shared_junction'
        rows.append(dict(junction_cluster=r.junction_cluster, form=r.form, samples=r.samples, species=r.species,
                         ape_end=r.carrier_end, ape_end_ortholog=orth, human_own_end=A, human_donor_end=B,
                         human_switch_from_telomere=sf, ape_switch_from_telomere=int(ape_switch),
                         carrier_probe_ident=cid, flank_donor_ident_ape=idF, flank_own_ident_ape=idO,
                         flank_donor_ident_chm13B=idF_B, flank_own_ident_chm13A=idO_A,
                         nativeA_same_locus=sA, nativeA_ident=iA, nativeB_same_locus=sB, nativeB_ident=iB,
                         carrier_probe_spans_chm13A=hiA is not None, carrier_probe_spans_chm13B=hiB is not None,
                         chm13A_span_ident=round(hiA.mlen / hiA.blen, 4) if hiA else np.nan,
                         chm13B_span_ident=round(hiB.mlen / hiB.blen, 4) if hiB else np.nan,
                         prior_call=r.call, diagnosis=cls))
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(OUT, 'ape_junction_hits_diagnosis.tsv'), sep='\t', index=False)
    print(D.diagnosis.value_counts().to_dict())
    print('clusters with any carrier-probe hit', D.junction_cluster.nunique())
    print(D[['junction_cluster', 'species', 'ape_end_ortholog', 'human_switch_from_telomere', 'ape_switch_from_telomere',
             'carrier_probe_ident', 'flank_donor_ident_ape', 'flank_own_ident_ape', 'nativeA_same_locus', 'nativeA_ident',
             'carrier_probe_spans_chm13A', 'chm13A_span_ident', 'carrier_probe_spans_chm13B', 'chm13B_span_ident', 'diagnosis']].to_string())


if __name__ == '__main__':
    main()
