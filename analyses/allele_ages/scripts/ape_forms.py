#!/usr/bin/env python3
"""Do human foreign subtelomere forms exist in the T2T great-ape genomes?

Inputs : work/allele_ages/apes/<species>.ends500k.fa (fetch_ape_ends.py), CHM13 v2.0 ends (origin_lib.chm13_end),
         $ATLAS_INPUT_DIR/results/subtel_foreign_tips.tsv, 100 kb tip windows (origin_lib.TipStore).
1. Orthology: the interior segment (300-500 kb from the terminus) of every ape end is aligned (mappy asm20) to the 46 CHM13 (chr1-22, X; Y excluded)
   outer-500 kb ends; the CHM13 end with most aligned matching bases (>= 20 kb) is the human ortholog.
2. Terminal profile: 5 kb windows of the ape end (first 40 kb after the telomere) aligned to all 46 CHM13 ends; best end per window.
3. Junction probes for every recurrent base-resolved junction cluster (>= 5 samples): a 6 kb probe centred on the switch of a
   representative carrier (the carrier whose junction is closest to the cluster median; probe = tip[switch-3000, switch+3000]),
   plus a native own-arm probe (CHM13 A at the same own-side coordinate +-3 kb) and a native donor probe (CHM13 B, the 6 kb
   aligned to the carrier probe's foreign half extended). Each probe is aligned to all ape ends (mappy asm20).
   'junction present' = a single alignment covering >= 1 kb on BOTH sides of the switch point (probe position 3000) with
   identity >= 0.90 and >= 90% of the probe; reported with the ape end, its human ortholog and identity.
   call: carrier-specific (carrier probe spans and beats every spanning native probe by >= 0.005 identity) |
   non-diagnostic (a native probe spans with similar identity) | absent.
Outputs: results/ape_end_orthology.tsv, results/ape_junction_probes.tsv, results/ape_forms_summary.tsv
"""
import os, sys, re
import numpy as np, pandas as pd, mappy as mp, pysam
sys.path.insert(0, os.path.dirname(__file__))
import origin_lib as OL

OUT = os.path.join(OL.REPO, 'analyses', 'allele_ages', 'results')
APE = os.path.join(OL.WORK, 'apes')
SPECIES = ['chimpanzee', 'bonobo', 'gorilla', 'sumatran_orangutan', 'bornean_orangutan']
ENDS = [f'chr{c}{a}' for c in list(range(1, 23)) + ['X'] for a in 'pq']


def ape_ends():
    out = {}
    for sp in SPECIES:
        fa = pysam.FastaFile(os.path.join(APE, f'{sp}.ends500k.fa'))
        for n in fa.references:
            s = fa.fetch(n).upper()
            out[n] = s if n.endswith('|L') else OL.rc(s)
    return out


def main():
    E = ape_ends()
    H = {e: OL.chm13_end(e) for e in ENDS}
    # 1. orthology
    with open(os.path.join(OL.WORK, 'chm13_ends_telofirst.fa'), 'w') as f:
        for e in ENDS:
            f.write(f'>{e}\n{H[e]}\n')
    al = mp.Aligner(os.path.join(OL.WORK, 'chm13_ends_telofirst.fa'), preset='asm20')
    orth = []
    for n, s in E.items():
        sc = {}
        for h in al.map(s[300_000:500_000]):
            sc[h.ctg] = sc.get(h.ctg, 0) + h.mlen
        best = max(sc, key=sc.get) if sc else None
        tl = OL.telo_len(s)
        # terminal profile, 5 kb windows from the end of the telomere
        prof = []
        for w in range(8):
            q = s[tl + w * 5000: tl + (w + 1) * 5000]
            hs = list(al.map(q))
            prof.append(max(hs, key=lambda h: h.mlen).ctg if hs else '.')
        sp, chrom, end = n.split('|')
        orth.append(dict(species=sp, ape_seq=chrom, ape_end=end, telo_len=tl, ortholog=best if best and sc[best] >= 20000 else 'none',
                         ortholog_matches=sc.get(best, 0) if best else 0, profile_5kb=','.join(prof),
                         terminal_best=prof[0], terminal_matches_ortholog=prof[0] == best))
    O = pd.DataFrame(orth); O.to_csv(os.path.join(OUT, 'ape_end_orthology.tsv'), sep='\t', index=False)
    omap = {f'{r.species}|{r.ape_seq}|{r.ape_end}': r.ortholog for r in O.itertuples()}
    print('ape ends', len(O), 'with ortholog', int((O.ortholog != 'none').sum()), file=sys.stderr)
    # 3. junction probes
    T = pd.read_csv(os.path.join(OL.IN, 'results', 'subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    J = T[(T.junction_cluster.notna()) & (T.switch_res == 'base')]
    cs = J.groupby('junction_cluster')['sample'].nunique()
    TS = OL.TipStore()
    with open(os.path.join(OL.WORK, 'ape_ends_all.fa'), 'w') as f:
        for n, s in E.items():
            f.write(f'>{n}\n{s}\n')
    aa = mp.Aligner(os.path.join(OL.WORK, 'ape_ends_all.fa'), preset='asm20')
    rows = []

    def best_hits(probe, center):
        out = {}
        for h in aa.map(probe):
            ident = h.mlen / h.blen
            spans = h.q_st <= center - 1000 and h.q_en >= center + 1000
            cov = (h.q_en - h.q_st) / len(probe)
            sp = h.ctg.split('|')[0]
            k = (sp,)
            rec = dict(end=h.ctg, ortholog=omap.get(h.ctg, 'none'), ident=round(ident, 4), cov=round(cov, 3), spans=spans,
                       r_st=h.r_st, r_en=h.r_en, mlen=h.mlen)
            if k not in out or (rec['spans'], rec['mlen']) > (out[k]['spans'], out[k]['mlen']):
                out[k] = rec
        return out

    for jc, n in cs[cs >= 5].sort_values(ascending=False).items():
        g = J[J.junction_cluster == jc]
        med = g.junction.median()
        r = g.iloc[(g.junction - med).abs().argsort().iloc[0]]
        s = TS.get(r.asm, r.contig, r.end)
        if s is None:
            continue
        sf = int(r.switch_final)
        if sf - 3000 < 0 or sf + 3000 > len(s):
            continue
        probe = s[sf - 3000: sf + 3000]
        A, B = 'chr' + r.own, 'chr' + r.foreign_arm
        jA = int(r.junction)
        nat_A = H[A][max(0, jA - 3000): jA + 3000]
        # donor-side native: align carrier foreign half to CHM13 B, take B coordinate of the switch
        hb = [h for h in mp.Aligner(seq=H[B][:300_000], preset='asm20').map(probe[:3000]) if h.strand == 1]
        jB = max(hb, key=lambda h: h.mlen).r_en if hb else None
        nat_B = H[B][max(0, jB - 3000): jB + 3000] if jB else None
        res = {'carrier': best_hits(probe, 3000), 'native_A': best_hits(nat_A, 3000) if len(nat_A) == 6000 else {},
               'native_B': best_hits(nat_B, 3000) if nat_B is not None and len(nat_B) == 6000 else {}}
        for sp in SPECIES:
            d = dict(junction_cluster=jc, form=r.form, samples=int(n), junction_median=med, rep_tip=r.q, species=sp,
                     chm13_B_switch_coord=jB)
            for kind, hh in res.items():
                h = hh.get((sp,))
                ok = bool(h and h['spans'] and h['ident'] >= 0.90 and h['cov'] >= 0.9)
                d[f'{kind}_present'] = ok
                for k in ('end', 'ortholog', 'ident', 'cov'):
                    d[f'{kind}_{k}'] = h[k] if h else None
            rows.append(d)
    P = pd.DataFrame(rows)
    # comparative criterion: the carrier junction is carrier-specific in a species only if the carrier probe spans and its
    # identity exceeds that of every native probe that also spans the switch point by >= MIN_DELTA; if a native probe spans
    # with similar identity, the switch lies inside sequence shared by A and B in that genome and the probe is non-diagnostic.
    MIN_DELTA = 0.005
    def classify(r):
        if not r.carrier_present:
            return 'absent'
        nat = [r[f'{k}_ident'] for k in ('native_A', 'native_B') if r[f'{k}_present']]
        if nat and r.carrier_ident - max(nat) < MIN_DELTA:
            return 'non-diagnostic'
        return 'carrier-specific'
    P['call'] = P.apply(classify, axis=1)
    P.to_csv(os.path.join(OUT, 'ape_junction_probes.tsv'), sep='\t', index=False)
    # summary per junction cluster
    Ssum = P.groupby(['junction_cluster', 'form', 'samples']).agg(
        species_carrier_junction=('carrier_present', 'sum'), species_nativeA_junction=('native_A_present', 'sum'),
        species_nativeB_junction=('native_B_present', 'sum'),
        carrier_species=('carrier_present', lambda x: ','.join(P.loc[x.index][x].species)),
        species_carrier_specific=('call', lambda x: int((x == 'carrier-specific').sum())),
        species_non_diagnostic=('call', lambda x: int((x == 'non-diagnostic').sum())),
        carrier_specific_species=('call', lambda x: ','.join(P.loc[x.index][x == 'carrier-specific'].species))).reset_index().sort_values('samples', ascending=False)
    Ssum.to_csv(os.path.join(OUT, 'ape_forms_summary.tsv'), sep='\t', index=False)
    print(Ssum.to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()
