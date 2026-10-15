"""Step 2b: home arm of every end.
  named contigs (CM*, CHM13, HG002) : arm from the chromosome name + side (arm_source = contig_name).
  unnamed contigs (JA*/JB*)         : 10 kb probes at 250 kb and 1 Mb from the tip (work/windows/probes; catalogue_01) are mapped
      to the full CHM13 v2.0 genome (minimap2 -x asm5, primary alignments). A probe is informative if MAPQ >= 20, aligned >= 5 kb and
      its tip-frame orientation is consistent with the arm (p arm: telomere->centromere = increasing coordinate). Arm = chromosome
      + (p if the probe lies in the first half of the chromosome else q). Confidence: high = two informative probes agree;
      medium = one informative probe (the other absent/uninformative); conflict = informative probes disagree (higher MAPQ x
      aligned length wins); fallback 'mosaic' = no informative probe -> the CHM13 end with the largest summed alignment score over
      the bins starting >= 50 kb from the telomere boundary (or all bins when the subtelomere is shorter). Accuracy of both rules is measured
      on named tips (1,000 random named tips carry probes; all named tips for the mosaic rule).
  Catalogue arm (catalogue_arm): the arm, or its pooled group for ends whose distal sequence is shared between chromosomes
      (acroP = 13p/14p/15p/21p/22p, PAR1 = Xp/Yp, PAR2 = Xq/Yq; catalogue_common.ARM_GROUPS). Probe agreement and conflicts are
      evaluated at group level. A mosaic fallback inside a pooled group gets confidence 'mosaic_group' (group-level accuracy is
      measured on named ends); other mosaic fallbacks are 'low'; ends with neither an informative probe nor any CHM13-aligned bin
      get home_arm 'unassigned' (confidence 'none'). Ends with 'low', 'conflict' or 'none' confidence do not define alleles.
Outputs: work/catalogue/arms.tsv (end_id, contig_arm, home_arm, arm_source, arm_confidence, probe_arm, probe_n_inf, probe_detail,
         mosaic_arm, mosaic_arm_score_frac), work/catalogue/arm_calibration.tsv
"""
import os, sys, subprocess, glob, pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t')
PR = pd.read_csv(f'{C.WIN}/probe_plan.tsv', sep='\t')
# ---- probe mapping
pfa = f'{C.CW}/probes.fa'; paf = f'{C.CW}/probes_vs_chm13.paf'
if not os.path.exists(paf):
    reg2 = {}
    for r in PR.itertuples(): reg2.setdefault((r.asm, r.region), []).append(f'{r.end_id}@{r.probe_offset}')
    with open(pfa, 'w') as f:
        for fn in sorted(glob.glob(f'{C.WIN}/probes/*.fa')):
            asm = os.path.basename(fn)[:-3]
            for n, s in C.readfa(fn):
                for q in reg2.get((asm, n), []): f.write(f'>{q}\n{s}\n')
    mmi = f'{C.CW}/chm13v2.0.asm5.mmi'
    if not os.path.exists(mmi):
        subprocess.run(f"minimap2 -x asm5 -t {C.THREADS} -d {mmi} {C.inp('ref/chm13v2.0.fa.gz')} 2> {C.CW}/mmi.log", shell=True, check=True)
    subprocess.run(f'minimap2 -x asm5 --secondary=no -t {C.THREADS} {mmi} {pfa} > {paf}.tmp 2> {paf}.log && mv {paf}.tmp {paf}', shell=True, check=True)
clen = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(C.inp('ref/chm13v2.0.fa.gz.fai'))}
side = dict(zip(W.end_id, W.side))
best = {}
for line in open(paf):
    f = line.split('\t'); q = f[0]; aln = int(f[3]) - int(f[2]); mq = int(f[11])
    if q not in best or aln * (mq + 1) > best[q][0]: best[q] = (aln * (mq + 1), f[5], int(f[7]), int(f[8]), f[4], mq, aln)
prow = []
for q, (sc, chrom, ts, te, st, mq, aln) in best.items():
    e, off = q.rsplit('@', 1); L = clen[chrom]; mid = (ts + te) / 2; arm = chrom[3:] + ('p' if mid < L / 2 else 'q')
    tipst = st if side[e] == 'L' else ('+' if st == '-' else '-')
    consistent = (arm.endswith('p') and tipst == '+') or (arm.endswith('q') and tipst == '-')
    dist = mid if arm.endswith('p') else L - mid
    prow.append(dict(end_id=e, offset=int(off), arm=arm, mapq=mq, aln=aln, consistent=consistent, dist_from_end=int(dist),
                     informative=mq >= 20 and aln >= 5000 and consistent, score=mq * aln))
PB = pd.DataFrame(prow)
# ---- mosaic fallback rule (summed AS per CHM13 end over proximal bins)
Bn = pd.read_csv(f'{C.CW}/bins.tsv.gz', sep='\t', dtype={'bin': str, 'cand': str})
Bn = Bn[~Bn.bin.isin(['S5', 'S10', 'S20', 't20'])]
s0 = dict(zip(W.end_id, W.junk_len + W.telo_len))
mos = {}
for e, g in Bn.groupby('end_id'):
    prox = g[g.qstart - s0[e] >= 50000]; prox = prox if len(prox) else g; sc = {}
    for c in prox.cand.dropna():
        for x in c.split(';'):
            t, a, _ = x.split(':'); sc[t] = sc.get(t, 0) + int(a)
    if sc:
        b = max(sorted(sc), key=sc.get); mos[e] = (C.arm_of(b), round(sc[b] / sum(sc.values()), 3))
rows = []; PBg = {e: g for e, g in PB.groupby('end_id')}
for w in W.itertuples():
    e = w.end_id; p = PBg.get(e, PB.iloc[0:0])
    inf = p[p.informative]; parm, conf = '', ''
    if len(inf):
        groups = set(inf.arm.map(C.group_of)); parm = inf.sort_values(['score', 'offset']).arm.iloc[-1]
        conf = ('high' if len(inf) >= 2 else 'medium') if len(groups) == 1 else 'conflict'
    ma, mf = mos.get(e, ('', np.nan))
    if isinstance(w.contig_arm, str) and w.contig_arm:
        home, src, hc = w.contig_arm, 'contig_name', 'named'
    elif parm:
        home, src, hc = parm, 'probe', conf
    else:
        home, src, hc = (ma, 'mosaic', ('mosaic_group' if C.group_of(ma) != ma else 'low')) if ma else ('unassigned', 'none', 'none')
    rows.append(dict(end_id=e, contig_arm=w.contig_arm if isinstance(w.contig_arm, str) else '', home_arm=home,
                     catalogue_arm=C.group_of(home), arm_source=src,
                     arm_confidence=hc, probe_arm=parm, probe_confidence=conf, probe_n_inf=len(inf), probe_n=len(p),
                     probe_detail=';'.join(f'{r.offset}:{r.arm}:{r.mapq}:{r.aln}:{int(r.consistent)}' for r in p.itertuples()),
                     mosaic_arm=ma, mosaic_arm_score_frac=mf))
A = pd.DataFrame(rows)
C.write_tsv_atomic(A, f'{C.CW}/arms.tsv'); PB.to_csv(f'{C.CW}/probes_best.tsv', sep='\t', index=False)
# ---- calibration on named tips (exact arm and catalogue group)
N = A[A.arm_source == 'contig_name'].copy(); cal = []
N['g'] = N.contig_arm.map(C.group_of); N['pg'] = N.probe_arm.map(C.group_of); N['mg'] = N.mosaic_arm.map(C.group_of)
def add(rule, conf, g, col, gcol):
    cal.append(dict(rule=rule, confidence=conf, n=len(g), correct_arm=int((g[col] == g.contig_arm).sum()), correct_group=int((g[gcol] == g.g).sum())))
pc = N[N.probe_confidence != '']
for conf, g in pc.groupby('probe_confidence'): add('probe', conf, g, 'probe_arm', 'pg')
add('probe', 'any', pc, 'probe_arm', 'pg')
cal.append(dict(rule='probe', confidence='none_informative', n=int(((N.probe_n > 0) & (N.probe_confidence == '')).sum()), correct_arm=0, correct_group=0))
m = N[N.mosaic_arm != '']
add('mosaic', 'all_named', m, 'mosaic_arm', 'mg')
add('mosaic', 'named_in_pooled_groups', m[m.g != m.contig_arm], 'mosaic_arm', 'mg')
add('mosaic', 'mosaic_call_in_pooled_group', m[m.mg != m.mosaic_arm], 'mosaic_arm', 'mg')
add('mosaic', 'mosaic_call_other_arms', m[m.mg == m.mosaic_arm], 'mosaic_arm', 'mg')
CA = pd.DataFrame(cal); CA['accuracy_arm'] = (CA.correct_arm / CA.n.clip(lower=1)).round(4); CA['accuracy_group'] = (CA.correct_group / CA.n.clip(lower=1)).round(4)
C.write_tsv_atomic(CA, f'{C.CW}/arm_calibration.tsv')
print(CA.to_string()); print(A.groupby(['arm_source', 'arm_confidence']).size().to_string())
