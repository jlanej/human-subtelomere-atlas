"""H2 step 1 (whole-genome extension) -- is the recipient's replaced terminal sequence present ELSEWHERE in the same haplotype?

The Release 2 test (scripts/reciprocal.py) can only search chromosome ends (tips and contig ends). HG002 Q100 v1.1 is a complete
diploid assembly, so for every HG002 foreign end with a located switch (switch_res base/window; place tables of
subtel/foreign/hg002/) the CHM13 recipient segment that was replaced (telomere-frame [telomere array end, sw_own_dist)) is
aligned genome-wide to each HG002 haplotype separately and to the full CHM13 v2.0 genome (paralogue background).
minimap2 -x asm20 -c -N 50 -p 0.5 (threads ATLAS_THREADS); alignments >= 1 kb.  Per query and target sequence, aligned
blocks are merged into loci (same contig, same strand, gaps < 20 kb); locus coverage (qcov) = union of query bases covered;
identity = matches / alignment columns (from cg/NM).  Locus position class: 'end:<arm>' if within 500 kb of a contig end
(arm from the contig name), else 'interstitial'.  Self loci excluded: the carrier's own end (same haplotype, same contig end)
and, in CHM13, the recipient end itself.
Calls per query (same haplotype as the carrier):
  present_specific   best non-self locus covers >= 50% at >= 98% identity AND cov*ident exceeds the best non-self CHM13
                     locus (paralogue background) by > 0.05
  paralogue_level    best non-self locus covers >= 50% at >= 98% but is not above the CHM13 paralogue background
  absent             no non-self locus reaching 50% coverage at >= 98% identity
Inputs : work/junctions/hg002v1.1.fasta.gz (scripts/fetch_hg002_genome.sh), ATLAS_INPUT_DIR/ref/chm13v2.0.fa.gz,
         subtel/foreign/hg002/{place.tips.tsv, chr16PAT_600k.tips.tsv}, subtel/foreign/ref/chm13_ends500k.fa
Outputs: results/hg002_lost_sequence_loci.tsv, results/hg002_lost_sequence_calls.tsv; large files in work/junctions/
Usage  : python scripts/lost_sequence_hg002.py
"""
import os, sys, re, subprocess, collections, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J
from reciprocal import chm13_telo_len

END_ZONE = 500_000

def split_haps():
    gz = os.path.join(J.WORK, 'hg002v1.1.fasta.gz')
    outs = {h: os.path.join(J.WORK, f'hg002v1.1.{h}.fa') for h in ('MATERNAL', 'PATERNAL')}
    if all(os.path.exists(o) and os.path.getsize(o) > 1e9 for o in outs.values()): return outs
    fh = {h: open(o + '.tmp', 'w') for h, o in outs.items()}; cur = None
    p = subprocess.Popen(['gzip', '-dc', gz], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        if line[0] == '>':
            n = line[1:].split()[0]; cur = fh['MATERNAL'] if n.endswith('_MATERNAL') else fh['PATERNAL'] if n.endswith('_PATERNAL') else None
        if cur is not None: cur.write(line)
    p.wait()
    for h, f in fh.items(): f.close(); os.replace(outs[h] + '.tmp', outs[h])
    return outs

def chm13_full():
    out = os.path.join(J.WORK, 'chm13v2.0.fa')
    if not (os.path.exists(out) and os.path.getsize(out) > 1e9):
        subprocess.run(f"gzip -dc '{J.inp('ref/chm13v2.0.fa.gz')}' > '{out}.tmp' && mv '{out}.tmp' '{out}'", shell=True, check=True)
    return out

def contig_lens(fa):
    fai = fa + '.fai'
    if not os.path.exists(fai): subprocess.run(['samtools', 'faidx', fa], check=True)
    return {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)}

def arm_of(ctg, pos, L):
    c = ctg.split('_')[0]
    if pos < END_ZONE: return f'end:{c[3:]}p'
    if pos > L - END_ZONE: return f'end:{c[3:]}q'
    return 'interstitial'

def loci(paf, lens):
    A = collections.defaultdict(list)
    for line in open(paf):
        p = line.rstrip('\n').split('\t')
        qs, qe, ql = int(p[2]), int(p[3]), int(p[1])
        if qe - qs < 1000: continue
        nm = int([t for t in p[12:] if t.startswith('NM:i:')][0][5:]); blen = int(p[10])
        A[p[0]].append(dict(ql=ql, qs=qs, qe=qe, st=p[4], t=p[5], ts=int(p[7]), te=int(p[8]), nm=nm, bl=blen))
    rows = []
    for q, al in A.items():
        al.sort(key=lambda a: (a['t'], a['st'], a['ts']))
        groups = []
        for a in al:
            g = groups[-1] if groups else None
            if g and g['t'] == a['t'] and g['st'] == a['st'] and a['ts'] - g['te'] < 20000: g['al'].append(a); g['te'] = max(g['te'], a['te'])
            else: groups.append(dict(t=a['t'], st=a['st'], ts=a['ts'], te=a['te'], al=[a]))
        for g in groups:
            ql = g['al'][0]['ql']; cov = np.zeros(ql, bool); m = c = 0
            for a in g['al']: cov[a['qs']:a['qe']] = True; m += a['bl'] - a['nm']; c += a['bl']
            L = lens[g['t']]; mid = (g['ts'] + g['te']) // 2
            rows.append(dict(query=q, target=g['t'], tstart=g['ts'], tend=g['te'], strand=g['st'], qcov=round(cov.mean(), 4), ident=round(m / c, 5) if c else 0,
                             blocks=len(g['al']), position=arm_of(g['t'], mid, L)))
    return pd.DataFrame(rows)

def mm2(ref, qfa, out):
    if not (os.path.exists(out) and os.path.getsize(out) > 0):
        subprocess.run(f"minimap2 -x asm20 -c -N 50 -p 0.5 -t {J.THREADS} '{ref}' '{qfa}' > '{out}.tmp' 2> '{out}.log' && mv '{out}.tmp' '{out}'", shell=True, check=True)
    return out

def main():
    H = pd.read_csv(J.inp('subtel/foreign/hg002/place.tips.tsv'), sep='\t')
    E16 = pd.read_csv(J.inp('subtel/foreign/hg002/chr16PAT_600k.tips.tsv'), sep='\t')
    H = pd.concat([H[H.q != 'HG002|chr16_PATERNAL|p'], E16], ignore_index=True)
    F = H[(H.t20_state == 'foreign') & H.switch_res.astype(str).str.startswith(('base', 'window')) & H.sw_own_dist.notna()].copy()
    telo = {}; qrows = []
    qfa = os.path.join(J.WORK, 'hg002_lost_segments.fa')
    with open(qfa, 'w') as f:
        for r in F.itertuples():
            R = r.own; t0 = telo.setdefault(R, chm13_telo_len('chr' + R)); Jd = int(r.sw_own_dist)
            seg = J.chm13_teloframe('chr' + R)[t0:Jd]
            if len(seg) < 2000: continue
            name = r.q.replace('|', '_'); f.write(f'>{name}\n{seg}\n')
            qrows.append(dict(query=name, q=r.q, hap=r.q.split('|')[1].split('_')[1], contig=r.q.split('|')[1], end=r.q.split('|')[2], recipient=R,
                              form=f"{R}>{r.t20_best_end[3:]}", switch_res=r.switch_res, sw_own_dist=Jd, lost_tf_start=t0, lost_len=len(seg)))
    Q = pd.DataFrame(qrows)
    haps = split_haps(); chm = chm13_full()
    L = []
    for h, fa in haps.items():
        lens = contig_lens(fa); D = loci(mm2(fa, qfa, os.path.join(J.WORK, f'hg002_lost_vs_{h}.paf')), lens); D['db'] = f'HG002_{h}'; L.append(D)
    lens = contig_lens(chm); D = loci(mm2(chm, qfa, os.path.join(J.WORK, 'hg002_lost_vs_chm13.paf')), lens); D['db'] = 'CHM13'; L.append(D)
    LO = pd.concat(L, ignore_index=True).merge(Q[['query', 'hap', 'contig', 'end', 'recipient']], on='query')
    # self loci
    LO['self'] = ((LO.db == 'HG002_' + LO.hap) & (LO.target == LO.contig) & (LO.position == 'end:' + LO.contig.str.split('_').str[0].str[3:] + LO.end)) | \
                 ((LO.db == 'CHM13') & (LO.position == 'end:' + LO.recipient))
    LO['score'] = LO.qcov * LO.ident
    LO.sort_values(['query', 'db', 'score'], ascending=[True, True, False]).to_csv(os.path.join(J.RES, 'hg002_lost_sequence_loci.tsv'), sep='\t', index=False)
    calls = []
    for r in Q.itertuples():
        d = r._asdict(); d.pop('Index')
        for tag, db in (('same', 'HG002_' + r.hap), ('other', 'HG002_' + ('PATERNAL' if r.hap == 'MATERNAL' else 'MATERNAL')), ('chm13', 'CHM13')):
            g = LO[(LO['query'] == r.query) & (LO.db == db) & ~LO['self']].sort_values('score', ascending=False)
            if len(g):
                b = g.iloc[0]; d.update({f'{tag}_best_target': b.target, f'{tag}_best_pos': f'{b.target}:{b.tstart}-{b.tend}', f'{tag}_best_class': b.position,
                                         f'{tag}_best_cov': b.qcov, f'{tag}_best_ident': b.ident, f'{tag}_loci_cov50_id98': int(((g.qcov >= 0.5) & (g.ident >= 0.98)).sum())})
            else: d.update({f'{tag}_best_cov': 0.0, f'{tag}_best_ident': 0.0, f'{tag}_loci_cov50_id98': 0})
            if tag == 'other':
                s = LO[(LO['query'] == r.query) & (LO.db == db) & (LO.target == r.contig.replace(r.hap, 'PATERNAL' if r.hap == 'MATERNAL' else 'MATERNAL')) &
                       (LO.position == 'end:' + r.recipient)].sort_values('score', ascending=False)
                d['other_hap_homologous_end_cov'] = s.qcov.iloc[0] if len(s) else 0.0; d['other_hap_homologous_end_ident'] = s.ident.iloc[0] if len(s) else 0.0
        sc = d['same_best_cov'] * d['same_best_ident']; pc = d['chm13_best_cov'] * d['chm13_best_ident']
        hit = d['same_best_cov'] >= 0.5 and d['same_best_ident'] >= 0.98
        d['call'] = 'present_specific' if hit and sc > pc + 0.05 else ('paralogue_level' if hit else 'absent')
        calls.append(d)
    C = pd.DataFrame(calls); C.to_csv(os.path.join(J.RES, 'hg002_lost_sequence_calls.tsv'), sep='\t', index=False)
    print(C[['q', 'form', 'switch_res', 'lost_len', 'same_best_class', 'same_best_cov', 'same_best_ident', 'chm13_best_class', 'chm13_best_cov', 'chm13_best_ident',
             'other_hap_homologous_end_cov', 'other_hap_homologous_end_ident', 'call']].to_string())
    print(C.call.value_counts().to_dict())

if __name__ == '__main__':
    main()
