"""Step 3 - per-read telomere length on the 92 HG002 Q100 v1.1 chromosome ends.

Reads: primary alignments (winnowmap to the diploid v1.1 assembly; T2T HG002 polishing mapping set) overlapping the outer
20 kb of every end, extracted with samtools -F 0x904 (ATLAS_INPUT_DIR/hg002/bams/<platform>.bam; see
hg002_fetch_region.sh in the analysis home and https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/).
Platforms: ont_ul_dorado (R10.4.1 UL, dorado sup) and revio_3cells (PacBio HiFi Revio) [primary]; hifi_dcv1.1 (HiFi
DeepConsensus) reported as a read-length contrast.
Per read (tip frame: chromosome end at position 0):
  boundary  = telomere/subtelomere boundary of the reference end (catalogue_telomere.find_block on the v1.1 end), projected
              into read coordinates through the alignment (<= 500 bp extrapolation if the alignment starts inside the boundary);
  anchor    = read aligned >= 5 kb inward of the boundary (unique subtelomere);
  complete  = read sequence (including soft clips) carries a canonical telomere block starting <= 1,000 bp from the read end
              on the chromosome-end side, >= 50 % canonical (CCCTAA)n coverage between block start and boundary, and >= 300 bp
              of telomere. Reads ending before the boundary or in non-telomeric sequence are not used;
  TL_read   = read bases between the start of the terminal telomere block and the boundary.
A read that is broken inside the telomere is indistinguishable from one that reaches the chromosome end, so per-read lengths
are lower bounds; this truncation is largest for short reads (HiFi) and is quantified by comparing platforms.
Haplotype analyses use MAPQ >= 20 (diploid reference: reads at ends where haplotypes are identical receive low MAPQ).
Inputs : hg002/ref_ends.fa, hg002/hg002v1.1.fasta.gz.fai, hg002/bams/*.bam, results/subtel_hg002_ends.tsv (allele states)
Outputs: work/telomere_cis/hg002_reads.tsv.gz (per read), results/hg002_end_telomere.tsv (per end x platform),
         results/hg002_fixed_window_states.tsv
"""
import os, sys, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, pysam
from tl_common import inp, RES, WORK, rc, find_block, canon_frac, CAN

PLATFORMS = ['ont_ul_dorado', 'revio_3cells', 'hifi_dcv1.1']
ANCHOR = 5000; MAXJUNK = 1000; MINTEL = 300; EXTRAP = 500
fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(inp('hg002', 'hg002v1.1.fasta.gz.fai'))}
ends = {}
name = None
for line in open(inp('hg002', 'ref_ends.fa')):
    if line[0] == '>': name = line[1:].strip(); ends[name] = []
    else: ends[name].append(line.strip())
E = {}
for n, parts in ends.items():
    s = ''.join(parts).upper(); ctg, rg = n.rsplit(':', 1); a, b = map(int, rg.split('-'))
    side = 'L' if a == 1 else 'R'; tf = s if side == 'L' else rc(s)
    blk = find_block(tf); B = blk['junk_len'] + blk['telo_len']
    clen = fai[ctg]; bref = B if side == 'L' else clen - 1 - B
    arm = ctg.split('_')[0].replace('chr', '') + ('p' if side == 'L' else 'q')
    hap = ctg.split('_')[1]
    E[(ctg, side)] = dict(end=f'HG002|{ctg}|{"p" if side == "L" else "q"}', contig=ctg, side=side, arm=arm, hap=hap, clen=clen,
                         ref_junk=blk['junk_len'], ref_telo_len=blk['telo_len'], B=B, bref=bref, tipseq=tf)
rows = []
for plat in PLATFORMS:
    bam = inp('hg002', 'bams', f'{plat}.bam')
    with pysam.AlignmentFile(bam, 'rb', check_sq=False) as bf:
        seen = set()
        for r in bf:
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            ctg = r.reference_name; clen = fai[ctg]
            for side in ('L', 'R'):
                if (ctg, side) not in E: continue
                e = E[(ctg, side)]
                if side == 'L' and r.reference_start >= 20000: continue
                if side == 'R' and r.reference_end <= clen - 20000: continue
                key = (plat, r.query_name, ctg, side)
                if key in seen: continue
                seen.add(key)
                q = r.query_sequence; ql = len(q); bref = e['bref']
                d = dict(platform=plat, read=r.query_name, end=e['end'], mapq=r.mapping_quality, read_len=ql, strand='-' if r.is_reverse else '+')
                # anchor inward of the boundary
                inward = (r.reference_end - bref) if side == 'L' else (bref - r.reference_start)
                d['anchor_bp'] = inward
                if inward < ANCHOR: d['status'] = 'no_anchor'; rows.append(d); continue
                # project boundary to query coordinate
                qpos = None
                if r.reference_start <= bref < r.reference_end:
                    for qp, rp in r.get_aligned_pairs(matches_only=True):
                        if rp >= bref: qpos = qp - (rp - bref) if side == 'L' else qp; break
                else:
                    gap = (r.reference_start - bref) if side == 'L' else (bref - r.reference_end + 1)
                    if 0 < gap <= EXTRAP:
                        qpos = r.query_alignment_start - gap if side == 'L' else r.query_alignment_end - 1 + gap
                if qpos is None: d['status'] = 'boundary_not_covered'; rows.append(d); continue
                qb = qpos if side == 'L' else ql - 1 - qpos     # read bases on the chromosome-end side of the boundary
                d['qb'] = qb
                if qb < MINTEL: d['status'] = 'ends_before_telomere'; rows.append(d); continue
                tf = q if side == 'L' else rc(q)
                blk = find_block(tf[:qb + 200], maxjunk=MAXJUNK)
                if blk is None or blk['junk_len'] > MAXJUNK: d['status'] = 'no_terminal_telomere'; rows.append(d); continue
                j = blk['junk_len']; tl = qb - j; cf = canon_frac(tf[j:qb])
                d.update(read_junk=j, tl=tl, canon_frac=round(cf, 3))
                d['status'] = 'complete' if (tl >= MINTEL and cf >= 0.5) else 'low_canonical'
                rows.append(d)
Rd = pd.DataFrame(rows)
Rd.to_csv(os.path.join(WORK, 'hg002_reads.tsv.gz'), sep='\t', index=False, compression='gzip')
# ---- per end x platform summary
H = pd.read_csv(inp('results', 'subtel_hg002_ends.tsv'), sep='\t')[['q', 't20_state', 't20_best_end']].rename(columns={'q': 'end'})
Ed = pd.DataFrame([{k: v for k, v in e.items() if k != 'tipseq'} for e in E.values()]).merge(H, on='end', how='left')
out = []
for e in Ed.itertuples():
    for plat in PLATFORMS:
        g = Rd[(Rd.end == e.end) & (Rd.platform == plat)]
        c = g[g.status == 'complete']; c20 = c[c.mapq >= 20]
        d = dict(end=e.end, platform=plat, n_reads=len(g), n_complete=len(c), n_complete_mq20=len(c20),
                 status_counts=';'.join(f'{k}:{v}' for k, v in g.status.value_counts().items()))
        for lab, x in (('', c20), ('_allmq', c)):
            if len(x):
                d.update({f'median_tl{lab}': x.tl.median(), f'mean_tl{lab}': x.tl.mean(), f'p75_tl{lab}': x.tl.quantile(.75),
                          f'p90_tl{lab}': x.tl.quantile(.9), f'max_tl{lab}': x.tl.max()})
        if len(c20):
            d['median_read_len'] = c20.read_len.median()
            d['frac_minus_strand'] = (c20.strand == '-').mean()
        out.append(d)
O = Ed.drop(columns=['B', 'bref']).merge(pd.DataFrame(out), on='end')
O.to_csv(os.path.join(RES, 'hg002_end_telomere.tsv'), sep='\t', index=False, float_format='%.5g')
# ---- telomere-length-independent allele state of HG002 ends: S10 = ref tip[s0, s0+10 kb) (same rules as 01b)
#      outer 40 kb of each end are range-read from the public v1.1 FASTA (bgzip + .fai/.gzi) so that S20 = tip[s0, s0+20 kb) exists
S3 = os.environ.get('ATLAS_S3_HTTPS', 'https://human-pangenomics.s3.amazonaws.com')
f40 = os.path.join(WORK, 'hg002_ends40k.fa'); reg = os.path.join(WORK, 'hg002_ends40k.regions')
if not os.path.exists(f40):
    with open(reg, 'w') as fo:
        for (ctg, side), e in sorted(E.items()):
            fo.write(f'{ctg}:1-40000\n' if side == 'L' else f'{ctg}:{e["clen"] - 39999}-{e["clen"]}\n')
    subprocess.run(f'cd {WORK} && samtools faidx {S3}/T2T/HG002/assemblies/hg002v1.1.fasta.gz -r {reg} -o {f40}', shell=True, check=True)
s40 = {}; name = None
for line in open(f40):
    if line[0] == '>': name = line[1:].strip(); s40[name] = []
    else: s40[name].append(line.strip())
for n, parts in s40.items():
    s = ''.join(parts).upper(); ctg, rg = n.rsplit(':', 1); side = 'L' if rg.startswith('1-') else 'R'
    tf = s if side == 'L' else rc(s); assert tf[:20000] == E[(ctg, side)]['tipseq'], n
    E[(ctg, side)]['tip40'] = tf
REF = inp('subtel', 'foreign', 'ref', 'chm13_ends500k.fa'); qfa = os.path.join(WORK, 'hg002_fixed.fa'); paf = os.path.join(WORK, 'hg002_fixed.paf')
with open(qfa, 'w') as fo:
    for e in E.values():
        for tag, L in (('S5', 5000), ('S10', 10000), ('S20', 20000)):
            sub = e['tip40'][e['B']:e['B'] + L]
            if len(sub) >= 0.9 * L: fo.write(f'>{e["end"]}@{tag}\n{sub}\n')
subprocess.run(f'minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 -t {os.environ.get("ATLAS_THREADS", "3")} {REF} {qfa} > {paf} 2> {paf}.log', shell=True, check=True)
best = collections.defaultdict(dict)
for line in open(paf):
    p = line.split('\t')
    if int(p[3]) - int(p[2]) < 0.5 * int(p[1]): continue
    AS = int(next(t[5:] for t in p[12:] if t.startswith('AS:i:')))
    if AS > best[p[0]].get(p[5], -1): best[p[0]][p[5]] = AS
srows = []
for e in E.values():
    rec = dict(end=e['end'], arm=e['arm'], hap=e['hap'], ref_telo_len=e['ref_telo_len'], s0=e['B'])
    for tag in ('S5', 'S10', 'S20'):
        d = best.get(f'{e["end"]}@{tag}')
        if not d: rec[f'{tag}_state'] = 'none'; continue
        be = max(d, key=d.get); b = d[be]; oAS = d.get(f'chr{e["arm"]}', 0)
        rec[f'{tag}_state'] = 'own' if be == f'chr{e["arm"]}' else ('tie' if oAS >= 0.98 * b else 'foreign')
        rec[f'{tag}_best_end'] = be.replace('chr', '')
    srows.append(rec)
S = pd.DataFrame(srows).merge(H, on='end', how='left')
S.to_csv(os.path.join(RES, 'hg002_fixed_window_states.tsv'), sep='\t', index=False)
print(Rd.groupby(['platform', 'status']).size().unstack(fill_value=0).to_string())
print(O.groupby('platform')[['n_complete', 'n_complete_mq20', 'median_tl', 'median_read_len']].median().to_string())
print(pd.crosstab(S.t20_state, S.S20_state)); print(pd.crosstab(S.S10_state, S.S20_state))
