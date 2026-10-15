"""Step 1a: plan telomere-anchored 100 kb windows for every Release 2 telomeric tip (named and unnamed contigs).
Inputs : $ATLAS_INPUT_DIR/tipdup/tips_meta.tsv (20,144 tips), meta/r2_assemblies_index.csv (S3 paths),
         subtel/foreign/r2_100k/<asm>.fa (outer 100 kb windows fetched earlier; reused, not re-downloaded).
Window : contig bases [1, min(100000, contig_len)] for L ends, [contig_len - w + 1, contig_len] for R ends (1-based, inclusive).
Outputs: work/windows/r2_plan.tsv (one row per tip: end_id, region, url, reuse source), work/windows/fetch/<asm>.regions +
         work/windows/fetch/jobs.txt for the regions not already present locally.
"""
import os, sys, re, glob, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import catalogue_common as C
os.makedirs(f'{C.WIN}/fetch', exist_ok=True)
M = pd.read_csv(C.inp('tipdup/tips_meta.tsv'), sep='\t')
urls = {}
for line in open(C.inp('meta/r2_assemblies_index.csv')):
    u = line.strip().split(',')[-1]; b = os.path.basename(u)
    m = re.match(r'(.+?)_(hap1|hap2|mat|pat)_hprc_r2', b)
    if u.endswith('.fa.gz') and m: urls[f'{m.group(1)}.{m.group(2)}'] = u
have = {}
for fn in sorted(glob.glob(C.inp('subtel/foreign/r2_100k/*.fa'))):
    if not os.path.exists(fn[:-3] + '.ok'): continue
    asm = os.path.basename(fn)[:-3]
    for line in open(fn):
        if line[0] == '>': have[(asm, line[1:].strip())] = 'r2_100k'
def region(r):
    w = min(C.WINDOW, int(r.contig_len))
    return (1, w) if r.end == 'L' else (int(r.contig_len) - w + 1, int(r.contig_len))
rows = []
for r in M.itertuples():
    a, b = region(r); reg = f'{r.contig}:{a}-{b}'; key = '.'.join(r.asm.split('.')[:2])
    rows.append(dict(end_id=f'{r.contig}#{r.end}', q=r.q, asm=r.asm, contig=r.contig, side=r.end, contig_len=int(r.contig_len),
                     win_start=a, win_end=b, region=reg, s3_path=urls[key],
                     source='r2_100k' if (r.asm, reg) in have else 'fetch'))
P = pd.DataFrame(rows)
assert P.end_id.is_unique
C.write_tsv_atomic(P, f'{C.WIN}/r2_plan.tsv')
jobs = []
for asm, g in P[P.source == 'fetch'].groupby('asm'):
    regs = sorted(set(g.region))
    fn = f'{C.WIN}/fetch/{asm}.regions'
    if not os.path.exists(fn) or open(fn).read() != '\n'.join(regs) + '\n':
        open(fn, 'w').write('\n'.join(regs) + '\n')
    jobs.append(f"{asm} {g.s3_path.iloc[0].replace('s3://human-pangenomics', C.S3)}")
open(f'{C.WIN}/fetch/jobs.txt', 'w').write('\n'.join(jobs) + '\n')
print('tips', len(P), 'reuse', (P.source == 'r2_100k').sum(), 'to fetch', (P.source == 'fetch').sum(),
      'unique regions to fetch', P[P.source == 'fetch'].region.nunique(), 'assemblies', len(jobs))

# ---- arm-inference probes: 10 kb at 250 kb and 1 Mb from the tip (tip frame) for every tip on an unnamed (JA/JB) contig and,
# as a calibration set, for 1,000 randomly chosen named tips (seed 20260715). A probe is planned only if it ends inside the
# contig's tip-side half (so that it cannot reach the other end).
PROBES = [(250_000, 10_000), (1_000_000, 10_000)]
os.makedirs(f'{C.WIN}/probes', exist_ok=True)
named = M[M.arm.notna()].sample(1000, random_state=20260715).q
sel = M[M.arm.isna() | M.q.isin(named)]
prow = []
for r in sel.itertuples():
    L = int(r.contig_len)
    for off, ln in PROBES:
        if off + ln > L // 2: continue
        a, b = (off + 1, off + ln) if r.end == 'L' else (L - off - ln + 1, L - off)
        prow.append(dict(end_id=f'{r.contig}#{r.end}', asm=r.asm, probe_offset=off, region=f'{r.contig}:{a}-{b}',
                         calibration=isinstance(r.arm, str)))
PR = pd.DataFrame(prow); C.write_tsv_atomic(PR, f'{C.WIN}/probe_plan.tsv')
pj = []
for asm, g in PR.groupby('asm'):
    open(f'{C.WIN}/probes/{asm}.regions', 'w').write('\n'.join(sorted(set(g.region))) + '\n')
    pj.append(f"{asm} {P[P.asm == asm].s3_path.iloc[0].replace('s3://human-pangenomics', C.S3)}")
open(f'{C.WIN}/probes/jobs.txt', 'w').write('\n'.join(pj) + '\n')
print('probe regions', len(PR), 'tips with >=1 probe', PR.end_id.nunique(), 'of', len(sel), 'assemblies', len(pj))
