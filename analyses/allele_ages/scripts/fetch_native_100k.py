#!/usr/bin/env python3
"""Region lists for the outer 100 kb of every chromosome-named Release 2 telomeric tip that has no 100 kb window yet
(mostly tips placed on their own arm, i.e. native alleles), for range-reading from the public HPRC S3 assemblies.
Inputs : $ATLAS_INPUT_DIR/tipdup/tips_meta.tsv, $ATLAS_INPUT_DIR/subtel/foreign/r2_100k/*.regions (already fetched),
         $ATLAS_INPUT_DIR/meta/r2_assemblies_index.csv (s3 paths of the Release 2 assemblies)
Outputs: work/allele_ages/r2_100k_native/<asm>.regions and jobs.txt ("asm url")
Then   : bash analyses/allele_ages/scripts/fetch_native_100k.sh
"""
import os, re, glob, pandas as pd
IN = os.environ['ATLAS_INPUT_DIR']; W = 100_000
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
OUT = os.path.join(REPO, 'work', 'allele_ages', 'r2_100k_native'); os.makedirs(OUT, exist_ok=True)
M = pd.read_csv(os.path.join(IN, 'tipdup', 'tips_meta.tsv'), sep='\t'); M = M[M.arm.notna()]
have = set()
for f in glob.glob(os.path.join(IN, 'subtel', 'foreign', 'r2_100k', '*.regions')):
    have |= {l.strip() for l in open(f)}
urls = {}
for line in open(os.path.join(IN, 'meta', 'r2_assemblies_index.csv')):
    u = line.strip().split(',')[-1]; b = os.path.basename(u)
    m = re.match(r'(.+?)_(hap1|hap2|mat|pat)_hprc_r2', b)
    if u.endswith('.fa.gz') and m:
        urls[f'{m.group(1)}.{m.group(2)}'] = u.replace('s3://human-pangenomics/', 'https://human-pangenomics.s3.amazonaws.com/')
jobs, n = [], 0
for asm, g in M.groupby('asm'):
    key = '.'.join(asm.split('.')[:2])
    regs = []
    for r in g.itertuples():
        w = min(W, r.contig_len)
        reg = f'{r.contig}:1-{w}' if r.end == 'L' else f'{r.contig}:{r.contig_len - w + 1}-{r.contig_len}'
        if reg not in have:
            regs.append(reg)
    if not regs:
        continue
    if key not in urls:
        print('no url', asm); continue
    p = os.path.join(OUT, f'{asm}.regions')
    if not os.path.exists(os.path.join(OUT, f'{asm}.ok')):
        open(p, 'w').write('\n'.join(regs) + '\n')
    jobs.append(f'{asm} {urls[key]}'); n += len(regs)
open(os.path.join(OUT, 'jobs.txt'), 'w').write('\n'.join(jobs) + '\n')
print('assemblies', len(jobs), 'regions', n)
