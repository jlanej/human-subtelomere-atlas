"""Step 1: stream every record of the Tan et al. 2026 archive and write a telomere-anchored window in the catalogue's tip frame
(telomere at position 0, reading inward; a correctly oriented telomere reads (CCCTAA)n).
Orientation: the telomere block (pipeline/catalogue_telomere.py, same caller and parameters as the catalogue) is called at the
  start of the record as delivered; if absent, at the start of the reverse complement; if absent in both, the record is kept as
  delivered with junk_len = telo_len = 0 and has_telomere = False (excluded from primary analyses).
Window: tip[0 : s0 + 100,000] with s0 = junk_len + telo_len, i.e. the first 100 kb of subtelomere from the telomere boundary
  (shorter if the record is shorter).
Outputs: work/tan2026/tips/<arm>.fa and <arm>.tsv (per archive member), then work/tan2026/tan_tips.fa and tan_windows.tsv.
Usage  : python tan2026_01_extract.py
"""
import os, sys, zipfile, gzip, hashlib, multiprocessing as mp, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
from catalogue_common import rc
from catalogue_telomere import find_block

def records(member):
    with zipfile.ZipFile(T.ZIP) as Z, Z.open(member) as f, gzip.open(f, 'rt') as g:
        name, buf = None, []
        for line in g:
            if line[0] == '>':
                if name is not None: yield name, ''.join(buf)
                name, buf = line[1:].strip().split()[0], []
            else: buf.append(line.strip())
        if name is not None: yield name, ''.join(buf)

def one(member):
    arm = os.path.basename(member).split('.')[0]; out = f'{T.TIPS}/{arm}'
    if os.path.exists(out + '.tsv'): return out
    rows = []
    with open(out + '.fa.tmp', 'w') as fo:
        for n, s in records(member):
            s = s.upper(); L = len(s); b = find_block(s[:70000]); ori = 'as_given'
            if b is None:
                b2 = find_block(rc(s[-70000:]))
                if b2 is not None: s, b, ori = rc(s), b2, 'reverse_complemented'
            has = b is not None; jl, tl = (b['junk_len'], b['telo_len']) if has else (0, 0)
            w = s[:jl + tl + T.SUB]
            fo.write(f'>{n}\n{w}\n')
            d = T.parse_name(n)
            rows.append(dict(end_id=n, sample=d['sample'], hap=d['hap'], cohort=d['cohort'], assembler=d['asm'], tan_arm=d['arm'],
                             contig=d['contig'], rec_start=int(d['start']), rec_end=int(d['end']), rec_ori=d['ori'], rec_len=L,
                             orient=ori, has_telomere=has, junk_len=jl, telo_len=tl, window_len=len(w),
                             window_sha256=hashlib.sha256(w.encode()).hexdigest()))
    os.replace(out + '.fa.tmp', out + '.fa'); pd.DataFrame(rows).to_csv(out + '.tsv.tmp', sep='\t', index=False)
    os.replace(out + '.tsv.tmp', out + '.tsv'); return out

def main():
    with zipfile.ZipFile(T.ZIP) as Z:
        mem = sorted((i for i in Z.infolist() if i.filename.endswith('.fasta.gz')), key=lambda i: -i.file_size)
    with mp.Pool(T.THREADS) as p:
        outs = p.map(one, [m.filename for m in mem], chunksize=1)
    outs = sorted(outs)
    D = pd.concat([pd.read_csv(o + '.tsv', sep='\t') for o in outs]).sort_values('end_id')
    with open(f'{T.TW}/tan_tips.fa.tmp', 'w') as fo:
        for o in outs:
            with open(o + '.fa') as fi:
                for line in fi: fo.write(line)
    os.replace(f'{T.TW}/tan_tips.fa.tmp', f'{T.TW}/tan_tips.fa')
    D.to_csv(f'{T.TW}/tan_windows.tsv', sep='\t', index=False)
    print('records', len(D), 'with telomere', int(D.has_telomere.sum()), D.orient.value_counts().to_dict())

if __name__ == '__main__':
    main()
