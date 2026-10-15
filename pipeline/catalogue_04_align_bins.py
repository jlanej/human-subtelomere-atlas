"""Step 2a: align tiled 5 kb bins of every end window (and the outer-20 kb subtelomere 't20') to the 48 CHM13 v2.0 500 kb ends.
Reference: work/catalogue/chm13_ends500k.fa = outer 500 kb of each CHM13 v2.0 chromosome end (chr1p..chrYq) in genome orientation
           (p ends: telomere at 0; q ends: telomere at 500,000), extracted here from $ATLAS_INPUT_DIR/ref/chm13v2.0.fa.gz.
Queries  : per end, s0 = junk_len + telo_len (telomere + outer junk, excluded). All queries are anchored at the telomere/subtelomere
           boundary s0, so no call depends on telomere length:
             bins k = sub[5k*k, 5k*(k+1)) with sub = tip[s0:] (the last bin kept if >= 2 kb);
             S5 / S10 / S20 = sub[0:5000] / sub[0:10000] / sub[0:20000] (boundary-anchored fixed windows);
           and, for comparison with the earlier calls only, 't20' = tip[s0:20000] (shrinks as the telomere grows).
Aligner  : minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 (as in the earlier calls). Per query and CHM13 end: best alignment score
           (AS) among alignments covering >= 50% of the query. The PAF is streamed and summarised (not stored).
Output   : work/catalogue/bins.tsv.gz: end_id, bin (k, S5, S10, S20 or t20), qstart, qend (tip frame; subtract s0 for boundary coordinates), best_end, best_AS, best_ident, best_dist0,
           best_dist1 (distance of the best alignment from that CHM13 end's telomere), best_strand, cand (all ends with AS >= 0.5 x best:
           'end:AS:ident' joined by ';'), n_ends
Usage    : python catalogue_04_align_bins.py
"""
import os, sys, subprocess, gzip, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
BIN, MINW, EL = 5000, 2000, 500_000
os.makedirs(C.CW, exist_ok=True)
REF = f'{C.CW}/chm13_ends500k.fa'

def make_ref():
    if os.path.exists(REF + '.fai'): return
    regs, names = [], []
    for line in open(C.inp('ref/chm13v2.0.fa.gz.fai')):
        c, L = line.split('\t')[:2]; L = int(L)
        if c == 'chrM': continue
        regs += [f'{c}:1-{EL}', f'{c}:{L - EL + 1}-{L}']; names += [f'{c}p', f'{c}q']
    open(f'{C.CW}/chm13_ends500k.regions', 'w').write('\n'.join(regs) + '\n')
    fa = subprocess.run(f"samtools faidx -r {C.CW}/chm13_ends500k.regions {C.inp('ref/chm13v2.0.fa.gz')}", shell=True, check=True,
                        capture_output=True, text=True).stdout
    it = iter(names); lines = [('>' + next(it)) if l.startswith('>') else l.upper() for l in fa.splitlines()]
    open(REF + '.tmp', 'w').write('\n'.join(lines) + '\n'); os.replace(REF + '.tmp', REF)
    subprocess.run(f'samtools faidx {REF}', shell=True, check=True)

def queries(ends_fa, W, out):
    s0 = dict(zip(W.end_id, W.junk_len + W.telo_len)); spans = {}
    with open(out, 'w') as f:
        for n, s in C.readfa(ends_fa):
            a0 = int(s0[n]); L = len(s)
            for k in range((L - a0 + BIN - 1) // BIN):
                a, b = a0 + k * BIN, min(a0 + (k + 1) * BIN, L)
                if b - a >= MINW: f.write(f'>{n}@{k}\n{s[a:b]}\n'); spans[f'{n}@{k}'] = (a, b)
            for tag, w in (('S5', 5000), ('S10', 10000), ('S20', 20000)):
                b = min(a0 + w, L)
                if b - a0 >= 1000: f.write(f'>{n}@{tag}\n{s[a0:b]}\n'); spans[f'{n}@{tag}'] = (a0, b)
            b = min(20000, L)
            if b - a0 >= 1000: f.write(f'>{n}@t20\n{s[a0:b]}\n'); spans[f'{n}@t20'] = (a0, b)
    return spans

def summarise(qn, al, span):
    best = {}
    for a in al:
        if a[2] - a[1] < 0.5 * a[0]: continue
        if a[4] not in best or a[7] > best[a[4]][7]: best[a[4]] = a
    end_id, k = qn.rsplit('@', 1)
    rec = dict(end_id=end_id, bin=k, qstart=span[0], qend=span[1])
    if not best:
        rec.update(best_end='', best_AS=0, best_ident='', best_dist0='', best_dist1='', best_strand='', cand='', n_ends=0); return rec
    b = max(best.values(), key=lambda a: (a[7], a[4]))
    d0, d1 = (b[5], b[6]) if b[4].endswith('p') else (EL - b[6], EL - b[5])
    cand = sorted((a for a in best.values() if a[7] >= 0.5 * b[7]), key=lambda a: (-a[7], a[4]))
    rec.update(best_end=b[4], best_AS=b[7], best_ident=round(b[8], 4), best_dist0=d0, best_dist1=d1, best_strand=b[3],
               cand=';'.join(f'{a[4]}:{a[7]}:{a[8]:.4f}' for a in cand), n_ends=len(best))
    return rec

def run(ends_fa, wtab, out, threads):
    make_ref()
    W = pd.read_csv(wtab, sep='\t'); qfa = out + '.queries.fa'; spans = queries(ends_fa, W, qfa)
    cmd = f'minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 -t {threads} {REF} {qfa} 2> {out}.mm2.log'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    cols = ['end_id', 'bin', 'qstart', 'qend', 'best_end', 'best_AS', 'best_ident', 'best_dist0', 'best_dist1', 'best_strand', 'cand', 'n_ends']
    seen = set()
    with gzip.open(out + '.tmp', 'wt') as fo:
        fo.write('\t'.join(cols) + '\n'); cur, al = None, []
        def flush():
            if cur is None: return
            r = summarise(cur, al, spans[cur]); seen.add(cur); fo.write('\t'.join(str(r[c]) for c in cols) + '\n')
        for line in p.stdout:
            f = line.split('\t')
            if f[0] != cur: flush(); cur, al = f[0], []
            AS = next(int(t[5:]) for t in f[12:] if t.startswith('AS:i:'))
            de = next(float(t[5:]) for t in f[12:] if t.startswith('de:f:'))
            al.append((int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), AS, 1 - de))
        flush()
        for qn in spans:                       # queries without any alignment
            if qn not in seen: r = summarise(qn, [], spans[qn]); fo.write('\t'.join(str(r[c]) for c in cols) + '\n')
    if p.wait() != 0: raise SystemExit('minimap2 failed')
    os.replace(out + '.tmp', out); os.remove(qfa)

if __name__ == '__main__':
    run(f'{C.WIN}/ends.fa', f'{C.WIN}/ends_windows.tsv', f'{C.CW}/bins.tsv.gz', C.THREADS)
