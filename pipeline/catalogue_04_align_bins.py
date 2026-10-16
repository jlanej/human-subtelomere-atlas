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
Shared code: tiling, PAF parsing and summarise() are atlas/subtel_atlas/bins.py (also used by the query tool).
"""
import os, sys, subprocess, gzip, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas.bins import BIN, MINW, EL, MM2_ARGS, COLS, tile, summarise, paf_alignment   # shared with the query tool
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
            for tag, a, b in tile(int(s0[n]), len(s)):
                f.write(f'>{n}@{tag}\n{s[a:b]}\n'); spans[f'{n}@{tag}'] = (a, b)
    return spans

def run(ends_fa, wtab, out, threads):
    make_ref()
    W = pd.read_csv(wtab, sep='\t'); qfa = out + '.queries.fa'; spans = queries(ends_fa, W, qfa)
    cmd = f"minimap2 {' '.join(MM2_ARGS)} -t {threads} {REF} {qfa} 2> {out}.mm2.log"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    cols = COLS
    seen = set()
    with gzip.open(out + '.tmp', 'wt') as fo:
        fo.write('\t'.join(cols) + '\n'); cur, al = None, []
        def flush():
            if cur is None: return
            r = summarise(cur, al, spans[cur]); seen.add(cur); fo.write('\t'.join(str(r[c]) for c in cols) + '\n')
        for line in p.stdout:
            qn, a = paf_alignment(line)
            if qn != cur: flush(); cur, al = qn, []
            al.append(a)
        flush()
        for qn in spans:                       # queries without any alignment
            if qn not in seen: r = summarise(qn, [], spans[qn]); fo.write('\t'.join(str(r[c]) for c in cols) + '\n')
    if p.wait() != 0: raise SystemExit('minimap2 failed')
    os.replace(out + '.tmp', out); os.remove(qfa)

STAMP = f'{C.CW}/bins_stamp.tsv'      # end_id, s0, window_sha256 of the windows the current bins.tsv.gz was computed from

def write_stamp(W):
    W.assign(s0=W.junk_len + W.telo_len)[['end_id', 's0', 'window_sha256']].to_csv(STAMP, sep='\t', index=False)

def update(threads, n_control=300, seed=20261015):
    """Recompute bins only for ends whose boundary (s0) or window changed since bins.tsv.gz was made (bins_stamp.tsv), plus a
    seeded control set of unchanged ends whose recomputed rows must equal the stored rows (minimap2 maps every query
    independently, so a subset run reproduces a full run). Rows are written sorted by end_id (stable within an end)."""
    import numpy as np
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t'); W['s0'] = W.junk_len + W.telo_len
    S = pd.read_csv(STAMP, sep='\t').set_index('end_id')
    m = W.set_index('end_id').join(S, rsuffix='_old')
    changed = sorted(m.index[(m.s0 != m.s0_old) | (m.window_sha256 != m.window_sha256_old)])
    rng = np.random.default_rng(seed); pool = sorted(set(W.end_id) - set(changed))
    control = sorted(rng.choice(pool, min(n_control, len(pool)), replace=False).tolist())
    todo = set(changed) | set(control)
    sub_fa = f'{C.CW}/bins_update.fa'
    with open(sub_fa, 'w') as f:
        for n, sq in C.readfa(f'{C.WIN}/ends.fa'):
            if n in todo: f.write(f'>{n}\n{sq}\n')
    run(sub_fa, f'{C.WIN}/ends_windows.tsv', f'{C.CW}/bins_update.tsv.gz', threads)
    old = pd.read_csv(f'{C.CW}/bins.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    new = pd.read_csv(f'{C.CW}/bins_update.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    key = lambda d: d.assign(_o=np.arange(len(d))).sort_values(['end_id', '_o']).drop(columns='_o').reset_index(drop=True)
    oc, nc = key(old[old.end_id.isin(control)]), key(new[new.end_id.isin(control)])
    same = oc.equals(nc)
    print(f'bins update: changed ends {len(changed)}, control ends {len(control)}, control rows identical: {same} ({len(nc)} rows)')
    if not same: raise SystemExit('control ends differ from stored bins; run a full step 04')
    out = pd.concat([old[~old.end_id.isin(changed)], new[new.end_id.isin(changed)]])
    out = key(out)
    with gzip.open(f'{C.CW}/bins.tsv.gz.tmp', 'wt') as fo: out.to_csv(fo, sep='\t', index=False)
    os.replace(f'{C.CW}/bins.tsv.gz.tmp', f'{C.CW}/bins.tsv.gz')
    write_stamp(W.drop(columns='s0'))
    pd.DataFrame(dict(end_id=changed)).to_csv(f'{C.CW}/bins_update_changed.tsv', sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'update':
        update(C.THREADS)
    else:
        run(f'{C.WIN}/ends.fa', f'{C.WIN}/ends_windows.tsv', f'{C.CW}/bins.tsv.gz', C.THREADS)
        write_stamp(pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t'))
