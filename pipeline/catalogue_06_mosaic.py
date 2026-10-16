"""Step 2c: describe every end as an ordered mosaic of CHM13 end segments (telomere -> centromere).
Inputs : work/catalogue/bins.tsv.gz (catalogue_04_align_bins.py), work/catalogue/arms.tsv (home arm per end; catalogue_05_arms.py),
         work/windows/ends.fa + ends_windows.tsv, work/catalogue/chm13_ends500k.fa.
Home ends: the CHM13 end(s) of the end's catalogue arm (catalogue_05_arms.py); for pooled groups (acroP, PAR1, PAR2) every CHM13
   end of the group counts as home.
Bin labels (5 kb tiled bins of the subtelomere = window minus the masked telomere/junk prefix s0):
   home    : a home CHM13 end is the best end or within 2% of the best alignment score (AS) -> shared/own sequence
   <end>   : the best CHM13 end is another end and the home end scores < 0.98 x best (a 'foreign' bin)
   none    : no alignment covering >= 50% of the bin
Smoothing: an isolated bin (single bin) whose two neighbours carry the same label Y is relabelled Y when Y scores >= 0.90 x best in
   that bin, and isolated 'none' bins between identical labels are absorbed. The number of relabelled bins is recorded.
Blocks : maximal runs of one label. Per block: source end, start/end in tip frame and in subtelomere coordinates (tip - s0), mean
   identity of the source end over the block's bins, CHM13 distance-from-telomere span of the source alignments.
Switches: each boundary between two aligned blocks is refined at base level as in the earlier foreign-end calls: the sequence
   within +-10 kb of the bin boundary (inside the two blocks) is aligned (mappy asm20) to both source ends; per-base error tracks
   (minimum over hits) give the changepoint x minimising errors(outer source before x) + errors(inner source after x). Base-level
   if both sides align at >= 95% identity over up to 5 kb and |x - boundary| <= 7.5 kb, else the bin boundary ('window', +-5 kb).
Fixed-window states S5/S10/S20 (boundary-anchored 5/10/20 kb; primary) and t20 (tip[s0:20 kb], as in the earlier calls; kept only
   for comparison - it shrinks as the telomere grows): own / tie (home within 2% of the best foreign end) / foreign / none.
Bins are tiled from the telomere/subtelomere boundary (catalogue_04_align_bins.py), so blocks, switches and structures are
   measured from the boundary and do not depend on telomere length. Note: a 'switch' is a change in which CHM13 end the sequence
   most resembles inside shared (paralogous) subtelomeric sequence, not necessarily a structural breakpoint.
Outputs: work/catalogue/mosaic_blocks.tsv, work/catalogue/mosaic_ends.tsv
Usage  : python catalogue_06_mosaic.py [threads]
Shared code: atlas/subtel_atlas/mosaic.py (labels, blocks, switch refinement, fixed-window states; also used by the query tool).
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas.mosaic import NEAR, SMOOTH, EL, BIN, FIX, RefEnds, mosaic_end   # shared with the query tool

def run(binsf, armsf, ends_fa, wtab, out_prefix, part=0, nparts=1):
    refs = RefEnds(C.readfa(f'{C.CW}/chm13_ends500k.fa'))
    W = pd.read_csv(wtab, sep='\t').set_index('end_id'); A = pd.read_csv(armsf, sep='\t').set_index('end_id')
    Bn = pd.read_csv(binsf, sep='\t', dtype={'bin': str, 'cand': str, 'best_end': str})
    Bn = Bn[Bn.end_id.isin(set(W.index[part::nparts]))]
    fixq = {t: {r.end_id: r for r in Bn[Bn.bin == t].itertuples()} for t in FIX}
    Bw = Bn[~Bn.bin.isin(FIX)].copy(); Bw['k'] = Bw.bin.astype(int); Bw = Bw.sort_values(['end_id', 'k'])
    seqs = {}; want = set(Bw.end_id)
    for n, s in C.readfa(ends_fa):
        if n in want: seqs[n] = s
    BR, ER = [], []
    for e, g in Bw.groupby('end_id', sort=True):
        w = W.loc[e]; s0 = int(w.junk_len + w.telo_len); grp = A.at[e, 'catalogue_arm']
        br, rec = mosaic_end(refs, e, g, {t: fixq[t].get(e) for t in FIX}, s0, grp, seqs[e])
        BR += br; ER.append(rec)
    pd.DataFrame(BR).to_csv(f'{out_prefix}.blocks.tsv', sep='\t', index=False)
    pd.DataFrame(ER).to_csv(f'{out_prefix}.ends.tsv', sep='\t', index=False)

if __name__ == '__main__':
    import subprocess
    if len(sys.argv) > 2 and sys.argv[1] == 'part':
        k, n = int(sys.argv[2]), int(sys.argv[3])
        run(f'{C.CW}/bins.tsv.gz', f'{C.CW}/arms.tsv', f'{C.WIN}/ends.fa', f'{C.WIN}/ends_windows.tsv', f'{C.CW}/mosaic_part{k}', k, n)
    else:
        n = C.THREADS
        # resumable: a part whose two output files exist is not recomputed (delete work/catalogue/mosaic_part* to force)
        bt = os.path.getmtime(f'{C.CW}/bins.tsv.gz')
        done = lambda k: all(os.path.exists(f'{C.CW}/mosaic_part{k}.{x}.tsv') and os.path.getmtime(f'{C.CW}/mosaic_part{k}.{x}.tsv') > bt
                             for x in ('blocks', 'ends'))
        ps = [subprocess.Popen([sys.executable, os.path.abspath(__file__), 'part', str(k), str(n)]) for k in range(n) if not done(k)]
        assert all(p.wait() == 0 for p in ps)
        for kind in ('blocks', 'ends'):
            D = pd.concat([pd.read_csv(f'{C.CW}/mosaic_part{k}.{kind}.tsv', sep='\t') for k in range(n)])
            D = D.sort_values(['end_id', 'block'] if kind == 'blocks' else ['end_id'])
            C.write_tsv_atomic(D, f'{C.CW}/mosaic_{kind}.tsv')
            for k in range(n): os.remove(f'{C.CW}/mosaic_part{k}.{kind}.tsv')
        print('ends', len(D))
