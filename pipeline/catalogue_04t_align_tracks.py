"""Step 2a': one alignment pass per end for the track-based mosaic (atlas/subtel_atlas/tracks.py).
Query = the subtelomere window[s0:] (s0 = junk_len + telo_len, optionally shifted by d for the boundary-shift test: window[s0+d:],
named <end_id>@<d>). Aligned to the 48 CHM13 v2.0 500 kb ends with minimap2 tracks.MM2_TRACK_ARGS; PAF (with cs) gzip-compressed.
Outputs: work/catalogue/tracks.paf.gz (catalogue); with 'shift': work/catalogue/tracks_shift.paf.gz (the ends of
         work/catalogue/stability_subset.tsv at the shifts of catalogue_12_stability_eval.SHIFTS).
Usage  : python catalogue_04t_align_tracks.py [shift]
"""
import os, sys, gzip, subprocess, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas.tracks import MM2_TRACK_ARGS
REF = f'{C.CW}/chm13_ends500k.fa'

def run(ends_fa, wtab, out, threads, shifts=(0,), only=None):
    W = pd.read_csv(wtab, sep='\t').set_index('end_id'); s0 = (W.junk_len + W.telo_len).astype(int)
    qfa = out + '.queries.fa'
    with open(qfa, 'w') as f:
        for n, s in C.readfa(ends_fa):
            if only is not None and n not in only: continue
            for d in shifts:
                a = int(s0[n]) + d
                if 0 <= a < len(s) - 1000: f.write(f'>{n}@{d}\n{s[a:]}\n')
    cmd = ['minimap2', *MM2_TRACK_ARGS, '-t', str(threads), REF, qfa]
    with gzip.open(out + '.tmp', 'wt', compresslevel=3) as fo:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open(out + '.mm2.log', 'w'), text=True, bufsize=1 << 20)
        for line in p.stdout: fo.write(line)
        if p.wait() != 0: raise SystemExit('minimap2 failed')
    os.replace(out + '.tmp', out)
    os.replace(qfa, out + '.queries.fa.done')        # kept (no deletion); overwritten by the next run

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'shift':
        import catalogue_12_stability_eval as SE
        S = pd.read_csv(f'{C.CW}/stability_subset.tsv', sep='\t')
        run(f'{C.WIN}/ends.fa', f'{C.WIN}/ends_windows.tsv', f'{C.CW}/tracks_shift.paf.gz', C.THREADS, shifts=SE.SHIFTS, only=set(S.end_id))
    else:
        run(f'{C.WIN}/ends.fa', f'{C.WIN}/ends_windows.tsv', f'{C.CW}/tracks.paf.gz', C.THREADS)
