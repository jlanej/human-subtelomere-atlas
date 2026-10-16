"""Step 2: classify every Tan et al. 2026 window with the catalogue's track-based mosaic method (catalogue revision 2).
  - Alignment: one minimap2 pass per record of the subtelomere window[s0:] against the 48 CHM13 v2.0 500 kb ends
    (pipeline/catalogue_04t_align_tracks.run, subtel_atlas.tracks.MM2_TRACK_ARGS; same reference and parameters as the catalogue).
    The PAF does not depend on paralog families and is cached: work/tan2026/tan_tracks.paf.gz (records aligned first, e.g. the
    same-haplotype partners) plus work/tan2026/tan_tracks_rest.paf.gz (all records not in the first file).
  - Labels: subtel_atlas.tracks.mosaic_end (per-base identity tracks, 100 bp grid anchored at the record tip, labels over the
    published paralog families data/catalogue/paralog_families.tsv, base-level boundaries, S5/S10/S20 = majority grid label).
    Home arm = the Tan arm label mapped to the catalogue arm group. Block columns follow catalogue_06t_mosaic_tracks.label
    (best_end = best single CHM13 end; source = label token, home blocks = best home end).
    Labels are cached together with the SHA-256 of paralog_families.tsv and recomputed when that table changes (no re-alignment).
Outputs: work/tan2026/tan_tracks_blocks.tsv, tan_tracks_ends.tsv (+ .families_sha256)
Usage  : python tan2026_02_mosaic.py [label]     ('label' = re-label from the cached PAFs only; fails if they are missing)
"""
import os, sys, gzip, multiprocessing as mp, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tan2026_common as T
import catalogue_common as C
import catalogue_04t_align_tracks as A4
import catalogue_06t_mosaic_tracks as M6
from subtel_atlas import tracks as TK

FA, WT = f'{T.TW}/tan_tips.fa', f'{T.TW}/tan_windows.tsv'
PAF, REST = f'{T.TW}/tan_tracks.paf.gz', f'{T.TW}/tan_tracks_rest.paf.gz'
FAMS = f'{T.REPO}/data/catalogue/paralog_families.tsv'
OUTB, OUTE = f'{T.TW}/tan_tracks_blocks.tsv', f'{T.TW}/tan_tracks_ends.tsv'
_G = {}

def covered(paf):
    s = set()
    with gzip.open(paf, 'rt') as fh:
        for line in fh: s.add(line.split('\t', 1)[0].rsplit('@', 1)[0])
    return s

def _lab(item):
    q, al = item; e = q.rsplit('@', 1)[0]; w = _G['W'][e]
    return TK.mosaic_end(e, al, w['s0'], w['L'], w['grp'], _G['fam'], M6.NAMES)

def align(label_only=False):
    W = pd.read_csv(WT, sep='\t')
    if not os.path.exists(PAF):
        if label_only: raise SystemExit(f'missing {PAF}')
        A4.run(FA, WT, PAF, T.THREADS)
    if not os.path.exists(REST):
        if label_only: raise SystemExit(f'missing {REST}')
        miss = set(W.end_id) - covered(PAF)
        if miss: A4.run(FA, WT, REST, T.THREADS, only=miss)
    return W

def label(W):
    fsha = C.sha256_file(FAMS)
    if os.path.exists(OUTE) and os.path.exists(OUTE + '.families_sha256') and open(OUTE + '.families_sha256').read().strip() == fsha:
        print('labels cached for families', fsha[:12]); return
    _G['W'] = {r.end_id: dict(s0=int(r.junk_len + r.telo_len), L=int(r.window_len), grp=C.group_of(r.tan_arm)) for r in W.itertuples()}
    _G['fam'] = TK.load_families(FAMS)
    BR, ER, seen = [], [], set()
    def items():
        for paf in (PAF, REST):
            if not os.path.exists(paf): continue
            for q, al in M6.groups(paf):
                e = q.rsplit('@', 1)[0]
                if e in seen: continue
                seen.add(e); yield q, al
        for e in W.end_id:                                   # records without any alignment
            if e not in seen: seen.add(e); yield f'{e}@0', []
    with mp.get_context('fork').Pool(T.THREADS) as pool:
        for bl, rec in pool.imap(_lab, items(), chunksize=16): BR += bl; ER.append(rec)
    B = pd.DataFrame(BR); E = pd.DataFrame(ER)
    B = B.rename(columns={'source': 'best_end'}); B['source'] = B.label
    B.loc[B.is_home, 'source'] = B.loc[B.is_home, 'best_end']
    C.write_tsv_atomic(B.sort_values(['end_id', 'block']), OUTB); C.write_tsv_atomic(E.sort_values('end_id'), OUTE)
    open(OUTE + '.families_sha256', 'w').write(fsha + '\n')
    print('labelled records', E.end_id.nunique())

if __name__ == '__main__':
    lo = len(sys.argv) > 1 and sys.argv[1] == 'label'
    label(align(label_only=lo))
