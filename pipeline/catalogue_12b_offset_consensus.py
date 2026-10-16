"""Comparison for the redesign: the previous 5 kb-bin method made phase-robust by consensus over OFFSETS tiling offsets (the
alternative design), versus the track-based method, on a seeded subsample of the stability subset at boundary shifts -500/0/+500.
Bins (old method: subtel_atlas.bins tile/summarise, minimap2 asm20) are tiled from s0 + d + offset; bin labels (old rules,
subtel_atlas.mosaic.label_bins) are mapped to paralog families (the published families); for every tip-frame 1 kb position the
consensus label is the majority over offsets whose bins cover it; runs < 2 kb are merged into the longer neighbour; core structure
via alleles.core. Output: work/catalogue/offset_consensus_compare.tsv"""
import os, sys, types, collections, subprocess, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C
from subtel_atlas.bins import MM2_ARGS, summarise, paf_alignment, BIN, MINW
from subtel_atlas import mosaic as MO, alleles as AL, tracks as TK
N_SUB, SEED, OFFSETS, SHIFTS = 300, 20261018, (0, 1000, 2000, 3000, 4000), (-500, 0, 500)
REF = f'{C.CW}/chm13_ends500k.fa'

def main():
    S = pd.read_csv(f'{C.CW}/stability_subset.tsv', sep='\t')
    sub = sorted(np.random.default_rng(SEED).choice(sorted(S.end_id), N_SUB, replace=False).tolist())
    W = pd.read_csv(f'{C.WIN}/ends_windows.tsv', sep='\t').set_index('end_id'); A = pd.read_csv(f'{C.CW}/arms.tsv', sep='\t').set_index('end_id')
    fam = TK.load_families(f'{C.CAT}/paralog_families.tsv')
    seqs = {n: s for n, s in C.readfa(f'{C.WIN}/ends.fa') if n in set(sub)}
    q, spans = [], {}
    for e in sub:
        s0 = int(W.at[e, 'junk_len'] + W.at[e, 'telo_len']); L = len(seqs[e])
        for d in SHIFTS:
            for o in OFFSETS:
                a0 = s0 + d + o; k = 0
                while a0 + k * BIN < L:
                    a, b = a0 + k * BIN, min(a0 + (k + 1) * BIN, L)
                    if b - a >= MINW: qn = f'{e}|{d}|{o}|{k}'; q.append((qn, seqs[e][a:b])); spans[qn] = (a, b)
                    k += 1
    qfa = f'{C.CW}/offset_consensus.queries.fa'
    with open(qfa, 'w') as f:
        for qn, s in q: f.write(f'>{qn}\n{s}\n')
    paf = f'{C.CW}/offset_consensus.paf'
    if not os.path.exists(paf):
        p = subprocess.run(['minimap2', *MM2_ARGS, '-t', str(C.THREADS), REF, qfa], capture_output=True, text=True, check=True)
        open(paf + '.tmp', 'w').write(p.stdout); os.replace(paf + '.tmp', paf)
    al = collections.defaultdict(list)
    for line in open(paf):
        qn, a = paf_alignment(line); al[qn].append(a)
    byk = collections.defaultdict(list)
    for qn in spans:
        e_, d_, o_, k_ = qn.split('|'); byk[(e_, int(d_), int(o_))].append(qn)
    rows = []
    for e in sub:
        s0 = int(W.at[e, 'junk_len'] + W.at[e, 'telo_len']); L = len(seqs[e]); homes = C.group_ends(A.at[e, 'catalogue_arm'])
        for d in SHIFTS:
            votes = collections.defaultdict(collections.Counter)
            for o in OFFSETS:
                bl = []
                for qn in byk[(e, d, o)]:
                    r = summarise(f"{e}@{qn.rsplit('|', 1)[1]}", al.get(qn, []), spans[qn]); r['qstart'], r['qend'] = spans[qn]; bl.append(r)
                g = pd.DataFrame(bl).sort_values('qstart'); g['best_end'] = g.best_end.replace('', np.nan)
                labs, _ = MO.label_bins(g, homes)
                for lab, cd, r, he in labs:
                    tok = lab if lab in ('home', 'none') else TK.family_of(fam, lab)
                    for x in range(-(-int(r.qstart) // 1000) * 1000, int(r.qend), 1000): votes[x][tok] += 1
            xs = sorted(x for x in votes if x >= s0 + d)
            runs = []
            for x in xs:
                tok = sorted(votes[x].items(), key=lambda t: (-t[1], t[0]))[0][0]
                if runs and runs[-1][0] == tok: runs[-1][2] = x + 1000
                else: runs.append([tok, x, x + 1000])
            while len(runs) > 1:
                ln = [r[2] - r[1] for r in runs]; i = int(np.argmin(ln))
                if ln[i] >= 2000: break
                k = i - 1 if (i == len(runs) - 1 or (i > 0 and ln[i - 1] >= ln[i + 1])) else i + 1
                if k < i: runs[k][2] = runs[i][2]
                else: runs[k][1] = runs[i][1]
                del runs[i]
                j = 1
                while j < len(runs):
                    if runs[j][0] == runs[j - 1][0]: runs[j - 1][2] = runs[j][2]; del runs[j]
                    else: j += 1
            sb = s0 + d
            blk = [types.SimpleNamespace(is_home=r[0] == 'home', source=r[0], sub_start=max(r[1], sb) - sb, sub_end=r[2] - sb) for r in runs]
            _, st, sw = AL.core(blk) if blk else (None, '', [])
            rows.append(dict(end_id=e, shift=d, structure_consensus=st))
    R = pd.DataFrame(rows)
    T = pd.read_csv(f'{C.CW}/stability_results.tsv', sep='\t')[['end_id', 'shift', 'structure']].rename(columns={'structure': 'structure_tracks'})
    R = R.merge(T, on=['end_id', 'shift'], how='left'); R.to_csv(f'{C.CW}/offset_consensus_compare.tsv', sep='\t', index=False)
    z = R[R['shift'] == 0].set_index('end_id'); nat = S.set_index('end_id').native
    for d in (-500, 500):
        x = R[R['shift'] == d].set_index('end_id').reindex(z.index)
        for m in ('structure_consensus', 'structure_tracks'):
            same = x[m] == z[m]; nn = ~nat.reindex(z.index).astype(bool)
            print(d, m, 'all', round(same.mean(), 4), 'non-native', round(same[nn].mean(), 4), int(nn.sum()), 'native', round(same[~nn].mean(), 4))

if __name__ == '__main__':
    main()
