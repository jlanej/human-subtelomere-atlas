"""Compare a (subset) rebuild with the released catalogue, end by end.
Reads the rebuild's intermediates (work/windows/ends_windows.tsv, work/catalogue/arms.tsv, mosaic_ends.tsv, sa_members.tsv,
haplogroup_members.tsv; i.e. steps 03-08b) and the released data/catalogue/ends.tsv given by --release, restricted to the ends of
the rebuild. Every compared column must be identical (numbers compared after rounding to the release precision).
Outputs (--out prefix): <prefix>.tsv (per column: n compared, n identical), <prefix>_mismatches.tsv, <prefix>.json.
Usage: python pipeline/verify_public_rebuild.py --release <path to released ends.tsv> [--out work/public_rebuild_check]"""
import os, sys, json, argparse, numpy as np, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W, CW = os.path.join(REPO, 'work', 'windows'), os.path.join(REPO, 'work', 'catalogue')
COLS = {  # rebuild table -> {rebuild column: release column}
    f'{W}/ends_windows.tsv': dict(window_sha256='window_sha256', window_len='window_len', junk_len='junk_len', telo_len='telo_len',
                                  win_start='win_start', win_end='win_end', window_extended='window_extended', assembly_url='assembly_url'),
    f'{CW}/arms.tsv': dict(home_arm='home_arm', catalogue_arm='catalogue_arm', arm_source='arm_source', arm_confidence='arm_confidence',
                           probe_arm='probe_arm', mosaic_arm='mosaic_arm'),
    f'{CW}/mosaic_ends.tsv': {**dict(profile='profile', structure='structure', switches_sub='switches_sub', switch_res='switch_res'),
                              **{f'{w}_{k}': f'{w}_{k}' for w in ('S5', 'S10', 'S20') for k in ('state', 'family', 'best_end')},
                              'S20_best_ident': 'S20_best_ident'},
    f'{CW}/sa_members.tsv': dict(structure='core_structure', switches='core_switches_sub', allele_id='structural_allele_id',
                                 assignment='sa_assignment', candidates='sa_candidates'),
    f'{CW}/haplogroup_members.tsv': dict(haplogroup_id='haplogroup_id', assignment='hg_assignment', ani_to_medoid='ani_to_medoid')}

def norm(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return ''
    if isinstance(v, (bool, np.bool_)): return str(bool(v))
    if isinstance(v, (int, np.integer)): return str(int(v))
    if isinstance(v, (float, np.floating)): return str(int(v)) if float(v).is_integer() else repr(round(float(v), 6))
    s = str(v)
    return {'True': 'True', 'False': 'False', 'nan': ''}.get(s, s)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--release', required=True); ap.add_argument('--out', default=os.path.join(REPO, 'work', 'public_rebuild_check'))
    a = ap.parse_args()
    R = pd.read_csv(a.release, sep='\t', low_memory=False).set_index('end_id')
    rows, bad = [], []
    for fn, cols in COLS.items():
        T = pd.read_csv(fn, sep='\t', low_memory=False)
        if 'shift' in T: T = T[T['shift'].fillna(0) == 0]
        T = T.drop_duplicates('end_id').set_index('end_id'); ids = [e for e in T.index if e in R.index]
        for c, rc in cols.items():
            n = same = 0
            for e in ids:
                x, y = norm(T.at[e, c]), norm(R.at[e, rc]); n += 1
                if x == y: same += 1
                else: bad.append(dict(table=os.path.basename(fn), column=rc, end_id=e, rebuild=x, release=y))
            rows.append(dict(table=os.path.basename(fn), column=rc, n=n, identical=same))
    S = pd.DataFrame(rows); S.to_csv(a.out + '.tsv', sep='\t', index=False); pd.DataFrame(bad, columns=['table', 'column', 'end_id', 'rebuild', 'release']).to_csv(a.out + '_mismatches.tsv', sep='\t', index=False)
    ends = pd.read_csv(f'{W}/ends_windows.tsv', sep='\t', usecols=['end_id', 'set'])
    summ = dict(n_ends=int(len(ends)), by_set=ends.set.value_counts().to_dict(), columns=int(len(S)), values=int(S.n.sum()),
                identical=int(S.identical.sum()), mismatching=int(len(bad)), all_identical=bool(len(bad) == 0))
    json.dump(summ, open(a.out + '.json', 'w'), indent=1); print(json.dumps(summ))
    sys.exit(0 if not bad else 1)

if __name__ == '__main__':
    main()
