"""Step 2 prerequisite - telomere-length-independent subtelomere allele states.

Problem: the published t20 state classifies tip[s0:20000] (s0 = terminal junk + telomere), so the amount of subtelomere that is
classified shrinks as the assembled telomere grows (20 kb - s0). Allele state and telomere length are therefore not
independent by construction. Here every named tip is re-classified on FIXED subtelomere-anchored windows that start at the
telomere boundary: S5 = tip[s0, s0+5 kb), S10 = tip[s0, s0+10 kb), S20 = tip[s0, s0+20 kb). A replica of the published query
(T20R = tip[s0_published, 20 kb)) is run in the same job to check the re-implementation against results/subtel_foreign_tips.tsv.
Classification rules are those of the published placement (subtel/foreign scripts): minimap2 -x asm20 -c --secondary=yes
-N 50 -p 0.2 against the outer 500 kb of the 48 CHM13 v2.0 ends; alignments covering >= 50 % of the query; best AS per end;
own = best end is own; tie = best end foreign but own AS >= 0.98 x best; foreign = otherwise; none = no alignment.
Inputs : 100 kb tip windows from the catalogue pipeline (work/windows/fetch/<asm>.fa, pipeline/catalogue_02_fetch_windows.sh)
         and ATLAS_INPUT_DIR/subtel/foreign/r2_100k/<asm>.fa (genome orientation; R ends reverse-complemented here),
         work/windows/r2_plan.tsv (tip -> region), results/tip_telomere_lengths.tsv (step 1), CHM13 ends
         ATLAS_INPUT_DIR/subtel/foreign/ref/chm13_ends500k.fa.
Outputs: work/telomere_cis/fixed_windows.{fa,paf} (intermediates), results/fixed_window_states.tsv
"""
import os, sys, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from tl_common import inp, RES, WORK, REPO, rc

THREADS = int(os.environ.get('ATLAS_THREADS', '3')); NEAR = 0.98
REF = inp('subtel', 'foreign', 'ref', 'chm13_ends500k.fa')
T = pd.read_csv(os.path.join(RES, 'tip_telomere_lengths.tsv'), sep='\t', low_memory=False)
T['s0'] = T.junk_len + T.telo_len; T['s0_pub'] = T.junk_len + T.telo_len_raw20
plan = pd.read_csv(os.path.join(REPO, 'work', 'windows', 'r2_plan.tsv'), sep='\t').set_index('q')
WIN = {'S5': 5000, 'S10': 10000, 'S20': 20000}
qfa = os.path.join(WORK, 'fixed_windows.fa'); paf = os.path.join(WORK, 'fixed_windows.paf')
def readfa(fn):
    name = None; buf = []
    for line in open(fn):
        if line[0] == '>':
            if name: yield name, ''.join(buf)
            name = line[1:].strip().split()[0]; buf = []
        else: buf.append(line.strip())
    if name: yield name, ''.join(buf)
qlen = {}; nmiss = 0; nshort = collections.Counter()
if not (os.path.exists(paf) and os.path.getsize(paf) > 0 and os.environ.get('REUSE_PAF')):
    with open(qfa, 'w') as fo:
        for asm, g in T.groupby('asm'):
            src = {}
            for fn in (os.path.join(REPO, 'work', 'windows', 'fetch', f'{asm}.fa'), inp('subtel', 'foreign', 'r2_100k', f'{asm}.fa')):
                if os.path.exists(fn):
                    for n, s in readfa(fn): src.setdefault(n, s)
            for r in g.itertuples():
                if r.q not in plan.index or plan.at[r.q, 'region'] not in src: nmiss += 1; continue
                s = src[plan.at[r.q, 'region']].upper(); s = s if r.end == 'L' else rc(s)
                s0 = int(r.s0)
                for tag, L in WIN.items():
                    sub = s[s0:s0 + L]
                    if len(sub) < L: nshort[tag] += 1
                    if len(sub) >= 0.9 * L: fo.write(f'>{r.q}@{tag}\n{sub}\n')
                sub = s[int(r.s0_pub):20000]
                if len(sub) >= 1000: fo.write(f'>{r.q}@T20R\n{sub}\n')
    print('tips without a 100 kb window:', nmiss, '| windows shorter than nominal:', dict(nshort))
    subprocess.run(f'minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.2 -t {THREADS} {REF} {qfa} > {paf} 2> {paf}.log', shell=True, check=True)
best = collections.defaultdict(dict)   # query -> end -> best AS
for line in open(paf):
    p = line.split('\t'); ql, qs, qe = int(p[1]), int(p[2]), int(p[3])
    if qe - qs < 0.5 * ql: continue
    AS = int(next(t[5:] for t in p[12:] if t.startswith('AS:i:'))); d = best[p[0]]
    if AS > d.get(p[5], -1): d[p[5]] = AS
rows = []
own = dict(zip(T.q, T.arm))
for q in T.q:
    rec = dict(q=q)
    for tag in list(WIN) + ['T20R']:
        d = best.get(f'{q}@{tag}')
        if not d: rec[f'{tag}_state'] = 'none'; continue
        be = max(d, key=d.get); b = d[be]; oe = f'chr{own[q]}'; oAS = d.get(oe, 0)
        st = 'own' if be == oe else ('tie' if oAS >= NEAR * b else 'foreign')
        rec.update({f'{tag}_state': st, f'{tag}_best_end': be.replace('chr', ''), f'{tag}_n_near': sum(v >= NEAR * b for v in d.values()) - 1,
                    f'{tag}_own_frac': round(oAS / b, 4)})
    rows.append(rec)
S = pd.DataFrame(rows).merge(T[['q', 't20_state', 'foreign_arm']], on='q')
S.to_csv(os.path.join(RES, 'fixed_window_states.tsv'), sep='\t', index=False)
agree = (S.T20R_state == S.t20_state).mean()
print('T20R replica agreement with published t20_state:', round(agree, 4), '| n', len(S))
print(pd.crosstab(S.t20_state, S.S10_state))
for tag in WIN: print(tag, S[f'{tag}_state'].value_counts().to_dict())
