#!/usr/bin/env python3
"""Root-cause diagnostics for the 1000 Genomes ONT population evidence (evidence/population/DIAGNOSTICS.md).

Diagnosis only: the frozen criteria v1 and subtel_atlas.readevidence are used unchanged; no verdict in claims.tsv is modified.
Needs the outputs of run_population.py (work/evidence/population: plan.json, reads_full.tsv, units.json.gz,
assembly_projection.json, slices/) and the end windows work/windows/ends.fa.

Analyses (tables in evidence/population/diagnostics/):
  conclusive_read_placement.tsv     every conclusive read aligned (minimap2 map-ont via mappy) to all of the individual's Release 2
                                    end windows; best end = most aligned matching bases
  confirmed_claims_annotation.tsv   confirmed claims: reads with a unique anchor, reads placed on the claimed arm, annotation
  contradicting_read_placement.tsv  contradicting reads of 'contradicting_only' claims placed the same way
  sa3q01_discordant_reads.tsv       SA-3q-01 conclusive reads: placement, and the claim model applied to the best-placing window
  sa3q01_window_projections.tsv     claim model applied to both 3q haplotype windows of each individual
  sa3q01_alt_donor_15q.tsv          the same reads and windows with donor CHM13 15q instead of 11p
  donor_window_edge_screen.tsv      claims whose assembly switch coincides with the end of the D window while CHM13 D and P stay
                                    homologous beyond it (P 'private' bins there are a window artefact)
  sa20p03_wide60kb_reads.tsv        chr20p outer 60 kb (8 carriers) evaluated with the SA-20p-03 models
  sa20p03_donor_slice_reads_placement.tsv  chr7p donor-slice reads of three carriers, placement on the person's windows
  claim_annotations.tsv             per-claim diagnostic notes and a recommendation (reporting only; verdicts unchanged)
Figures: fig_sa3q01_tracks.(png|pdf), fig_sa20p03_tracks.(png|pdf).
"""
import os, sys, json, gzip, time, collections, argparse, tempfile
import numpy as np, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subtel_atlas import readevidence as RE  # noqa: E402
import fetch_reads as FR  # noqa: E402

A = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
inp = os.environ.get('ATLAS_INPUT_DIR', '')
A.add_argument('--chm13', default=os.path.join(inp, 'ref', 'chm13v2.0.fa.gz'))
A.add_argument('--chm13-ends', default=os.path.join(inp, 'subtel', 'foreign', 'ref', 'chm13_ends500k.fa'))
A.add_argument('--genome', default=os.path.join(REPO, 'work', 'readevidence', 'chm13v2.map-ont.mmi'))
A.add_argument('--work', default=os.path.join(REPO, 'work', 'evidence', 'population'))
A.add_argument('--out', default=os.path.join(REPO, 'evidence', 'population', 'diagnostics'))
A.add_argument('-t', '--threads', type=int, default=5)
A.add_argument('--annotate-only', action='store_true', help='only rebuild claim_annotations.tsv from the existing tables')
a = A.parse_args()
os.makedirs(a.out, exist_ok=True)
import pysam, mappy as mp  # noqa: E402

W = a.work
crit = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json')); th = crit['thresholds']
ends = RE.Chm13Ends(a.chm13_ends); wf = pysam.FastaFile(os.path.join(REPO, 'work', 'windows', 'ends.fa'))
MASK = os.path.join(REPO, 'work', 'readevidence', 'masks')
P = json.load(open(os.path.join(W, 'plan.json'))); cm = {c['claim']: c for c in P['claims']}
C = pd.read_csv(os.path.join(REPO, 'evidence', 'population', 'claims.tsv'), sep='\t')
E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False); E2 = E.set_index('end_id')
R2 = E[E['set'] == 'HPRC_R2']
Rd = pd.read_csv(os.path.join(W, 'reads_full.tsv'), sep='\t', low_memory=False); car = Rd[Rd.role == 'carrier']
U = json.load(gzip.open(os.path.join(W, 'units.json.gz'), 'rt'))
asmj = json.load(open(os.path.join(W, 'assembly_projection.json')))
slf = {(x['sample'], x['region']): x['file'] for x in P['slices']}
TMP = tempfile.mkdtemp(prefix='popdiag_')
_ctx = {}


def getctx(cid, preset='map-hifi', model=None):
    k = (cid, preset)
    if k not in _ctx: _ctx[k] = RE.ClaimContext(model or cm[cid]['model'], ends, crit, preset=preset, mask_cache=MASK)
    return _ctx[k]


def proj(ctx, seq):
    r, u, _ = RE.evaluate_read(ctx, 'w', seq.upper(), side='L', meta=dict(platform='assembly'))
    return r, u


def get_reads(sample, region, names, path=None):
    names = set(names)
    return {n: s for n, s, _ in RE.iter_reads(path or os.path.join(W, slf[(sample, region)])) if n in names}


def place(sample, seqs):
    """aligned matching bases of each read per end window of the individual (all Release 2 ends); top 4"""
    ids = [e for e in R2[R2['sample'] == sample].end_id if e in wf.references]
    fa = os.path.join(TMP, f'{sample}.fa')
    with open(fa, 'w') as f:
        for e in ids: f.write(f'>{e}\n{wf.fetch(e)}\n')
    al = mp.Aligner(fa, preset='map-ont', best_n=10); out = {}
    for n, s in seqs.items():
        agg = collections.defaultdict(lambda: [0, 0])
        for h in al.map(s): agg[h.ctg][0] += h.mlen; agg[h.ctg][1] += h.blen
        out[n] = sorted([(c, m, b, round(m / b, 4)) for c, (m, b) in agg.items()], key=lambda x: -x[1])[:4]
    return out


def endid(s, arm, hap):
    return R2[(R2['sample'] == s) & (R2.home_arm == arm) & (R2.haplotype == hap)].end_id.iloc[0]


def lab(e):
    return f"{E2.loc[e, 'home_arm']}:{E2.loc[e, 'structural_allele_id']}:{E2.loc[e, 'haplotype']}"


def annotate():
    cfa = pd.read_csv(os.path.join(a.out, 'confirmed_claims_annotation.tsv'), sep='\t').set_index('claim_id')
    pcp = pd.read_csv(os.path.join(a.out, 'contradicting_read_placement.tsv'), sep='\t').groupby('claim_id')[
        ['on_claimed_end', 'on_arm_other_hap', 'other_end']].sum()
    pla = pd.read_csv(os.path.join(a.out, 'conclusive_read_placement.tsv'), sep='\t')
    off = pla[pla.concordant == False].groupby('claim_id').apply(lambda d: int((~d.placed_on_claimed_arm).sum()))  # noqa: E712
    sc = pd.read_csv(os.path.join(a.out, 'donor_window_edge_screen.tsv'), sep='\t'); edge = set(sc[sc.edge].claim_id)
    REC = {'robust': 'keep', 'shared_anchor_on_arm_placement': 'keep; note anchor in CHM13-shared sequence',
           'not_independently_supported': 'annotate: reads not attributable to the claimed end'}
    rows = []
    for c in C.itertuples():
        notes, rec = [], ''
        if c.claim_id in cfa.index:
            k = cfa.loc[c.claim_id, 'class']; an = cfa.loc[c.claim_id, 'annotation']
            notes.append(f'confirmed:{k}' + (f' ({an})' if isinstance(an, str) and an else '')); rec = REC[k]
        if c.claim_id in edge:
            notes.append('donor-window-edge switch: assembly switch is the end of the v1 D window; switch test does not separate the '
                         'allele from native donor-arm sequence')
            rec = rec or 'treat as not assessable under v1 (model limitation)'
        if c.claim_id in pcp.index:
            p = pcp.loc[c.claim_id]
            notes.append(f'contradicting reads placed: claimed end {int(p.on_claimed_end)}, other haplotype {int(p.on_arm_other_hap)}, '
                         f'other chromosome end {int(p.other_end)}')
            rec = rec or ('no contradiction of the claimed haplotype' if p.on_claimed_end == 0 else
                          'review: contradicting read(s) place on the claimed end')
        if c.assessability == 'discordant' and off.get(c.claim_id, 0) > 0:
            notes.append(f'{off[c.claim_id]} discordant conclusive read(s) place on another chromosome end')
            rec = rec or 'discordant reads not from the claimed end; no evidence against the assembly switch'
        if notes: rows.append(dict(claim_id=c.claim_id, allele=c.allele, verdict=c.verdict, verdict_class=c.verdict_class,
                                   assessability=c.assessability, diagnostic_notes=' | '.join(notes), recommendation=rec))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'claim_annotations.tsv'), sep='\t', index=False)
    print(f'claim_annotations.tsv: {len(rows)} claims', file=sys.stderr)


if a.annotate_only:
    annotate(); sys.exit(0)
t0 = time.time()
# ---------------------------------------------------------------- 1. placement of every conclusive read
conc = car[car.verdict == 'conclusive']; PLA = []
for s, g in conc.groupby('sample'):
    seqs = {}
    for reg, gg in g.groupby('slice_region'): seqs.update(get_reads(s, reg, gg.read))
    pl = place(s, seqs)
    for r in g.itertuples():
        top = pl.get(r.read, []); arm = r.claim_id.split('|')[1]; t1 = top[0] if top else None
        on = [t for t in top if E2.loc[t[0], 'home_arm'] == arm]
        PLA.append(dict(claim_id=r.claim_id, sample=s, read=r.read, concordant=r.switch_concordant, anchor_tier=r.anchor_tier,
                        anchor_unique_bp=r.anchor_unique_bp, best_end=t1[0] if t1 else None,
                        best_arm=E2.loc[t1[0], 'home_arm'] if t1 else None, best_sa=E2.loc[t1[0], 'structural_allele_id'] if t1 else None,
                        best_m=t1[1] if t1 else 0, best_id=t1[3] if t1 else None, arm_m=on[0][1] if on else 0,
                        arm_end=on[0][0] if on else None, placed_on_claimed_arm=bool(t1 and E2.loc[t1[0], 'home_arm'] == arm)))
PLA = pd.DataFrame(PLA); PLA.to_csv(os.path.join(a.out, 'conclusive_read_placement.tsv'), sep='\t', index=False)
cc = PLA[PLA.concordant == True]; g = cc.groupby('claim_id')  # noqa: E712
cf = C[C.verdict == 'confirmed'][['claim_id', 'allele', 'n_conclusive_concordant', 'other_haplotype', 'allele_specific_confirmed']].copy()
cf['n_unique_anchor'] = cf.claim_id.map(g.apply(lambda d: (d.anchor_tier == 'unique').sum())).fillna(0).astype(int)
cf['n_placed_on_arm'] = cf.claim_id.map(g.placed_on_claimed_arm.sum()).fillna(0).astype(int)
cf['n_unique_and_on_arm'] = cf.claim_id.map(g.apply(lambda d: ((d.anchor_tier == 'unique') & d.placed_on_claimed_arm).sum())).fillna(0).astype(int)


def ann(r):
    x = []
    if r.n_placed_on_arm < 2: x.append('reads_place_on_other_end')
    if r.n_unique_anchor < 2: x.append('shared_anchor_only' if r.n_unique_anchor == 0 else 'one_unique_anchor')
    if not r.allele_specific_confirmed: x.append('reads_shared_with_other_allele')
    return ';'.join(x)


cf['annotation'] = cf.apply(ann, axis=1)
cf['class'] = np.where((cf.n_placed_on_arm < 2) | ~cf.allele_specific_confirmed, 'not_independently_supported',
                       np.where(cf.n_unique_and_on_arm >= 2, 'robust', 'shared_anchor_on_arm_placement'))
cf.to_csv(os.path.join(a.out, 'confirmed_claims_annotation.tsv'), sep='\t', index=False)
print(f'[1] placement {time.time() - t0:.0f} s', cf['class'].value_counts().to_dict(), file=sys.stderr)

# ---------------------------------------------------------------- 2. contradicting reads of contradicting_only claims
co = set(C[C.assessability == 'contradicting_only'].claim_id); crd = car[car.claim_id.isin(co) & (car.verdict == 'contradicts')]
Ci = C.set_index('claim_id'); PC = []
for s, g2 in crd.groupby('sample'):
    seqs = {}
    for reg, gg in g2.groupby('slice_region'): seqs.update(get_reads(s, reg, gg.read))
    pl = place(s, seqs)
    for r in g2.itertuples():
        mine = set(Ci.loc[r.claim_id, 'end_ids'].split(',')); arm = Ci.loc[r.claim_id, 'arm']; top = pl.get(r.read, []); t = top[0] if top else None
        PC.append(dict(claim_id=r.claim_id, read=r.read, best_end=t[0] if t else None, best_arm=E2.loc[t[0], 'home_arm'] if t else None,
                       on_claimed_end=bool(t and t[0] in mine), on_arm_other_hap=bool(t and t[0] not in mine and E2.loc[t[0], 'home_arm'] == arm),
                       other_end=bool(t and E2.loc[t[0], 'home_arm'] != arm)))
pd.DataFrame(PC).to_csv(os.path.join(a.out, 'contradicting_read_placement.tsv'), sep='\t', index=False)
print(f'[2] contradicting {time.time() - t0:.0f} s', pd.DataFrame(PC)[['on_claimed_end', 'on_arm_other_hap', 'other_end']].sum().to_dict(), file=sys.stderr)

# ---------------------------------------------------------------- 3. SA-3q-01
dis = car[car.claim_id.str.endswith('SA-3q-01') & (car.verdict == 'conclusive')]
rows = []
for cid in sorted(dis.claim_id.unique()):
    s = cid.split('|')[0]; ctx = getctx(cid)
    for e in R2[(R2['sample'] == s) & (R2.home_arm == '3q')].end_id:
        r, u = proj(ctx, wf.fetch(e))
        rows.append(dict(claim_id=cid, end_id=e, allele=E2.loc[e, 'structural_allele_id'], model_end=e == cm[cid]['model_end'],
                         switch_pass=r.get('switch_pass'), segments=r.get('segments'), switch_ci_sub_lo=r.get('switch_ci_sub_lo'),
                         switch_ci_sub_hi=r.get('switch_ci_sub_hi')))
for e in E[(E['set'] == 'HG002_Q100') & (E.home_arm == '3q')].end_id:
    r, u = proj(getctx('HG00133|3q|SA-3q-01'), wf.fetch(e))
    rows.append(dict(claim_id='HG002 (model of HG00133|3q|SA-3q-01)', end_id=e, allele=E2.loc[e, 'structural_allele_id'], model_end=False,
                     switch_pass=r.get('switch_pass'), segments=r.get('segments'), switch_ci_sub_lo=r.get('switch_ci_sub_lo'),
                     switch_ci_sub_hi=r.get('switch_ci_sub_hi')))
pd.DataFrame(rows).to_csv(os.path.join(a.out, 'sa3q01_window_projections.tsv'), sep='\t', index=False)
pa = PLA[PLA.claim_id.str.endswith('SA-3q-01')].merge(dis[['claim_id', 'read', 'switch_ci_sub_lo', 'switch_ci_sub_hi', 'read_len']], on=['claim_id', 'read'])
out = []
for r in pa.itertuples():
    ctx = getctx(r.claim_id); rb, _ = proj(ctx, wf.fetch(r.best_end)); rc, _ = proj(ctx, wf.fetch(cm[r.claim_id]['model_end']))
    out.append(dict(claim_id=r.claim_id, read=r.read, read_len=r.read_len, read_switch=f'{r.switch_ci_sub_lo:.0f}-{r.switch_ci_sub_hi:.0f}',
                    anchor=f'{r.anchor_tier}/{r.anchor_unique_bp}', best_window=f'{r.best_arm} {r.best_sa} ({E2.loc[r.best_end, "haplotype"]})',
                    best_m=r.best_m, claimed_3q_m=r.arm_m,
                    best_window_switch=f"{rb.get('switch_ci_sub_lo')}-{rb.get('switch_ci_sub_hi')}" if rb.get('switch_pass') else 'no switch',
                    claimed_window_switch=f"{rc.get('switch_ci_sub_lo')}-{rc.get('switch_ci_sub_hi')}",
                    match=bool(rb.get('switch_pass') and RE._ovl(r.switch_ci_sub_lo, r.switch_ci_sub_hi, rb['switch_ci_sub_lo'], rb['switch_ci_sub_hi'],
                                                                 th['concordance_tol_bp']))))
pd.DataFrame(out).to_csv(os.path.join(a.out, 'sa3q01_discordant_reads.tsv'), sep='\t', index=False)
# alternative donor 15q (diagnostic model only)
m0 = dict(cm['HG00133|3q|SA-3q-01']['model'])
h = max(mp.Aligner(seq=ends['chr15q'][:300000], preset='asm20').map(ends['chr11p'][m0['d_lo']:m0['d_hi']]), key=lambda x: x.mlen)
m15 = dict(m0, D='chr15q', d_lo=int(h.r_st), d_hi=int(h.r_en), claim_id='alt_D_15q')
c15h = getctx('alt15h', 'map-hifi', m15); c15o = getctx('alt15o', 'map-ont', m15); alt = []
for e in [cm['HG00133|3q|SA-3q-01']['model_end'], endid('HG00133', '6p', 'hap2'), 'HG002#2#chr3_MATERNAL#R']:
    r, _ = proj(c15h, wf.fetch(e)); alt.append(dict(item=lab(e) if e in E2.index else e, kind='window', switch_pass=r.get('switch_pass'),
                                                  segments=r.get('segments'), lo=r.get('switch_ci_sub_lo'), hi=r.get('switch_ci_sub_hi')))
for (s, reg), gg in dis.groupby(['sample', 'slice_region']):
    side = 'L' if gg.slice_end.iloc[0].endswith('p') else 'R'
    for n, sq in get_reads(s, reg, gg.read).items():
        r, _, _ = RE.evaluate_read(c15o, n, sq, side=side, meta=dict(platform='ONT_R9.4.1'))
        alt.append(dict(item=f'{s} {n}', kind='read', switch_pass=r.get('switch_pass'), segments=r.get('segments'), lo=r.get('switch_ci_sub_lo'),
                        hi=r.get('switch_ci_sub_hi')))
pd.DataFrame(alt).to_csv(os.path.join(a.out, 'sa3q01_alt_donor_15q.tsv'), sep='\t', index=False)
print(f'[3] SA-3q-01 {time.time() - t0:.0f} s', file=sys.stderr)

# ---------------------------------------------------------------- 4. D-window-edge screen
scr = []
for c in P['claims']:
    if not c.get('model_valid'): continue
    m = c['model']; pj = asmj.get(f"{c['claim']}||{c['model_end']}") or {}
    if not pj.get('switch_pass') or pj.get('switch_ci_P_lo') is None: continue
    dw, pw = RE.claim_windows(m, th); bl = RE.PairMap(ends, m['D'], m['P'], dw, pw, th).summary()['blocks']
    if not bl: continue
    last = max(bl, key=lambda b: b['d_end']); off = last['p_start'] - last['d_start']; dend = dw[1]
    qo = max(0, dend - 2000); to = max(0, last['p_end'] - 10000); ext = 0
    for hh in mp.Aligner(seq=ends[m['P']][to:last['p_end'] + 40000], preset='asm20').map(ends[m['D']][qo:dend + 20000]):
        if hh.strand == 1 and hh.mlen / hh.blen >= 0.85 and abs((hh.r_st + to) - (hh.q_st + qo) - off) <= 5000: ext = max(ext, hh.q_en + qo - dend)
    sp = (pj['switch_ci_P_lo'] + pj['switch_ci_P_hi']) / 2
    scr.append(dict(claim_id=c['claim'], allele=c['allele'], D=m['D'], P=m['P'], d_win_end=dend, block_p_end=last['p_end'], asm_switch_P=sp,
                    homology_beyond_D_window=ext, edge=bool(abs(sp - last['p_end']) <= 3000 and ext >= 5000)))
SC = pd.DataFrame(scr).merge(C[['claim_id', 'verdict', 'verdict_class', 'assessability']], on='claim_id')
SC.to_csv(os.path.join(a.out, 'donor_window_edge_screen.tsv'), sep='\t', index=False)
print(f'[4] edge screen {time.time() - t0:.0f} s: {int(SC.edge.sum())} of {len(SC)}', file=sys.stderr)

# ---------------------------------------------------------------- 5. SA-20p-03 wider chr20p fetch and donor-slice reads
D20 = ['HG00235', 'HG00146', 'NA18960', 'HG02155', 'NA18570', 'NA18940', 'HG01530', 'NA21102']
IDX = os.path.join(REPO, 'work', 'readevidence', 'index'); jobs = []
for s in D20:
    url = FR.ONT1KG.format(s=s); out = os.path.join(W, 'diag', 'slices', s, 'chr20p_0_60000.bam')
    FR.fetch(url, FR.index_for(url, IDX, [IDX]), 'chr20:1-60000', out, a.chm13)
    cid = [c for c in C[C.allele == 'SA-20p-03'].claim_id if c.startswith(s + '|')][0]; jobs.append((cid, out))
allrows, Q = [], []
for cid, bam in jobs:
    ctx = getctx(cid, 'map-ont'); reads = list(RE.iter_reads(bam, min_len=5000))
    rows, units, q = RE.stage1(ctx, reads, side='L', dataset='wide60k', platform='ONT_R9.4.1', tag=f'{cid}|'); allrows.append((cid, rows, units)); Q += q
G = RE.map_genome(Q, a.genome, a.threads); W20 = []
for cid, rows, units in allrows:
    ctx = getctx(cid, 'map-ont'); RE.stage2(ctx, rows, G)
    for r in rows:
        RE.concordance(r, asmj.get(f"{cid}||{cm[cid]['model_end']}"), th); r['claim'] = cid; r['_units'] = units[r['read']]; W20.append(r)
pd.DataFrame([{k: r.get(k) for k in ['claim', 'read', 'read_len', 'terminal', 'segments', 'switch_pass', 'donor_interval_units',
                                    'donor_interval_P_purity', 'anchor_pass', 'anchor_tier', 'verdict']} for r in W20]).to_csv(
    os.path.join(a.out, 'sa20p03_wide60kb_reads.tsv'), sep='\t', index=False)
RC = []
for s in ['HG00235', 'HG00146', 'NA18960']:
    cid = [c for c in C[C.allele == 'SA-20p-03'].claim_id if c.startswith(s + '|')][0]
    g3 = car[(car.claim_id == cid) & (car.slice_end == 'chr7p') & (car.verdict == 'consistent')]; seqs = {}
    for reg, gg in g3.groupby('slice_region'): seqs.update(get_reads(s, reg, gg.read))
    pl = place(s, seqs)
    for r in g3.itertuples():
        top = pl.get(r.read, []); t = top[0] if top else None
        RC.append(dict(sample=s, claim_id=cid, read=r.read, read_len=r.read_len, switch=f'{r.switch_ci_sub_lo:.0f}-{r.switch_ci_sub_hi:.0f}',
                       anchor_chrom=r.anchor_chrom, anchor_mapq=r.anchor_mapq, best=lab(t[0]) if t else None, best_m=t[1] if t else 0))
pd.DataFrame(RC).to_csv(os.path.join(a.out, 'sa20p03_donor_slice_reads_placement.tsv'), sep='\t', index=False)
print(f'[5] SA-20p-03 {time.time() - t0:.0f} s', file=sys.stderr)

# ---------------------------------------------------------------- 6. figures
import matplotlib as mpl  # noqa: E402
mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
mpl.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'D': '#d95f02', 'P': '#1f78b4', 'o': '#9e9e9e'}


def read_units(row):
    return [(p - row.s0_read, al[0]) for p, al, k, dp, pp in U.get(f'{row.job}|{row.read}', [])]


def win_units(cid, e):
    r, u = proj(getctx(cid), wf.fetch(e)); s0 = r['s0_read']
    return [(v['read_pos'] - s0, v['allele'][0]) for v in u], (r.get('switch_ci_sub_lo'), r.get('switch_ci_sub_hi')) if r.get('switch_pass') else None


def draw(tracks, xmax, dlabel, plabel, fn, title, width=7.2):
    h = 0.21 * len(tracks) + 1.2; fig, ax = plt.subplots(figsize=(width, h)); y = 0; yt, yl = [], []
    for kind, lb, un, ci in tracks:
        if kind == 'header':
            y += 0.3; ax.text(0, -y + 0.05, lb, ha='left', va='bottom', fontsize=7, fontweight='bold',
                              transform=mpl.transforms.blended_transform_factory(fig.transFigure, ax.transData)); y += 0.7; continue
        if ci and ci[0] is not None:
            if kind == 'window': ax.add_patch(mpl.patches.Rectangle((ci[0], -y - 0.38), max(ci[1] - ci[0], 250), 0.76, color='#e0e0e0', zorder=0))
            else: ax.plot([(ci[0] + ci[1]) / 2] * 2, [-y - 0.42, -y + 0.42], color='k', lw=1.2, zorder=3)
        for al in ('o', 'P', 'D'):
            xs = [p for p, b in un if b == al and 0 <= p <= xmax]; ax.scatter(xs, [-y] * len(xs), marker='|', s=26, color=COL[al], lw=0.8, zorder=2)
        yt.append(-y); yl.append(lb); y += 1
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_xlim(0, xmax); ax.set_ylim(-y - 0.1, 0.3)
    ax.set_xlabel('distance from the telomere boundary (kb)'); ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f'{v / 1000:.0f}'))
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    hs = [mpl.lines.Line2D([], [], marker='|', ls='', color=COL['D'], ms=8, mew=1.5, label=dlabel),
          mpl.lines.Line2D([], [], marker='|', ls='', color=COL['P'], ms=8, mew=1.5, label=plabel),
          mpl.lines.Line2D([], [], marker='|', ls='', color=COL['o'], ms=8, mew=1.5, label='neither'),
          mpl.patches.Patch(color='#e0e0e0', label='assembly switch, 95% CI'), mpl.lines.Line2D([], [], color='k', lw=1.2, label='read switch')]
    fig.tight_layout(rect=(0, 0.06, 1, 0.97)); fig.legend(handles=hs, loc='lower center', ncol=5, frameon=False, handletextpad=0.3, columnspacing=1.0)
    fig.text(0.01, 0.995, title, ha='left', va='top', fontsize=8); fig.savefig(fn + '.png', dpi=300); fig.savefig(fn + '.pdf'); plt.close(fig)


T = []


def group(title, cid, windows, reads):
    T.append(('header', title, None, None))
    for lb, e in windows: un, ci = win_units(cid, e); T.append(('window', lb, un, ci))
    for row in reads.itertuples(): T.append(('read', f'read {row.read[:6]} ({row.read_len / 1000:.0f} kb)', read_units(row), (row.switch_ci_sub_lo, row.switch_ci_sub_hi)))


group('HG00133 (3q: SA-3q-01 / SA-3q-02 native)', 'HG00133|3q|SA-3q-01', [('3q hap1, SA-3q-01 (claimed)', endid('HG00133', '3q', 'hap1')),
      ('3q hap2, SA-3q-02', endid('HG00133', '3q', 'hap2')), ('6p hap2, SA-6p-118', endid('HG00133', '6p', 'hap2'))], dis[dis['sample'] == 'HG00133'])
group('HG00099 (3q: SA-3q-01 / SA-3q-39)', 'HG00099|3q|SA-3q-01', [('3q hap2, SA-3q-01 (claimed)', endid('HG00099', '3q', 'hap2')),
      ('3q hap1, SA-3q-39', endid('HG00099', '3q', 'hap1')), ('6p hap2, SA-6p-01', endid('HG00099', '6p', 'hap2'))], dis[dis['sample'] == 'HG00099'])
group('NA19043 (3q: SA-3q-01 / SA-3q-02 native)', 'NA19043|3q|SA-3q-01', [('3q hap2, SA-3q-01 (claimed)', endid('NA19043', '3q', 'hap2')),
      ('20q hap1, SA-20q-01', endid('NA19043', '20q', 'hap1'))], dis[dis['sample'] == 'NA19043'])
group('HG002 Q100 (reference for the claim model)', 'HG00133|3q|SA-3q-01', [('3q maternal, SA-3q-01', 'HG002#2#chr3_MATERNAL#R'),
      ('3q paternal, SA-3q-04', 'HG002#1#chr3_PATERNAL#R')], dis.iloc[0:0])
draw(T, 50000, 'CHM13 11p allele (donor)', 'CHM13 3q allele (recipient)', os.path.join(a.out, 'fig_sa3q01_tracks'),
     'SA-3q-01 "discordant" reads switch where the person\'s 6p or 20q window switches, not where the claimed 3q window does')
T = []; cB = 'HG00235|20p|SA-20p-03'; RCd = pd.DataFrame(RC)
T.append(('header', 'HG00235 (20p hap1: SA-20p-03; 20p hap2 not assembled)', None, None))
for lb, e in [('20p hap1, SA-20p-03 (claimed)', endid('HG00235', '20p', 'hap1')), ('7p hap1, SA-7p-01 (native 7p)', endid('HG00235', '7p', 'hap1')),
              ('7p hap2, SA-7p-01 (native 7p)', endid('HG00235', '7p', 'hap2')), ('20q hap2, SA-20q-06 (20p-like block)', endid('HG00235', '20q', 'hap2'))]:
    un, ci = win_units(cB, e); T.append(('window', lb, un, ci))
dsl = car[(car.claim_id == cB) & (car.slice_end == 'chr7p') & (car.verdict == 'consistent')]
for best, tag in [('20p:SA-20p-03:hap1', 'placed on SA-20p-03'), ('7p:SA-7p-01:hap1', 'placed on 7p')]:
    for rr in RCd[(RCd['sample'] == 'HG00235') & (RCd.best == best)].head(2).itertuples():
        row = dsl[dsl.read == rr.read].iloc[0]
        T.append(('read', f'chr7p-slice read, {tag} ({row.read_len / 1000:.0f} kb)', read_units(row), (row.switch_ci_sub_lo, row.switch_ci_sub_hi)))
wu = lambda r: [(v['read_pos'] - r['s0_read'], v['allele'][0]) for v in r['_units']]  # noqa: E731
for r in [r for r in W20 if r['claim'] == cB and r['verdict'] == 'contradicts'][:3]: T.append(('read', f"chr20p read, contradicts ({r['read_len'] / 1000:.0f} kb)", wu(r), None))
cB2 = 'HG00146|20p|SA-20p-03'; T.append(('header', 'HG00146 (20p hap1: SA-20p-01 native; hap2: SA-20p-03)', None, None))
for lb, e in [('20p hap2, SA-20p-03 (claimed)', endid('HG00146', '20p', 'hap2')), ('20p hap1, SA-20p-01 (native)', endid('HG00146', '20p', 'hap1')),
              ('9p hap1, SA-9p-04 (20p-like block)', endid('HG00146', '9p', 'hap1'))]:
    un, ci = win_units(cB2, e); T.append(('window', lb, un, ci))
for r in [r for r in W20 if r['claim'] == cB2 and r['verdict'] == 'contradicts'][:3]: T.append(('read', f"chr20p read, contradicts ({r['read_len'] / 1000:.0f} kb)", wu(r), None))
draw(T, 60000, 'CHM13 7p allele (donor)', 'CHM13 20p allele (recipient)', os.path.join(a.out, 'fig_sa20p03_tracks'),
     'SA-20p-03: native 7p haplotypes show the same 7p-to-20p allele switch as the claimed 20p window')
annotate()
print(f'done {time.time() - t0:.0f} s', file=sys.stderr)
