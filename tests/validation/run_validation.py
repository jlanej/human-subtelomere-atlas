"""Validation of the subtel-atlas query tool against catalogue v0.1 (results in tests/validation/, summary in VALIDATION.md).

Stages (python tests/validation/run_validation.py <stage>; each stage is resumable and writes repo-relative outputs):
  select   stratified query set: every end of 60 Release 2 samples (12 per superpopulation; seed 20261015) + per arm and tier up to
           3 ends of structural alleles and 2 of haplogroups of each tier + 40 ends with a telomere block reaching 20 kb + all
           CHM13 and HG002 Q100 ends -> tests/validation/query_set.tsv
  query    re-query the selected catalogue windows (work/windows/ends.fa; telomere re-called by the tool; arm supplied from the catalogue
           so that mosaic/allele/haplogroup steps are compared at the same arm; the tool's own mosaic-rule arm is also recorded)
           -> work/qtool/val/chunks/*.json
  compare  (a) self-consistency vs ends.tsv -> self_consistency.tsv, self_consistency_summary.tsv
  loso     (b) leave-one-sample-out for the 60 samples (+ CHM13, HG002) -> loso.tsv, loso_summary.tsv
  indep    (c) CHM13 v2.0 and HG002 Q100 whole assemblies and the WashU pedigree ends in contig mode; HG002 Q100 vs HG002 Release 2;
           pedigree transmissions -> independent_*.tsv
  probes   arm probe rule (tool implementation, mappy asm5) on catalogue probe sequences -> probe_rule.tsv
  reads    experimental read mode on HG002 HiFi / ONT UL reads at chromosome ends -> reads_hg002.tsv
  mappy_backend  mappy backend vs minimap2 backend on 400 ends -> backend_mappy.tsv
  shift    boundary-shift robustness on 600 stratified ends, with re-alignment: (A) boundary moved -500 / +500 bp (header override);
           (B) telomere lengthened (+250 / +500 bp of (CCCTAA)n inserted) or shortened (-250 / -500 bp trimmed, window kept at
           100 kb by N-padding its inner end), boundary re-called -> boundary_shift.tsv, boundary_shift_summary.tsv
Inputs: work/windows (catalogue build), work/qtool/refs (subtel-atlas fetch-refs --dir work/qtool/refs), $ATLAS_INPUT_DIR (pedigree
ends, HG002 read slices, CHM13 v2.0 FASTA; see pipeline/config.sh), work/qtool/hg002/hg002v1.1.fasta.gz (public HG002 Q100 v1.1).
"""
import os, sys, json, glob, time, gzip, subprocess, collections
import numpy as np, pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import query as Q, alleles as AL, haplo as H, sketch as S, arms as AR
from subtel_atlas.common import group_of, readfa

OUT = os.path.join(REPO, 'tests', 'validation')
VW = os.path.join(REPO, 'work', 'qtool', 'val_r2'); os.makedirs(VW + '/chunks', exist_ok=True)
REFS = os.environ.get('SUBTEL_ATLAS_REFS', os.path.join(REPO, 'work', 'qtool', 'refs'))
WIN = os.path.join(REPO, 'work', 'windows')
INPUT = os.environ.get('ATLAS_INPUT_DIR', '')
THREADS = int(os.environ.get('ATLAS_THREADS', '3'))
SEED = 20261015
rd = lambda f, **k: pd.read_csv(f, sep='\t', low_memory=False, **k)
E = rd(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'))
EI = E.set_index('end_id')
SAt = rd(os.path.join(REPO, 'data', 'catalogue', 'alleles_structural.tsv')); HGt = rd(os.path.join(REPO, 'data', 'catalogue', 'haplogroups.tsv'))


def fetch_windows(ids):
    import pysam
    fa = pysam.FastaFile(os.path.join(WIN, 'ends.fa'))
    return {e: fa.fetch(e) for e in ids}


# ------------------------------------------------------------------------------------------------ select
def select():
    rng = np.random.default_rng(SEED)
    r2 = E[E.set == 'HPRC_R2']; samp = []
    for sp, g in sorted(r2.groupby('superpopulation')):
        s = sorted(g['sample'].unique()); s = [x for x in s if x != 'HG002']
        samp += list(rng.choice(s, 12, replace=False))
    rows = [dict(end_id=e, reason='loso_sample') for e in r2[r2['sample'].isin(samp)].end_id]
    have = {r['end_id'] for r in rows}
    sa_t = dict(zip(SAt.allele_id, SAt.tier)); hg_t = dict(zip(HGt.haplogroup_id, HGt.tier))
    pool = E[(E.set != 'pedigree_aux') & ~E.end_id.isin(have)].copy()
    pool['sa_tier'] = pool.structural_allele_id.map(sa_t); pool['hg_tier'] = pool.haplogroup_id.map(hg_t)
    for (arm, t), g in sorted(pool.dropna(subset=['sa_tier']).groupby(['catalogue_arm', 'sa_tier'])):
        for e in g.sample(min(3, len(g)), random_state=SEED).end_id:
            if e not in have: rows.append(dict(end_id=e, reason=f'sa_tier_{t}')); have.add(e)
    for (arm, t), g in sorted(pool.dropna(subset=['hg_tier']).groupby(['catalogue_arm', 'hg_tier'])):
        g = g[~g.end_id.isin(have)]
        for e in g.sample(min(2, len(g)), random_state=SEED).end_id: rows.append(dict(end_id=e, reason=f'hg_tier_{t}')); have.add(e)
    cap = pool[(pool.junk_len + pool.telo_len >= 20000) & ~pool.end_id.isin(have)]
    for e in cap.sample(min(40, len(cap)), random_state=SEED).end_id: rows.append(dict(end_id=e, reason='telomere_20kb_cap')); have.add(e)
    for e in E[E.set.isin(['CHM13', 'HG002_Q100'])].end_id:
        if e not in have: rows.append(dict(end_id=e, reason='t2t')); have.add(e)
    D = pd.DataFrame(rows).merge(E[['end_id', 'set', 'sample', 'superpopulation', 'catalogue_arm']], on='end_id')
    D.to_csv(os.path.join(OUT, 'query_set.tsv'), sep='\t', index=False)
    print('query set', len(D), D.reason.value_counts().to_dict(), 'samples', len(samp))


# ------------------------------------------------------------------------------------------------ query
def query(chunk=400):
    D = rd(os.path.join(OUT, 'query_set.tsv')); ref = Q.Reference(REFS)
    ids = list(D.end_id); chunk = int(chunk)
    for k in range(0, len(ids), chunk):
        fn = f'{VW}/chunks/c{k // chunk:03d}.json'
        if os.path.exists(fn): continue
        part = ids[k:k + chunk]; seqs = fetch_windows(part); ends = []
        for e in part:
            r = EI.loc[e]; x = Q.window_end(e, f'arm={r.home_arm}' if r.home_arm != 'unassigned' else '', seqs[e], name_arms=False)
            x['arm_hint_source'] = 'catalogue_' + r.arm_source if x['arm_hint'] else ''
            ends.append(x)
        t = time.time(); recs = Q.analyse(ends, ref, 'minimap2', THREADS)
        json.dump(dict(seconds=time.time() - t, n=len(part), records=Q._clean(recs)), open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
        print(fn, len(part), round(time.time() - t, 1), 's', flush=True)


def load_chunks(pattern='chunks/c*.json'):
    recs, sec = [], 0.0
    for fn in sorted(glob.glob(f'{VW}/{pattern}')):
        d = json.load(open(fn)); recs += d['records']; sec += d['seconds']
    return recs, sec


# ------------------------------------------------------------------------------------------------ (a) self-consistency
nz = lambda v: '' if v is None or (isinstance(v, float) and v != v) else str(v)
REPS_HG = set(HGt[HGt.representative_in_fasta].haplogroup_id)
TIER_SA = dict(zip(SAt.allele_id, SAt.tier)); NAT_SA = dict(zip(SAt.allele_id, SAt.nativity_class if 'nativity_class' in SAt else np.where(SAt.is_native, 'native', 'non_native'))); TIER_HG = dict(zip(HGt.haplogroup_id, HGt.tier))


def compare_rows(recs, idmap=None, cols_prefix=''):
    """Per-end comparison of tool records with ends.tsv rows (idmap: tool query_id -> catalogue end_id)."""
    out = []
    for r in recs:
        eid = (idmap or {}).get(r['query_id'], r['query_id'])
        if eid not in EI.index: continue
        e = EI.loc[eid]; s0c = int(e.junk_len + e.telo_len)
        sa_c, hg_c = nz(e.structural_allele_id), nz(e.haplogroup_id)
        o = dict(end_id=eid, set=e['set'], sample=e['sample'], catalogue_arm=e.catalogue_arm, arm_source=e.arm_source,
                 s0_catalogue=s0c, s0_tool=r['s0'], s0_match=s0c == r['s0'], arm_tool=r['catalogue_arm'], arm_match=r['catalogue_arm'] == e.catalogue_arm,
                 mosaic_arm_tool=r['mosaic_arm'], mosaic_arm_catalogue=nz(e.mosaic_arm), mosaic_arm_match=nz(r['mosaic_arm']) == nz(e.mosaic_arm),
                 mosaic_group_true=group_of(r['mosaic_arm']) == e.catalogue_arm if e.arm_source == 'contig_name' and r['mosaic_arm'] else None)
        for t in ('S5', 'S10', 'S20'):        # state and, for foreign windows, the paralog family
            o[f'{t}_match'] = nz(r[f'{t}_state']) == nz(e[f'{t}_state']) and nz(r.get(f'{t}_family')) == nz(e.get(f'{t}_family'))
        o['S20_tool'], o['S20_catalogue'] = r['S20_state'], e.S20_state
        o['S20_family_tool'], o['S20_family_catalogue'] = nz(r.get('S20_family')), nz(e.get('S20_family'))
        o['window_len_match'] = int(r['window_len']) == int(e.window_len)
        o['form_match'] = nz(r['form']) == nz(e.form)
        o['boundary_stable_allele_match'] = nz(r.get('boundary_stable_allele')) == nz(e.get('boundary_stable_allele'))
        o['shift_alleles_match'] = nz(r.get('allele_minus500')) == nz(e.get('allele_minus500')) and nz(r.get('allele_plus500')) == nz(e.get('allele_plus500'))
        o['boundary_stable_allele_catalogue'] = nz(e.get('boundary_stable_allele'))
        o['profile_match'] = nz(r['profile']) == nz(e.profile)
        o['core_structure_match'] = nz(r['core_structure']) == nz(e.core_structure)
        o['core_switches_match'] = nz(r['core_switches_sub']) == nz(e.core_switches_sub).replace('.0', '')
        o.update(sa_catalogue=sa_c, sa_tool=nz(r['structural_allele_id']), sa_match=sa_c == nz(r['structural_allele_id']),
                 sa_assignment_catalogue=nz(e.sa_assignment), sa_assignment_tool=r['sa_assignment'], sa_tier=TIER_SA.get(sa_c, ''),
                 sa_nativity=NAT_SA.get(sa_c, ''),
                 sa_candidates_tool=nz(r['sa_candidates']), novelty_flags=nz(r['novelty_flags']))
        o.update(hg_catalogue=hg_c, hg_tool=nz(r['haplogroup_id']), hg_nearest_tool=nz(r.get('hg_nearest_id')), hg_match=hg_c == nz(r['haplogroup_id']),
                 hg_catalogue_represented=hg_c in REPS_HG, hg_tier=TIER_HG.get(hg_c, ''), hg_ani=r['hg_ani'], hg_status=r['hg_status'],
                 hg_assignment_catalogue=nz(e.hg_assignment), ani_to_medoid_catalogue=e.ani_to_medoid, hg_aln_ident=r['hg_aln_ident'], hg_aln_frac=r['hg_aln_frac'])
        out.append(o)
    return pd.DataFrame(out)


def rate(x): x = pd.Series(x).dropna(); return (int(x.sum()), int(len(x)), round(float(x.mean()), 4) if len(x) else np.nan)


def summarise_cmp(C, label):
    rows = []
    def add(stratum, sub):
        for m in ['s0_match', 'window_len_match', 'arm_match', 'mosaic_arm_match', 'mosaic_group_true', 'S5_match', 'S10_match', 'S20_match',
                  'form_match', 'profile_match', 'core_structure_match', 'core_switches_match', 'sa_match', 'hg_match',
                  'boundary_stable_allele_match', 'shift_alleles_match']:
            k, n, f = rate(sub[m].astype('boolean') if m in sub else [])
            rows.append(dict(analysis=label, stratum=stratum, metric=m, n_agree=k, n=n, agreement=f))
    add('all', C)
    add('boundary_identical', C[C.s0_match])
    for t in ['A', 'B', 'C']: add(f'sa_tier_{t}', C[C.s0_match & (C.sa_tier == t)])
    for t in ['A', 'B', 'C']: add(f'hg_tier_{t}', C[C.s0_match & (C.hg_tier == t)])
    for t in sorted(set(C.sa_nativity) - {''}): add(f'sa_nativity={t}', C[C.s0_match & (C.sa_nativity == t)])
    if 'boundary_stable_allele_catalogue' in C:
        add('catalogue_boundary_stable_allele', C[C.s0_match & (C.boundary_stable_allele_catalogue == 'True')])
        add('catalogue_boundary_unstable_allele', C[C.s0_match & (C.boundary_stable_allele_catalogue == 'False')])
    for a, g in C[C.s0_match].groupby('sa_assignment_catalogue'): add(f'sa_assignment={a}', g)
    add('hg_catalogue_represented', C[C.s0_match & C.hg_catalogue_represented])
    add('hg_catalogue_not_represented', C[C.s0_match & ~C.hg_catalogue_represented & (C.hg_catalogue != '')])
    for st, g in C.groupby('set'): add(f'set={st}', g)
    return pd.DataFrame(rows)


def compare():
    recs, sec = load_chunks(); C = compare_rows(recs)
    D = rd(os.path.join(OUT, 'query_set.tsv'))[['end_id', 'reason', 'superpopulation']]
    C = C.merge(D, on='end_id', how='left')
    C.to_csv(os.path.join(OUT, 'self_consistency.tsv'), sep='\t', index=False)
    Sm = summarise_cmp(C, 'self_consistency'); Sm.to_csv(os.path.join(OUT, 'self_consistency_summary.tsv'), sep='\t', index=False)
    a = Sm[Sm.stratum.isin(['all', 'boundary_identical'])]
    print(a.to_string()); print('query seconds', round(sec, 1), 'ends', len(recs), 's/end', round(sec / max(len(recs), 1), 3))


# ------------------------------------------------------------------------------------------------ (b) leave-one-sample-out
def _jacc(sk):
    import scipy.sparse as sp
    allh = np.unique(np.concatenate(sk)); rows = np.concatenate([np.full(len(s), i) for i, s in enumerate(sk)])
    cols = np.concatenate([np.searchsorted(allh, s) for s in sk])
    X = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(sk), len(allh)))
    I = (X @ X.T).toarray().astype(np.float64); sz = np.array([len(s) for s in sk], float)
    return I / np.maximum(sz[:, None] + sz[None, :] - I, 1)


def loso_reference(ref, sample, sk_cache, win_cache):
    """Reference with `sample`'s ends removed: structural alleles re-clustered from the remaining defining ends (IDs mapped to the
    original allele by majority of members); haplogroup medoids from `sample` replaced by the remaining member with the highest mean
    ANI to the other remaining clustered members; haplogroups seen in < 3 remaining samples dropped from the representatives."""
    import copy
    R = copy.copy(ref); Dd = E[(E.sa_assignment == 'full') & (E['sample'] != sample)]
    alle = {}
    for (arm, st), g in Dd.groupby(['catalogue_arm', 'core_structure'], sort=True):
        sws = [[int(float(x)) for x in str(v).split(',')] if isinstance(v, str) and v else [] for v in g.core_switches_sub]
        for c in AL.cluster_structure(list(g.end_id), sws):
            orig = g.set_index('end_id').loc[c['members'], 'structural_allele_id'].value_counts()
            aid = sorted(orig.index[orig == orig.max()])[0]
            alle.setdefault(arm, []).append(dict(allele_id=aid, structure=st, sw=c['sw'], n=len(c['members'])))
    R.alleles = alle
    core = E[(E.set != 'pedigree_aux') & E.haplogroup_id.notna() & (E['sample'] != sample)]
    ns = core.groupby('haplogroup_id')['sample'].nunique()
    med = {}
    for hid, h in HGt.set_index('haplogroup_id').iterrows():
        if ns.get(hid, 0) < 3: continue
        mend = h.medoid_end_id
        if EI.at[mend, 'sample'] == sample:
            mem = sorted(core[(core.haplogroup_id == hid) & (core.hg_assignment == 'full')].end_id)
            if not mem: continue
            need = [m for m in mem if m not in sk_cache]
            if need:
                for m, w in fetch_windows(need).items():
                    s0 = int(EI.at[m, 'junk_len'] + EI.at[m, 'telo_len'])
                    win_cache[m] = w[s0:s0 + H.HG_LEN]; sk_cache[m] = S.sketch(win_cache[m], H.K, H.SCALE)
            if len(mem) == 1: mend = mem[0]
            else:
                A = S.ani_from_jaccard(_jacc([sk_cache[m] for m in mem]), H.K); np.fill_diagonal(A, 1.0)
                mean = (A.sum(1) - 1) / (len(mem) - 1); mend = sorted(range(len(mem)), key=lambda j: (-mean[j], mem[j]))[0]; mend = mem[mend]
        if mend in ref.rep:
            w, s0 = ref.rep[mend]; sub = w[s0:s0 + H.HG_LEN]
        else:
            if mend not in win_cache:
                w = fetch_windows([mend])[mend]; s0 = int(EI.at[mend, 'junk_len'] + EI.at[mend, 'telo_len'])
                win_cache[mend] = w[s0:s0 + H.HG_LEN]
            sub = win_cache[mend]
        med.setdefault(h.arm, []).append((hid, sub, mend))
    R.medoids = med
    return R


def loso_eval(ref_l, recs, sample, seqs):
    ex_sa = {a['allele_id'] for v in ref_l.alleles.values() for a in v}
    ex_hg = {m[0] for v in ref_l.medoids.values() for m in v}
    out = []
    for r in recs:
        e = EI.loc[r['end_id_cat']]
        sw = [int(x) for x in str(r['core_switches_sub']).split(',')] if r['core_switches_sub'] else []
        sub = seqs[r['end_id_cat']][r['s0']:r['s0'] + H.HG_LEN]
        a = Q.assign(ref_l, r['catalogue_arm'], r['core_structure'] or '', sw, r['sub_len'], bool(r['full_length']), sub)
        sa_t = a['sa']['allele_id'] if a['sa'] else ''
        hg_t = a['hg']['id'] if a['hg'] and a['hg']['ani'] >= H.HG_ANI else ''
        sa_c, hg_c = nz(e.structural_allele_id), nz(e.haplogroup_id)
        exp_sa = sa_c if sa_c in ex_sa else ''; exp_hg = hg_c if hg_c in ex_hg else ''
        out.append(dict(held_out_sample=sample, end_id=r['end_id_cat'], set=e['set'], catalogue_arm=e.catalogue_arm, sa_catalogue=sa_c,
                        sa_tier=TIER_SA.get(sa_c, ''), sa_assignment_catalogue=nz(e.sa_assignment), sa_still_defined=sa_c in ex_sa,
                        sa_expected=exp_sa, sa_full_ref=nz(r['structural_allele_id']), sa_loso=sa_t, sa_loso_correct=sa_t == exp_sa,
                        hg_catalogue=hg_c, hg_tier=TIER_HG.get(hg_c, ''), hg_still_represented=hg_c in ex_hg, hg_expected=exp_hg,
                        hg_full_ref=nz(r['haplogroup_id']), hg_loso=hg_t, hg_loso_ani=round(a['hg']['ani'], 5) if a['hg'] else np.nan,
                        hg_loso_correct=hg_t == exp_hg, s0_match=r['s0'] == int(e.junk_len + e.telo_len)))
    return out


def loso_outcome(exp, got):
    if exp and got == exp: return 'correct_assignment'
    if exp and not got: return 'false_novel'
    if exp and got != exp: return 'wrong_allele'
    if not exp and not got: return 'correct_novel'
    return 'assigned_to_other'


def loso(extra=''):
    recs, _ = load_chunks(); D = rd(os.path.join(OUT, 'query_set.tsv'))
    samples = sorted(D[D.reason == 'loso_sample']['sample'].unique())
    by_s = collections.defaultdict(list)
    for r in recs:
        r['end_id_cat'] = r['query_id']; by_s[EI.at[r['query_id'], 'sample']].append(r)
    ref = Q.Reference(REFS); skc, wc = {}, {}; out = []
    fn = os.path.join(OUT, 'loso.tsv')
    for i, smp in enumerate(samples):
        t = time.time(); Rl = loso_reference(ref, smp, skc, wc)
        rs = [r for r in by_s[smp] if EI.at[r['query_id'], 'set'] == 'HPRC_R2']
        seqs = fetch_windows([r['query_id'] for r in rs])
        out += loso_eval(Rl, rs, smp, seqs)
        print(i, smp, len(rs), round(time.time() - t, 1), 's', flush=True)
    L = pd.DataFrame(out)
    L['sa_outcome'] = [loso_outcome(a, b) for a, b in zip(L.sa_expected, L.sa_loso)]
    L['hg_outcome'] = [loso_outcome(a, b) for a, b in zip(L.hg_expected, L.hg_loso)]
    L.to_csv(fn, sep='\t', index=False); loso_summary(L)


def loso_summary(L=None):
    L = L if L is not None else rd(os.path.join(OUT, 'loso.tsv'), keep_default_na=False)
    rows = []
    for kind in ('sa', 'hg'):
        for strat, g in [('all', L)] + [(f'tier_{t}', L[L[f'{kind}_tier'] == t]) for t in 'ABC'] + \
                        [('still_defined' if kind == 'sa' else 'still_represented', L[L[f'{kind}_still_defined' if kind == 'sa' else f'{kind}_still_represented'].astype(str) == 'True']),
                         ('removed_by_holdout', L[(L[f'{kind}_catalogue'] != '') & (L[f'{kind}_still_defined' if kind == 'sa' else f'{kind}_still_represented'].astype(str) != 'True')])]:
            vc = g[f'{kind}_outcome'].value_counts()
            ok = int(vc.get('correct_assignment', 0) + vc.get('correct_novel', 0))
            rows.append(dict(allele_type='structural_allele' if kind == 'sa' else 'haplogroup', stratum=strat, n=len(g), n_correct=ok,
                             accuracy=round(ok / len(g), 4) if len(g) else np.nan,
                             **{o: int(vc.get(o, 0)) for o in ['correct_assignment', 'correct_novel', 'false_novel', 'wrong_allele', 'assigned_to_other']}))
    Sm = pd.DataFrame(rows); Sm.to_csv(os.path.join(OUT, 'loso_summary.tsv'), sep='\t', index=False); print(Sm.to_string())


# ------------------------------------------------------------------------------------------------ (c) independent assemblies, (d) runtime
def run_cli(inp, out, extra=()):
    """subtel-atlas query on a whole file via the CLI; returns (wall seconds, peak RSS bytes of the tool or its minimap2 child, log)."""
    import resource
    cmd = [sys.executable, '-m', 'subtel_atlas', 'query', inp, '--refs', REFS, '-t', str(THREADS), '--format', 'json', '-o', out] + list(extra)
    env = dict(os.environ, PYTHONPATH=os.path.join(REPO, 'atlas'))
    t = time.time(); p = subprocess.run(cmd, capture_output=True, text=True, env=env); wall = time.time() - t
    if p.returncode: raise SystemExit(p.stderr[-3000:])
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss          # bytes on macOS (kilobytes on Linux)
    rss = rss if sys.platform == 'darwin' else rss * 1024
    return wall, rss, '\n'.join(l for l in p.stderr.splitlines() if 'sysmon' not in l)


def _hgid(q): c, side = q.rsplit('#', 1); return f"HG002#{1 if c.endswith('PATERNAL') else 2}#{c}#{side}"
def _pedid(q):
    name, side = q.rsplit('#', 1); hap, chrom = name.split('|')[:2]; smp, h = hap.split('_'); return f'{smp}#{h}#{chrom}#{side}'


def indep(which='chm13,hg002,pedigree'):
    which = which.split(','); runs = []; ref = Q.Reference(REFS); allC = []
    jobs = {'chm13': (os.path.join(INPUT, 'ref', 'chm13v2.0.fa.gz'), lambda q: f'CHM13#0#{q}', 'CHM13 v2.0 (T2T), whole assembly'),
            'hg002': (os.path.join(REPO, 'work', 'qtool', 'hg002', 'hg002v1.1.fasta.gz'), _hgid, 'HG002 Q100 v1.1, whole diploid assembly'),
            'pedigree': (os.path.join(VW, 'pedigree_ends250k.fa'), _pedid, 'WashU pedigree v1.1, outer 250 kb of each chromosome end (validation-only set)')}
    if 'pedigree' in which and not os.path.exists(jobs['pedigree'][0]):
        with open(jobs['pedigree'][0], 'w') as f:
            for fn in ('ends250k.fa', 'ends250k_PAN028_h2.fa'):
                for n, sq in readfa(os.path.join(INPUT, 'subtel', 'pedigree', fn)): f.write(f'>{n}\n{sq}\n')
    for k in which:
        inp, idf, desc = jobs[k]; out = os.path.join(VW, f'indep_{k}.json')
        if not os.path.exists(out + '.time'):
            wall, rss, log = run_cli(inp, out); json.dump(dict(wall=wall, rss=rss, log=log), open(out + '.time', 'w'))
        tm = json.load(open(out + '.time')); d = json.load(open(out)); recs = d['ends']
        idmap = {r['query_id']: idf(r['query_id']) for r in recs}
        C = compare_rows(recs, idmap); C['dataset'] = k; allC.append(C)
        exp = E[E.set == {'chm13': 'CHM13', 'hg002': 'HG002_Q100', 'pedigree': 'pedigree_aux'}[k]]
        found = set(idmap.values()) & set(exp.end_id)
        runs.append(dict(dataset=k, description=desc, n_input_records=d['n_input_records'], n_telomeric_ends_found=len(recs),
                         n_catalogue_ends=len(exp), n_matched_to_catalogue=len(found), n_unmatched_tool_ends=len(set(idmap.values()) - set(exp.end_id)),
                         wall_seconds=round(tm['wall'], 1), peak_rss_gb=round(tm['rss'] / 1e9, 2), threads=THREADS, aligner=d['aligner']))
        # leave-one-sample-out for the T2T genomes (they define alleles in the catalogue)
        if k in ('chm13', 'hg002'):
            smp = 'CHM13' if k == 'chm13' else 'HG002'
            for r in recs: r['end_id_cat'] = idmap[r['query_id']]
            rs = [r for r in recs if r['end_id_cat'] in EI.index]
            Rl = loso_reference(ref, smp, {}, {}); L = pd.DataFrame(loso_eval(Rl, rs, smp, fetch_windows([r['end_id_cat'] for r in rs])))
            L['sa_outcome'] = [loso_outcome(a, b) for a, b in zip(L.sa_expected, L.sa_loso)]
            L['hg_outcome'] = [loso_outcome(a, b) for a, b in zip(L.hg_expected, L.hg_loso)]
            L.to_csv(os.path.join(OUT, f'independent_{k}_loso.tsv'), sep='\t', index=False)
            runs[-1].update(loso_sa_accuracy=round(float(L.sa_outcome.isin(['correct_assignment', 'correct_novel']).mean()), 4),
                            loso_hg_accuracy=round(float(L.hg_outcome.isin(['correct_assignment', 'correct_novel']).mean()), 4))
        if k == 'pedigree': pedigree_tx(recs, idmap)
    C = pd.concat(allC); C.to_csv(os.path.join(OUT, 'independent_ends.tsv'), sep='\t', index=False)
    Sm = pd.concat([summarise_cmp(g, f'independent_{k}') for k, g in C.groupby('dataset')])
    Sm = Sm[Sm.stratum.isin(['all', 'boundary_identical'])]; Sm.to_csv(os.path.join(OUT, 'independent_summary.tsv'), sep='\t', index=False)
    R_ = pd.DataFrame(runs); fn = os.path.join(OUT, 'runtime.tsv')
    if os.path.exists(fn):
        old = rd(fn); R_ = pd.concat([old[~old.dataset.isin(R_.dataset)], R_])
    R_.to_csv(fn, sep='\t', index=False); print(R_.to_string()); print(Sm.pivot_table(index='metric', columns=['analysis', 'stratum'], values='agreement').to_string())


def hg002_vs_r2(recs, idmap):
    """HG002 Q100 (tool, whole assembly) vs the HG002 Release 2 assemblies in the catalogue: per arm, the pair of calls."""
    q = collections.defaultdict(list)
    for r in recs: q[r['catalogue_arm']].append((nz(r['structural_allele_id']) or 'novel', nz(r['haplogroup_id']) or 'none', r['S20_state']))
    r2 = E[(E.set == 'HPRC_R2') & (E['sample'] == 'HG002')]; rows = []
    if not len(r2):
        print('HG002 has no Release 2 ends in this catalogue release: Q100-vs-Release-2 comparison not possible'); return
    for arm in sorted(set(q) | set(r2.catalogue_arm)):
        g = r2[r2.catalogue_arm == arm]
        a = sorted(x[0] for x in q.get(arm, [])); b = sorted(nz(x) or 'novel' for x in g.structural_allele_id)
        c = sorted(x[1] for x in q.get(arm, [])); d = sorted(nz(x) or 'none' for x in g.haplogroup_id)
        rows.append(dict(arm=arm, n_q100=len(a), n_r2=len(b), sa_q100=';'.join(a), sa_r2=';'.join(b), sa_concordant=a == b,
                         hg_q100=';'.join(c), hg_r2=';'.join(d), hg_concordant=c == d,
                         s20_q100=';'.join(sorted(x[2] for x in q.get(arm, []))), s20_r2=';'.join(sorted(map(str, g.S20_state))),
                         s20_concordant=sorted(x[2] for x in q.get(arm, [])) == sorted(map(str, g.S20_state))))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'independent_hg002_q100_vs_r2.tsv'), sep='\t', index=False)


def pedigree_tx(recs, idmap):
    TR = rd(os.path.join(REPO, 'work', 'catalogue', 'pedigree_transmission_alleles.tsv'))
    t = {idmap[r['query_id']]: r for r in recs}; rows = []
    for x in TR.itertuples():
        c, p_ = t.get(x.child_id), t.get(x.parent_id)
        if c is None or p_ is None: continue
        rows.append(dict(child_end=x.child_id, parent_end=x.parent_id, intact=x.intact, sa_child=nz(c['structural_allele_id']), sa_parent=nz(p_['structural_allele_id']),
                         sa_same=nz(c['structural_allele_id']) == nz(p_['structural_allele_id']) and nz(c['structural_allele_id']) != '',
                         hg_child=nz(c['haplogroup_id']), hg_parent=nz(p_['haplogroup_id']),
                         hg_same=nz(c['haplogroup_id']) == nz(p_['haplogroup_id']) and nz(c['haplogroup_id']) != '',
                         s20_same=c['S20_state'] == p_['S20_state'], core_structure_same=c['core_structure'] == p_['core_structure'],
                         sa_same_catalogue=nz(x.child_sa) == nz(x.parent_sa) and nz(x.child_sa) != '', hg_same_catalogue=nz(x.child_hg) == nz(x.parent_hg) and nz(x.child_hg) != ''))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'independent_pedigree_transmissions.tsv'), sep='\t', index=False)


# ------------------------------------------------------------------------------------------------ probe rule
def probes(n='300'):
    import mappy
    n = int(n); A = rd(os.path.join(REPO, 'work', 'catalogue', 'arms.tsv')); PR = rd(os.path.join(WIN, 'probe_plan.tsv'))
    W = rd(os.path.join(WIN, 'ends_windows.tsv'), usecols=['end_id', 'side'])
    A = A[A.probe_n > 0].merge(W, on='end_id')
    pick = pd.concat([A[A.arm_source == 'probe'].sample(n, random_state=SEED), A[A.arm_source == 'contig_name'].sample(n, random_state=SEED)])
    pr = PR[PR.end_id.isin(pick.end_id)]; seq = {}
    for asm, g in pr.groupby('asm'):
        want = set(g.region); have = {nm: sq for nm, sq in readfa(os.path.join(WIN, 'probes', f'{asm}.fa')) if nm in want}
        for r in g.itertuples():
            if r.region in have: seq.setdefault(r.end_id, []).append((int(r.probe_offset), have[r.region].upper()))
    mmi = os.path.join(REPO, 'work', 'catalogue', 'chm13v2.0.asm5.mmi')
    t = time.time(); gen = mappy.Aligner(mmi, preset='asm5'); tl = time.time() - t
    ends = [dict(end_id=r.end_id, side=r.side, probes=sorted(seq.get(r.end_id, []))) for r in pick.itertuples()]
    t = time.time(); mm = Q.probe_calls_minimap2(ends, mmi, THREADS, Q.B.find_minimap2()); tm = time.time() - t
    rows, clen = [], {}
    t = time.time()
    for r, e in zip(pick.itertuples(), ends):
        parm, pconf, _ = Q._probe_call(e, gen, clen)
        mparm, mpconf, _ = mm.get(r.end_id, ('', '', []))
        rows.append(dict(probe_arm_minimap2=mparm, probe_conf_minimap2=mpconf,
                         identical_minimap2=(mparm == nz(r.probe_arm)) and (mpconf == nz(r.probe_confidence)), end_id=r.end_id, arm_source=r.arm_source, contig_arm=nz(r.contig_arm), probe_arm_catalogue=nz(r.probe_arm),
                         probe_conf_catalogue=nz(r.probe_confidence), probe_arm_tool=parm, probe_conf_tool=pconf,
                         identical=(parm == nz(r.probe_arm)) and (pconf == nz(r.probe_confidence)),
                         group_correct_vs_name=group_of(parm) == group_of(nz(r.contig_arm)) if nz(r.contig_arm) and parm else None))
    P = pd.DataFrame(rows); P.to_csv(os.path.join(OUT, 'probe_rule.tsv'), sep='\t', index=False)
    print('index load s', round(tl, 1), 'mappy map s', round(time.time() - t, 1), 'minimap2 s', round(tm, 1))
    print(P.groupby('arm_source')[['identical', 'identical_minimap2']].agg(['sum', 'size']))
    print(P.dropna(subset=['group_correct_vs_name']).group_correct_vs_name.agg(['sum', 'size']))


# ------------------------------------------------------------------------------------------------ experimental read mode
def reads(per_end='4', datasets='revio_30hr:HiFi_Revio,ont_ul_dorado:ONT_UL'):
    import pysam
    per_end = int(per_end); ref = Q.Reference(REFS); out = []
    for ds in datasets.split(','):
        fnb, label = ds.split(':'); bam = os.path.join(INPUT, 'hg002', 'bams', f'{fnb}.bam'); cnt = collections.Counter(); ends, truth = [], {}
        with pysam.AlignmentFile(bam, 'rb') as f:
            L = dict(zip(f.references, f.lengths))
            for r in f.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence or r.query_length < 10000: continue
                c = r.reference_name; side = 'L' if r.reference_start < 50000 else ('R' if r.reference_end > L[c] - 50000 else None)
                if side is None: continue
                tid = f"HG002#{1 if c.endswith('PATERNAL') else 2}#{c}#{side}"
                if tid not in EI.index or cnt[tid] >= per_end: continue
                sq = Q.rc(r.query_sequence) if r.is_reverse else r.query_sequence
                e = Q.read_end(r.query_name, sq)
                if e is None: continue
                cnt[tid] += 1; e['end_id'] = f'{label}:{r.query_name}'; truth[e['end_id']] = tid; ends.append(e)
        for variant in ('arm_inferred', 'arm_given'):
            xs = []
            for e in ends:
                x = dict(e)
                if variant == 'arm_given': x['arm_hint'], x['arm_hint_source'] = EI.at[truth[e['end_id']], 'home_arm'], 'truth'
                xs.append(x)
            t = time.time(); recs = Q.analyse(xs, ref, 'minimap2', THREADS); sec = time.time() - t
            C = compare_rows(recs, truth); C['dataset'] = label; C['variant'] = variant; C['read_len'] = [len(x['seq']) for x in xs][:len(C)]
            C['seconds_per_read'] = sec / max(len(xs), 1); out.append(C)
            print(label, variant, len(xs), 'reads', len(cnt), 'ends', round(sec, 1), 's', flush=True)
    C = pd.concat(out); C.to_csv(os.path.join(OUT, 'reads_hg002.tsv'), sep='\t', index=False)
    Sm = C.groupby(['dataset', 'variant']).agg(n_reads=('end_id', 'size'), n_ends=('end_id', 'nunique'), arm=('arm_match', 'mean'),
                                                 S20=('S20_match', 'mean'), core_structure=('core_structure_match', 'mean'), sa=('sa_match', 'mean'),
                                                 hg=('hg_match', 'mean'), median_read_bp=('read_len', 'median')).round(4).reset_index()
    Sm.to_csv(os.path.join(OUT, 'reads_hg002_summary.tsv'), sep='\t', index=False); print(Sm.to_string())


# ------------------------------------------------------------------------------------------------ mappy backend
def mappy_backend(n='400'):
    n = int(n); recs, _ = load_chunks(); byid = {r['query_id']: r for r in recs}
    D = rd(os.path.join(OUT, 'query_set.tsv')); ids = list(D[D.end_id.isin(byid)].sample(n, random_state=SEED).end_id)
    seqs = fetch_windows(ids); ref = Q.Reference(REFS); ends = []
    for e in ids:
        r = EI.loc[e]; x = Q.window_end(e, f'arm={r.home_arm}' if r.home_arm != 'unassigned' else '', seqs[e], name_arms=False); ends.append(x)
    t = time.time(); mr = Q.analyse(ends, ref, 'mappy', THREADS); sec = time.time() - t
    rows = []
    for m in mr:
        a = byid[m['query_id']]
        rows.append(dict(end_id=m['query_id'], **{f'{k}_agree': nz(a[k]) == nz(m[k]) for k in
                    ['S5_state', 'S10_state', 'S20_state', 'S20_best_end', 'profile', 'core_structure', 'core_switches_sub', 'structural_allele_id', 'haplogroup_id', 'mosaic_arm']}))
    M_ = pd.DataFrame(rows); M_.to_csv(os.path.join(OUT, 'backend_mappy.tsv'), sep='\t', index=False)
    sm = M_.drop(columns='end_id').mean().round(4); sm['seconds_per_end_mappy'] = round(sec / n, 3)
    sm.rename('value').rename_axis('metric').reset_index().to_csv(os.path.join(OUT, 'backend_mappy_summary.tsv'), sep='\t', index=False); print(sm.to_string())


# ------------------------------------------------------------------------------------------------ boundary-shift robustness
def shift(n='600'):
    n = int(n); recs, _ = load_chunks(); byid = {r['query_id']: r for r in recs}
    D = rd(os.path.join(OUT, 'query_set.tsv')); D = D[D.end_id.isin(byid)]
    D['nat'] = [NAT_SA.get(nz(byid[e]['structural_allele_id']), 'unassigned') for e in D.end_id]
    D['full'] = [bool(byid[e]['full_length']) for e in D.end_id]
    D = D[D.full]; nn = D[D.nat != 'native']; na = D[D.nat == 'native']
    pick = pd.concat([nn.sample(min(n // 2, len(nn)), random_state=SEED), na.sample(n - min(n // 2, len(nn)), random_state=SEED)])
    seqs = fetch_windows(list(pick.end_id)); ref = Q.Reference(REFS)
    ins = 'CCCTAA' * 200
    variants = [('boundary', -500), ('boundary', 500), ('telomere', -500), ('telomere', -250), ('telomere', 250), ('telomere', 500)]
    rows = []; t0 = time.time()
    for kind, d in variants:
        ends = []
        for e in pick.end_id:
            r0 = byid[e]; J, T = int(r0['junk_len']), int(r0['telo_len']); s = seqs[e]; arm = f"arm={EI.at[e, 'home_arm']}"
            if kind == 'boundary':
                x = Q.window_end(e, f'{arm} junk_len={J} telo_len={T + d}', s, name_arms=False)
            else:
                if d > 0: s2 = s[:J] + ins[:d] + s[J:]
                else:
                    if T + d < 1000: continue
                    s2 = s[:J] + s[J - d:] + 'N' * (-d)
                x = Q.window_end(e, arm, s2[:len(s)], name_arms=False)
            ends.append(x)
        out = Q.analyse(ends, ref, 'minimap2', THREADS)
        for m in out:
            a = byid[m['query_id']]
            pred = a.get('allele_minus500' if d < 0 else 'allele_plus500') if kind == 'boundary' else None
            rows.append(dict(end_id=m['query_id'], variant=f'{kind}{d:+d}', nativity=NAT_SA.get(nz(a['structural_allele_id']), 'unassigned'),
                             s0_ref=a['s0'], s0_variant=m['s0'],
                             S20_stable=(m['S20_state'], nz(m['S20_family'])) == (a['S20_state'], nz(a['S20_family'])),
                             core_structure_stable=nz(m['core_structure']) == nz(a['core_structure']),
                             allele_stable=nz(m['structural_allele_id']) == nz(a['structural_allele_id']),
                             haplogroup_stable=nz(m['haplogroup_id']) == nz(a['haplogroup_id']),
                             allele_ref=nz(a['structural_allele_id']), allele_variant=nz(m['structural_allele_id']),
                             predicted_allele=nz(pred) if pred is not None else '',
                             realigned_equals_predicted=(nz(m['structural_allele_id']) == nz(pred)) if pred is not None else None,
                             tool_flag_boundary_stable_allele=nz(a.get('boundary_stable_allele'))))
        print(kind, d, len(out), round(time.time() - t0, 1), 's', flush=True)
    R_ = pd.DataFrame(rows); R_.to_csv(os.path.join(OUT, 'boundary_shift.tsv'), sep='\t', index=False)
    sm = []
    for (v, nat), g in list(R_.groupby(['variant', 'nativity'])) + [((v, 'all'), g) for v, g in R_.groupby('variant')]:
        sm.append(dict(variant=v, nativity=nat, n=len(g), **{k: round(float(g[k].mean()), 4) for k in
                       ['S20_stable', 'core_structure_stable', 'allele_stable', 'haplogroup_stable']},
                       realigned_equals_predicted=round(float(g.realigned_equals_predicted.dropna().astype(bool).mean()), 4) if g.realigned_equals_predicted.notna().any() else np.nan))
    S_ = pd.DataFrame(sm).sort_values(['variant', 'nativity']); S_.to_csv(os.path.join(OUT, 'boundary_shift_summary.tsv'), sep='\t', index=False)
    print(S_.to_string())


if __name__ == '__main__':
    st = sys.argv[1]
    globals()[st](*sys.argv[2:])
