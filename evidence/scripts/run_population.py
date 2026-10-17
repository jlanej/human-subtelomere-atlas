#!/usr/bin/env python3
"""Read-level evidence in independent 1000 Genomes ONT reads (frozen criteria v1; evidence/CRITERIA.md).

Claims: for each HPRC Release 2 individual also sequenced in the 1000 Genomes ONT Vienna panel (Schloissnig, Pani et al. 2025,
Nature, doi:10.1038/s41586-025-09290-7; R9.4.1 reads aligned to CHM13 v2.0), every arm with a non-native (foreign or
within_family) tier A/B structural allele. A claim is individual x arm x allele; an individual carrying the same allele on both
haplotypes has one claim (model taken from haplotype 1, the other end is reported). Controls: individuals of the panel that are
strictly native (raw structure home/none) in both haplotypes at a claimed arm; each control read set is tested against the model
of every allele claimed at that arm.

Reads: primary alignments in the outer 5 kb of the recipient CHM13 end, plus the CHM13 donor end at the start of the claimed block
(tip coordinates [d_lo - 2 kb, d_lo + 5 kb]) because some carrier reads have their primary alignment on the donor end; reads are
deduplicated by name per claim. The read verdict and the claim rule are those of subtel_atlas.readevidence (unchanged, unphased
claim rule).

Steps (each resumable):
  plan      -> <work>/plan.json
  fetch     -> <work>/slices/<sample>/<end>_<lo>_<hi>.bam and <work>/fetch_manifest.tsv (<= 4 concurrent range requests)
  evaluate  -> <work>/reads_full.tsv, <work>/units.json.gz, <work>/genome.sam
  report    -> evidence/population/claims.tsv, reads.tsv, controls.tsv, alleles.tsv
  bams      -> <work>/bam/<claim>.chm13.bam and <claim>.own.bam (conclusive reads; read groups name the source CRAM URL)
Inputs: ATLAS_INPUT_DIR (CHM13 v2.0 FASTA under ref/, 500 kb CHM13 ends under subtel/foreign/ref/), or explicit options.
"""
import argparse, os, sys, json, gzip, csv, time, collections, subprocess, hashlib, tempfile
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp_
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subtel_atlas import readevidence as RE  # noqa: E402
import fetch_reads as FR  # noqa: E402

CRAM = FR.ONT1KG
PLATFORM = 'ONT_R9.4.1'
FLANK = 5000
G = {}


def strict(s):
    return isinstance(s, str) and all(t in ('home', 'none') for t in s.split('|'))


def model_ok(m, th):
    """the claim's D and P windows (readevidence.claim_windows) are non-empty; a model whose raw donor blocks give d_lo > d_hi
    (donor span not collinear on D) has no defined D window and is not assessable under v1"""
    if not m: return False
    dw, pw = RE.claim_windows(m, th)
    return dw[1] > dw[0] and pw[1] > pw[0]


def tip_to_region(end, lo, hi, fai):
    c = end[:-1]; L = fai[c]
    return f'{c}:{lo + 1}-{hi}' if end.endswith('p') else f'{c}:{L - hi + 1}-{L - lo}'


# --------------------------------------------------------------------------------------------------------------------- plan
def plan(a):
    import pandas as pd
    cat = os.path.join(REPO, 'data', 'catalogue')
    E = pd.read_csv(os.path.join(cat, 'ends.tsv'), sep='\t', low_memory=False)
    A = pd.read_csv(os.path.join(cat, 'alleles_structural.tsv'), sep='\t')
    ov = [l.strip() for l in open(a.samples) if l.strip()]
    E = E[(E['sample'].isin(ov)) & (E['set'] == 'HPRC_R2')].copy()
    E['sa_tier'] = E.structural_allele_id.map(A.set_index('allele_id').tier)
    sel = E[E.nativity_class.isin(['foreign', 'within_family']) & E.sa_tier.isin(['A', 'B'])]
    models = {c['end_id']: c for c in RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'),
                                                              end_ids=sel.end_id.tolist())}
    claims = []
    for (s, arm, sa), g in sel.groupby(['sample', 'home_arm', 'structural_allele_id']):
        g = g.sort_values(['haplotype', 'end_id']); ids = [e for e in g.end_id if e in models]
        other = E[(E['sample'] == s) & (E.home_arm == arm) & (~E.end_id.isin(g.end_id))]
        if len(g) == 2: oh = 'same_allele'
        elif not len(other): oh = 'not_assembled'
        elif len(other) > 1: oh = 'multiple_ends'
        else:
            o = other.iloc[0]
            oh = ('native_strict' if strict(o.structure) else 'native_core') if o.nativity_class == 'native' else \
                ('nonnative_other_allele' if isinstance(o.nativity_class, str) else 'unclassified')
        exp = {'same_allele': 1.0, 'native_strict': 0.5, 'native_core': 0.5}.get(oh)
        claims.append(dict(claim=f'{s}|{arm}|{sa}', sample=s, arm=arm, allele=sa, tier=g.sa_tier.iloc[0],
                           nativity_class=g.nativity_class.iloc[0], end_ids=list(g.end_id), model_end=ids[0] if ids else None,
                           model=models.get(ids[0]) if ids else None, n_haplotypes=len(g), other_haplotype=oh,
                           other_end=other.end_id.iloc[0] if len(other) == 1 else None,
                           other_allele=other.structural_allele_id.iloc[0] if len(other) == 1 else None, expected_fraction=exp))
    if a.only:   # restrict to the listed claims (e.g. to re-derive one verdict): their arms' controls are kept
        keep = set(a.only.split(',')); claims = [c for c in claims if c['claim'] in keep]
        if not claims: raise SystemExit(f'no claim among {sorted(keep)}')
    arms = sorted({c['arm'] for c in claims}); controls = []
    for arm in arms:
        for s in ov:
            g = E[(E.home_arm == arm) & (E['sample'] == s)]
            if len(g) == 2 and all(strict(x) for x in g.structure): controls.append(dict(sample=s, arm=arm, end_ids=list(g.end_id)))
    if a.no_controls: controls = []
    # allele models used for controls: first claim (sorted) of each (arm, allele)
    th = RE.load_criteria(a.criteria)['thresholds']
    for c in claims: c['model_valid'] = model_ok(c['model'], th)
    amod = {}
    for c in sorted(claims, key=lambda c: c['claim']):
        if c['model_valid']: amod.setdefault(f"{c['arm']}|{c['allele']}", c['claim'])
    fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.chm13 + '.fai')}
    need = collections.defaultdict(list)
    for c in claims:
        m = c['model']
        if not m: continue
        need[(c['sample'], m['R'])].append((0, FLANK)); need[(c['sample'], m['D'])].append((max(0, m['d_lo'] - 2000), m['d_lo'] + FLANK))
    for k in controls: need[(k['sample'], 'chr' + k['arm'])].append((0, FLANK))
    slices = []
    for (s, end), iv in sorted(need.items()):
        cur = None
        for lo, hi in sorted(iv):
            if cur and lo <= cur[1] + 10000: cur[1] = max(cur[1], hi)
            else:
                if cur: slices.append(dict(sample=s, end=end, lo=cur[0], hi=cur[1]))
                cur = [lo, hi]
        slices.append(dict(sample=s, end=end, lo=cur[0], hi=cur[1]))
    for x in slices:
        x['region'] = tip_to_region(x['end'], x['lo'], x['hi'], fai); x['url'] = CRAM.format(s=x['sample'])
        x['file'] = os.path.join('slices', x['sample'], f"{x['end']}_{x['lo']}_{x['hi']}.bam")
    os.makedirs(a.work, exist_ok=True)
    json.dump(dict(criteria_sha=RE.criteria_hash(RE.load_criteria(a.criteria)), claims=claims, controls=controls,
                   allele_models=amod, slices=slices), open(os.path.join(a.work, 'plan.json'), 'w'), indent=1, default=str)
    print(f'{len(claims)} claims ({sum(c["model"] is None for c in claims)} without model), {len(controls)} control sample-arms, '
          f'{len(amod)} allele models, {len(slices)} slices', file=sys.stderr)


# --------------------------------------------------------------------------------------------------------------------- fetch
def fetch(a):
    P = json.load(open(os.path.join(a.work, 'plan.json')))
    idx_dirs = [os.path.join(REPO, 'work', 'readevidence', 'index')] + a.index_dir

    def go(x):
        idx = FR.index_for(x['url'], idx_dirs[0], idx_dirs)
        out = os.path.join(a.work, x['file']); t0 = time.time()
        FR.fetch(x['url'], idx, x['region'], out, a.chm13)
        return dict(sample=x['sample'], end=x['end'], url=x['url'], region=x['region'], file=x['file'], n_reads=FR.count(out),
                    sha256=FR.sha(out), seconds=round(time.time() - t0, 1))
    for s in sorted({x['sample'] for x in P['slices']}): FR.index_for(CRAM.format(s=s), idx_dirs[0], idx_dirs)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(go, P['slices']))
    with open(os.path.join(a.work, 'fetch_manifest.tsv'), 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]), delimiter='\t'); w.writeheader(); w.writerows(res)
    print(f'{len(res)} slices, {sum(r["n_reads"] for r in res)} records', file=sys.stderr)


# --------------------------------------------------------------------------------------------------------------------- evaluate
def build_jobs(P):
    sl = collections.defaultdict(list)
    for x in P['slices']: sl[(x['sample'], x['end'])].append(x)
    cm = {c['claim']: c for c in P['claims']}; jobs = []
    for c in P['claims']:
        m = c['model']
        if not c.get('model_valid'): continue
        dlo, dhi = max(0, m['d_lo'] - 2000), m['d_lo'] + FLANK      # the claim's own donor interval (tip frame on D)
        ss = sl[(c['sample'], m['R'])] + ([] if m['D'] == m['R'] else
                                          [x for x in sl[(c['sample'], m['D'])] if x['lo'] < dhi and dlo < x['hi']])
        jobs.append(dict(claim=c['claim'], model_claim=c['claim'], role='carrier', sample=c['sample'], slices=ss))
    for k in P['controls']:
        for am, mc in P['allele_models'].items():
            if am.split('|')[0] != k['arm']: continue
            jobs.append(dict(claim=f"{k['sample']}|{k['arm']}|{am.split('|')[1]}", model_claim=mc, role='control_native_strict',
                             sample=k['sample'], slices=sl[(k['sample'], 'chr' + k['arm'])]))
    jobs = [j for j in jobs if j['slices']]
    for i, j in enumerate(jobs): j['job'] = i
    return jobs, cm


def _ctx(mc, preset='map-ont'):
    k = (mc, preset)
    if k not in G['ctx']:
        G['ctx'][k] = RE.ClaimContext(G['cm'][mc]['model'], G['ends'], G['crit'], preset=preset, mask_cache=G['mask_cache'])
    return G['ctx'][k]


def _cache(j):
    import hashlib
    k = hashlib.sha1(json.dumps([j['claim'], j['role'], j['model_claim'], [x['file'] for x in j['slices']], G['min_len']]).encode()).hexdigest()[:16]
    return os.path.join(G['work'], 'stage1', f'{k}.pkl')


def work1(j):
    import pickle
    fn = _cache(j)
    if os.path.exists(fn):
        jj, out, Q = pickle.load(open(fn, 'rb')); return j, out, Q
    j, out, Q = _work1(j)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + '.tmp', 'wb') as f: pickle.dump((j, out, Q), f)
    os.replace(fn + '.tmp', fn)
    return j, out, Q


def _work1(j):
    ctx = _ctx(j['model_claim']); out = []; Q = []; seen = set()
    for x in j['slices']:
        side = 'L' if x['end'].endswith('p') else 'R'
        for name, seq, rec in RE.iter_reads(os.path.join(G['work'], x['file']), min_len=G['min_len']):
            if name in seen: continue
            seen.add(name)
            meta = dict(dataset='1KG_ONT_' + j['sample'], platform=PLATFORM, role=j['role'], sample=j['sample'], job=j['job'],
                        cram_url=x['url'], slice_region=x['region'], slice_end=x['end'],
                        read_group=rec.get_tag('RG') if rec is not None and rec.has_tag('RG') else '',
                        input_contig=rec.reference_name, input_pos=rec.reference_start, input_mapq=rec.mapping_quality,
                        input_has_sa=rec.has_tag('SA'), input_strand='-' if rec.is_reverse else '+')
            r, u, t = RE.evaluate_read(ctx, name, seq, side=side, meta=meta)
            r['claim_id'] = j['claim']; r['model_claim'] = j['model_claim']
            r['_w'] = f"W|{j['sample']}|{name}"; r['_a'] = f"A|{j['sample']}|{name}|{r['anchor_from']}"
            Q.append((r['_w'], t))
            if len(t) - r['anchor_from'] >= ctx.th['anchor_min_bp']: Q.append((r['_a'], t[r['anchor_from']:]))
            out.append((r, [(v['read_pos'], v['allele'][0], v['kind'][0], v['d_pos'], v['p_pos']) for v in u]))
    return j, out, Q


def init(cm, ends, crit, mask_cache, min_len, work):
    G.update(cm=cm, ends=ends, crit=crit, mask_cache=mask_cache, min_len=min_len, work=work, ctx={})


def evaluate(a):
    P = json.load(open(os.path.join(a.work, 'plan.json')))
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    assert RE.criteria_hash(crit) == P['criteria_sha'] == 'a89c5eaf023ed032', 'criteria changed'
    ends = RE.Chm13Ends(a.chm13_ends); jobs, cm = build_jobs(P)
    mask_cache = os.path.join(REPO, 'work', 'readevidence', 'masks')
    for R in sorted({c['model']['R'] for c in P['claims'] if c['model']}): RE.shared_mask_cached(ends, R, th, mask_cache)
    print(f'{len(jobs)} jobs; criteria {crit["status"]} {RE.criteria_hash(crit)}', file=sys.stderr); t0 = time.time()
    with mp_.get_context('fork').Pool(a.threads, initializer=init, initargs=(cm, ends, crit, mask_cache, a.min_len, a.work)) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(work1, sorted(jobs, key=lambda j: -len(j['slices'])), chunksize=1)):
            res.append(r)
            if i % 100 == 0: print(f'  stage 1: {i}/{len(jobs)} jobs, {time.time() - t0:.0f} s', file=sys.stderr, flush=True)
    print(f'stage 1: {time.time() - t0:.0f} s', file=sys.stderr)
    Q = {}
    for _, out, q in res:
        for k, s in q: Q.setdefault(k, s)
    Gm = RE.map_genome(list(Q.items()), a.genome, a.threads, sam_out=os.path.join(a.work, 'genome.sam'))
    print(f'genome mapping of {len(Q)} queries: {time.time() - t0:.0f} s', file=sys.stderr)
    init(cm, ends, crit, mask_cache, a.min_len, a.work)
    import pysam
    wf = pysam.FastaFile(a.windows); asm = {}
    for c in P['claims']:
        if not c.get('model_valid'): continue
        ctx = _ctx(c['claim'], 'map-hifi')
        for e in c['end_ids']:
            asm[(c['claim'], e)] = RE.assembly_projection(ctx, wf.fetch(e).upper()) if e in wf.references else None
    rows = []; units = {}
    for j, out, _ in res:
        ctx = _ctx(j['model_claim']); c = cm[j['model_claim']]
        for r, u in out:
            RE.finalize(ctx, r, Gm.get(r['_w']), Gm.get(r['_a']))
            # concordance: against the model end's assembly projection (claim rule); the second end of a homozygous claim and
            # the catalogue switch are reported alongside
            if len(c['end_ids']) > 1:
                r2 = dict(r); RE.concordance(r2, asm.get((c['claim'], c['end_ids'][1])), th)
                r['switch_concordant_hap2'] = r2.get('switch_concordant')
            RE.concordance(r, asm.get((c['claim'], c['model_end'])), th)
            r.update(claim_pattern=ctx.c['pattern'], D=ctx.c['D'], P=ctx.c['P'], R=ctx.c['R'],
                     assembly_switch_sub=ctx.c.get('switch_sub'), assembly_p_switch=ctx.c.get('p_switch'))
            rows.append(r); units[f"{j['job']}|{r['read']}"] = u
    cols = ['role', 'sample', 'model_claim', 'claim_pattern', 'D', 'P', 'R', 'assembly_switch_sub', 'assembly_p_switch', 'job',
            'cram_url', 'slice_region', 'slice_end', 'read_group', 'input_strand', 'switch_concordant_hap2'] + RE.READ_COLS
    with open(os.path.join(a.work, 'reads_full.tsv'), 'w') as f: RE.write_tsv(rows, f, cols)
    with gzip.open(os.path.join(a.work, 'units.json.gz'), 'wt') as f: json.dump(units, f)
    asmj = {f'{k[0]}||{k[1]}': v for k, v in asm.items()}
    json.dump(asmj, open(os.path.join(a.work, 'assembly_projection.json'), 'w'), default=str)
    print(f'done: {len(rows)} read rows in {time.time() - t0:.0f} s', file=sys.stderr)


# --------------------------------------------------------------------------------------------------------------------- main
def main():
    inp = os.environ.get('ATLAS_INPUT_DIR', '')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['plan', 'fetch', 'evaluate', 'report', 'bams', 'all'])
    ap.add_argument('--samples', default=os.path.join(REPO, 'work', 'ont1kg', 'overlap_samples.txt'))
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--chm13', default=os.path.join(inp, 'ref', 'chm13v2.0.fa.gz'))
    ap.add_argument('--chm13-ends', default=os.path.join(inp, 'subtel', 'foreign', 'ref', 'chm13_ends500k.fa'))
    ap.add_argument('--genome', default=os.path.join(REPO, 'work', 'readevidence', 'chm13v2.map-ont.mmi'))
    ap.add_argument('--windows', default=os.path.join(REPO, 'work', 'windows', 'ends.fa'))
    ap.add_argument('--work', default=os.path.join(REPO, 'work', 'evidence', 'population'))
    ap.add_argument('--out', default=os.path.join(REPO, 'evidence', 'population'))
    ap.add_argument('--index-dir', action='append', default=[])
    ap.add_argument('--only', help='comma-separated claim ids (<sample>|<arm>|<allele>) to plan; default all')
    ap.add_argument('--no-controls', action='store_true', help='plan without control sample-arms (with --only)')
    ap.add_argument('--min-len', type=int, default=5000)
    ap.add_argument('-j', '--jobs', type=int, default=4, help='concurrent range requests (fetch)')
    ap.add_argument('-t', '--threads', type=int, default=5)
    a = ap.parse_args()
    steps = ['plan', 'fetch', 'evaluate', 'report', 'bams'] if a.step == 'all' else [a.step]
    for s in steps:
        if s in ('report', 'bams'):
            import population_report as PR
            getattr(PR, s)(a)
        else:
            globals()[s](a)


if __name__ == '__main__':
    main()
