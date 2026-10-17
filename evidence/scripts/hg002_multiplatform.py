#!/usr/bin/env python3
"""HG002 Q100 non-native ends: spanning reads from every range-readable public long-read alignment of the HG002 Q100 v1.1 project,
tested with the frozen read-evidence criteria v1 (evidence/CRITERIA.md; decision logic unchanged, atlas/subtel_atlas/readevidence.py).

Read capture (per dataset and contig end; primary alignments, -F 0x900):
  region = outer max(5 kb, s0 + min(switch, 100 kb) + 30 kb) bp of the contig end, where s0 = telomere + junk length of the end window
  and the switch is the largest catalogue switch of any claim at the same arm. Reads are then split into
    'end'      : primary alignment overlaps the outer 5 kb (the calibration read population; all reads >= 5 kb are evaluated);
    'proximal' : primary alignment lies proximal of the outer 5 kb (telomere-reaching reads whose tip is soft-clipped or
                 supplementary); evaluated only if they carry a terminal telomere by the module's own caller (reads without one can
                 be neither conclusive nor contradicting).
Roles (as in run_calibration.py): carrier (claimed contig end = haplotype-assignable set) and HG002 other-haplotype reads at the same
arm (strictly native = negative control; core-native / other / same structure reported only).

Claim-level verdict (CRITERIA.md section 4) on the pooled carrier reads of all primary (non-derived) datasets after removing duplicate
molecules (same read UUID / PacBio movie+ZMW in more than one dataset); ONT HERRO-corrected reads are derived from ONT ultra-long reads
and are reported separately, never counted. Per-dataset counts and verdicts are also reported.

Stages: fetch | eval | report | all. Outputs: evidence/hg002/{claims.tsv,claims_by_dataset.tsv,reads.tsv,controls_by_dataset.tsv,
fetch_manifest.tsv}; work/evidence/hg002/ (slices, genome SAM, evidence BAMs; git-ignored).
"""
import argparse, os, sys, csv, json, gzip, time, hashlib, subprocess, urllib.request, re, collections
import multiprocessing as mp_
from concurrent.futures import ThreadPoolExecutor

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subtel_atlas import readevidence as RE  # noqa: E402

Q100 = 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/'
# name: (relative path, mappy preset, platform label, derived)
DATASETS = {
    'ont_ul_dorado':   ('ont_r10_ul_dorado/hg002v1.1_ont_r10_ul_dorado.pri.bam', 'map-ont', 'ONT_R10_UL', False),
    'ont_lc24_ulk':    ('ont_lc24_ulk/hg002v1.1_ont_lc24_ulk.pri.bam', 'map-ont', 'ONT_R10_UL', False),
    'ont_duplex':      ('ont_r10_duplex/hg002v1.1_ont_r10_duplex.pri.bam', 'map-ont', 'ONT_R10_duplex', False),
    'ont_q28':         ('ont_epi2me_q28/hg002v1.1_ont_epi2me_q28.pri.bam', 'map-ont', 'ONT_R10', False),
    'ont_lc24_apk':    ('ont_lc24_apk/hg002v1.1_ont_lc24_apk.pri.bam', 'map-ont', 'ONT_R10', False),
    'ont_q28_lc24':    ('ont_q28_lc24/hg002v1.1_ont_q28_lc24.bam', 'map-ont', 'ONT_R10', False),
    'ont_r941_recall': ('ont_methylation_recall_oct2024/03_08_22_R941_HG002.pri.bam', 'map-ont', 'ONT_R9.4.1', False),
    'ont_herro':       ('ont_HERRO_corrected/hg002v1.1_ont_herro_corrected.pri.bam', 'map-ont', 'ONT_HERRO', True),
    'hifi_revio_3cells': ('hifi_revio_3cell/hg002v1.1_hifi_revio_3cells.bam', 'map-hifi', 'HiFi_Revio', False),
    'hifi_revio_pbmay24': ('hifi_revio_pbmay24/hg002v1.1_hifi_revio_pbmay24.pri.bam', 'map-hifi', 'HiFi_Revio', False),
    'hifi_revio_30hr': ('hifi_revio_30hr/HG002.hifi_revio_30hr.pri.bam', 'map-hifi', 'HiFi_Revio', False),
    'hifi_vega_24hr':  ('hifi_vega_24hr_dec24/HG002.hifi_vega_24hr_17.5kb.pri.bam', 'map-hifi', 'HiFi_Vega', False),
    'hifi_sprq_30hr':  ('sprq_revio_30hr_nov24/HG002.sprq_revio_30hr.pri.bam', 'map-hifi', 'HiFi_Revio_SPRQ', False),
    'hifi_sprq_24hr':  ('sprq_revio_24hr_nov24/HG002.sprq_revio_24hr.pri.bam', 'map-hifi', 'HiFi_Revio_SPRQ', False),
    'hifi_pb_sprq':    ('pb_hifi_sprq/hg002v1.1_pb_hifi_sprq.bam', 'map-hifi', 'HiFi_Revio_SPRQ', False),
    'hifi_dc':         ('hifi_dcv1.1/hg002v1.1_hifi_dcv1.1.pri.bam', 'map-hifi', 'HiFi_DeepConsensus', False),
}
END_FLANK = 5000
EXTRA = 30000
CAP = 100000
G = {}


def log(*a): print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


def sha(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def remote_size(url):
    req = urllib.request.Request(url, method='HEAD')
    with urllib.request.urlopen(req, timeout=60) as r: return int(r.headers.get('Content-Length', -1))


def index_for(url, cache, local_dirs=()):
    """Local copy of the public .bai/.crai (re-used from local_dirs when the byte size matches the remote; else downloaded)."""
    ext = '.crai' if url.endswith('.cram') else '.bai'; base = os.path.basename(url) + ext
    rs = remote_size(url + ext)
    for d in local_dirs:
        c = os.path.join(d, base)
        if os.path.exists(c) and os.path.getsize(c) == rs: return c
    os.makedirs(cache, exist_ok=True); fn = os.path.join(cache, base)
    if not (os.path.exists(fn) and os.path.getsize(fn) == rs):
        tmp = fn + '.part'; urllib.request.urlretrieve(url + ext, tmp); os.replace(tmp, fn)
    return fn


def fetch(url, idx, region, out, ref=None, tries=3):
    if os.path.exists(out) and os.path.exists(out + '.bai'): return out
    os.makedirs(os.path.dirname(out), exist_ok=True); tmp = out + '.part.bam'
    cmd = ['samtools', 'view', '-b', '-F', '0x900', '-o', tmp] + (['-T', ref] if ref else []) + ['-X', url, idx, region]
    for i in range(tries):
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(out))
        if r.returncode == 0: break
        if i == tries - 1: raise RuntimeError(f'{region}: {r.stderr[-300:]}')
        time.sleep(5)
    subprocess.run(['samtools', 'index', tmp], check=True)
    os.replace(tmp, out); os.replace(tmp + '.bai', out + '.bai')
    return out


def count(fn):
    return int(subprocess.run(['samtools', 'view', '-c', fn], capture_output=True, text=True, check=True).stdout)


def molecule_keys(name):
    """Molecule identity for de-duplication across datasets: ONT read UUIDs (duplex 'a;b' -> both), PacBio movie/ZMW."""
    if name.startswith('m') and '/' in name: return {'/'.join(name.split('/')[:2])}
    return set(re.split(r'[;,]', name))


def run_of(name, rec):
    if rec is not None and rec.has_tag('RG'): return rec.get_tag('RG')
    if name.startswith('m') and '/' in name: return name.split('/')[0]      # PacBio movie
    return ''


def load_ends(cat):
    import pandas as pd
    E = pd.read_csv(os.path.join(cat, 'ends.tsv'), sep='\t', low_memory=False)
    return E[E.set == 'HG002_Q100'].copy()


def arm_flanks(claims):
    """outer bp to fetch per recipient arm: s0 + min(switch, CAP) + EXTRA (largest over claims at the arm)"""
    fl = {}
    for c in claims:
        sw = c.get('switch_sub'); sw = CAP if sw is None else min(sw, CAP)
        fl[c['recipient_arm']] = max(fl.get(c['recipient_arm'], END_FLANK), int(c.get('s0') or 0) + sw + EXTRA)
    return fl


def slice_path(work, ds, contig, side): return os.path.join(work, 'slices', ds, f'{contig}_{side}.bam')


def stage_fetch(a, claims, H):
    fl = arm_flanks(claims)
    arms = set(fl)
    ends = H[H.home_arm.isin(arms)]
    jobs = []
    for ds in a.datasets:
        url = Q100 + DATASETS[ds][0]
        idx = index_for(url, os.path.join(a.work, 'index'), a.index_dir)
        for r in ends.itertuples():
            n = fl[r.home_arm]; L = int(r.contig_len)
            reg = f'{r.contig}:1-{n}' if r.side == 'L' else f'{r.contig}:{L - n + 1}-{L}'
            jobs.append((ds, url, idx, reg, slice_path(a.work, ds, r.contig, r.side), r.end_id))

    def go(j):
        ds, url, idx, reg, out, eid = j
        t0 = time.time(); fetch(url, idx, reg, out)
        return dict(dataset=ds, end_id=eid, url=url, index=os.path.basename(idx), region=reg, file=os.path.relpath(out, REPO),
                    n_records=count(out), sha256=sha(out), seconds=round(time.time() - t0, 1))
    res = []
    with ThreadPoolExecutor(a.jobs) as ex:
        for i, r in enumerate(ex.map(go, jobs)):
            res.append(r)
            if i % 25 == 0: log(f'fetched {i + 1}/{len(jobs)} {r["dataset"]} {r["region"]} {r["n_records"]} rec')
    man = os.path.join(a.out, 'fetch_manifest.tsv')
    with open(man, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]), delimiter='\t'); w.writeheader(); w.writerows(res)
    log(f'{len(res)} slices, {sum(r["n_records"] for r in res)} records -> {man}')


# ---------------------------------------------------------------------------------------------------------------------- eval
def build_jobs(claims, H, a, fl):
    from run_calibration import same_block
    byid = {c['claim_id']: c for c in claims}
    jobs = []
    for c in claims:
        other = H[(H.home_arm == c['recipient_arm']) & (H.end_id != c['end_id'])]
        for ds in a.datasets:
            rel, preset, plat, derived = DATASETS[ds]
            base = dict(claim=c['claim_id'], set=ds, preset=preset, platform=plat, sample='HG002', url=Q100 + rel, derived=derived)
            jobs.append(dict(base, role='carrier', bam=slice_path(a.work, ds, c['contig'], c['side']), side=c['side'],
                             contig=c['contig'], contig_len=int(c['contig_len']), source_end=c['end_id']))
            for o in other.itertuples():
                oc = byid.get(o.end_id)
                strict = isinstance(o.structure, str) and all(t in ('home', 'none') for t in o.structure.split('|'))
                role = ('hg002_otherhap_native' if strict else 'hg002_otherhap_corenative') if o.nativity_class == 'native' else \
                    ('hg002_otherhap_same' if same_block(c, oc) else 'hg002_otherhap_other')
                jobs.append(dict(base, role=role, bam=slice_path(a.work, ds, o.contig, o.side), side=o.side, contig=o.contig,
                                 contig_len=int(o.contig_len), source_end=o.end_id))
    for i, j in enumerate(jobs): j['job'] = i
    return jobs


def _ctx(cid, preset):
    k = (cid, preset)
    if k not in G['ctx']:
        G['ctx'][k] = RE.ClaimContext(G['claims'][cid], G['ends'], G['crit'], preset=preset, mask_cache=G['mask_cache'])
    return G['ctx'][k]


def _region_of(bam):
    """region string of a slice (from the manifest written at fetch time)"""
    return G['regions'].get(os.path.relpath(bam, REPO), '')


def work1(j):
    import pysam
    ctx = _ctx(j['claim'], j['preset']); th = ctx.th; out = []; Q = []; n_prox = n_prox_tel = n_short = 0
    motifs = RE.telomere_motifs(th, j['platform'])
    L = j['contig_len']
    for name, seq, rec in RE.iter_reads(j['bam'], min_len=0):
        if len(seq) < G['min_len']: n_short += 1; continue
        at_end = (rec.reference_start < END_FLANK) if j['side'] == 'L' else (rec.reference_end > L - END_FLANK)
        capture = 'end' if at_end else 'proximal'
        if not at_end:
            n_prox += 1
            t, o, tb = RE.orient(seq, j['side'], th, motifs)
            if not (tb and tb['junk_len'] <= th['telomere_max_offset'] and tb['telo_len'] >= th['telomere_min_bp']): continue
            n_prox_tel += 1
        meta = dict(dataset=j['set'], platform=j['platform'], role=j['role'], sample=j['sample'], job=j['job'], source_end=j['source_end'],
                    capture=capture, run=run_of(name, rec), source_url=j['url'], region=_region_of(j['bam']),
                    input_strand='-' if rec.is_reverse else '+', input_end=rec.reference_end, derived=j['derived'])
        meta.update(RE.assembly_fields(rec, ctx.c))
        r, u, t = RE.evaluate_read(ctx, name, seq, side=j['side'], meta=meta)
        r['_w'] = f"W|{j['set']}|{name}"; r['_a'] = f"A|{j['job']}|{name}"
        Q.append((r['_w'], t))
        if len(t) - r['anchor_from'] >= th['anchor_min_bp']: Q.append((r['_a'], t[r['anchor_from']:]))
        out.append((r, u))
    stats = dict(job=j['job'], n_short=n_short, n_proximal=n_prox, n_proximal_telomeric=n_prox_tel)
    return j, out, Q, stats


def init(claims, ends, crit, mask_cache, min_len, regions):
    G.update(claims=claims, ends=ends, crit=crit, mask_cache=mask_cache, min_len=min_len, regions=regions, ctx={})


def stage_eval(a, claims, H):
    import pickle
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    ends = RE.Chm13Ends(a.chm13_ends)
    fl = arm_flanks(claims)
    regions = {r['file']: r['region'] for r in csv.DictReader(open(os.path.join(a.out, 'fetch_manifest.tsv')), delimiter='\t')}
    jobs = build_jobs(claims, H, a, fl)
    cmap = {c['claim_id']: c for c in claims}
    mask_dir = os.path.join(REPO, 'work', 'readevidence', 'masks')
    for c in claims: RE.shared_mask_cached(ends, c['R'], th, mask_dir)
    log(f'{len(jobs)} jobs; criteria {crit["status"]} {RE.criteria_hash(crit)}')
    t0 = time.time()
    # longest jobs first
    jobs.sort(key=lambda j: -os.path.getsize(j['bam']))
    res = []
    with mp_.get_context('fork').Pool(a.threads, initializer=init, initargs=(cmap, ends, crit, mask_dir, a.min_len, regions)) as pool:
        for i, x in enumerate(pool.imap_unordered(work1, jobs, chunksize=1)):
            res.append(x)
            if i % 50 == 0: log(f'stage 1: {i + 1}/{len(jobs)} jobs, {time.time() - t0:.0f} s')
    log(f'stage 1 done: {time.time() - t0:.0f} s')
    Q = {}
    for _, out, q, _ in res:
        for k, s in q: Q.setdefault(k, s)
    sam = os.path.join(a.work, 'genome.sam')
    log(f'genome mapping of {len(Q)} queries ({sum(len(s) for s in Q.values()) / 1e6:.0f} Mb)')
    Gm = RE.map_genome(list(Q.items()), a.genome, a.threads, sam_out=sam)
    log(f'genome mapping done: {time.time() - t0:.0f} s')
    with open(os.path.join(a.work, 'stage.pkl'), 'wb') as f:
        pickle.dump(dict(res=[(j, out, st) for j, out, _, st in res], Gm=Gm), f)
    log('stage pickle written')


def tech(platform):
    """sequencing chemistry of a platform label: ONT R10.4.1, ONT R9.4.1 or PacBio HiFi"""
    return 'ONT_R9.4.1' if platform == 'ONT_R9.4.1' else ('PacBio_HiFi' if platform.startswith('HiFi') else 'ONT_R10.4.1')


def chm13_coord(ends_len, P, x):
    """CHM13 v2.0 chromosome coordinate (1-based) of tip-frame position x on end P"""
    if x is None or x == '' or (isinstance(x, float) and x != x): return None
    x = int(x); chrom = P[:-1]
    return x + 1 if P.endswith('p') else ends_len[chrom] - x


def reason_for(v, rows, c, asm):
    if v['verdict'] == 'confirmed': return ''
    sp = [r for r in rows if r.get('switch_pass')]
    if v['n_contradicts'] and v['n_conclusive_concordant'] >= 2:
        return f"{v['n_contradicts']} contradicting read(s) on the claimed haplotype contig block confirmation (phased rule)"
    if v['n_conclusive_concordant'] == 1: return 'only 1 conclusive switch-concordant read'
    if v['n_conclusive_discordant'] and not v['n_conclusive_concordant']:
        if c.get('switch_sub') is None: return 'conclusive reads, switch beyond the assembly window (concordance undefined in v1)'
        return 'conclusive reads are switch-discordant'
    if not sp:
        el = [r for r in rows if r.get('eligible')]
        if c.get('switch_sub') is None:
            return 'no read passes the switch test; switch beyond the assembly window (concordance undefined in v1)'
        if not el: return 'no telomeric read long enough to cross the switch with a 10 kb anchor (s0 + switch + 10 kb)'
        md = max([r.get('seg_D_units') or 0 for r in el], default=0)
        mp = max([r.get('seg_D_purity') or 0 for r in el], default=0)
        return (f'{len(el)} eligible read(s), none passes the switch test (max donor-segment units {md}, max donor purity {mp}; '
                'too few diagnostic units or mosaic donor block)')
    tp = [r for r in sp if r.get('terminal')]
    if not tp: return 'switch-passing reads lack a terminal telomere'
    ap = [r for r in tp if r.get('anchor_pass')]
    if not ap: return f'{len(tp)} telomeric switch-passing read(s), none with a MAPQ-60 anchor of >= 10 kb on the recipient chromosome'
    return 'switch-passing anchored reads fail integrity'


def stage_report(a, claims, H):
    import pickle, pandas as pd, pysam
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    ends = RE.Chm13Ends(a.chm13_ends)
    cmap = {c['claim_id']: c for c in claims}
    mask_dir = os.path.join(REPO, 'work', 'readevidence', 'masks')
    init(cmap, ends, crit, mask_dir, a.min_len, {})
    S = pickle.load(open(os.path.join(a.work, 'stage.pkl'), 'rb')); Gm = S['Gm']
    wf = pysam.FastaFile(a.windows); asm = {}
    for c in claims:
        asm[c['claim_id']] = RE.assembly_projection(_ctx(c['claim_id'], 'map-hifi'), wf.fetch(c['end_id']).upper()) \
            if c['end_id'] in wf.references else None
    fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.chm13_fai)}
    rows = []; stats = []
    for j, out, st in S['res']:
        ctx = _ctx(j['claim'], j['preset']); stats.append(dict(st, claim=j['claim'], role=j['role'], dataset=j['set']))
        for r, u in out:
            RE.finalize(ctx, r, Gm.get(r['_w']), Gm.get(r['_a'])); RE.concordance(r, asm[j['claim']], th)
            c = ctx.c
            r.update(claim_pattern=c['pattern'], D=c['D'], P=c['P'], R=c['R'])
            r['switch_ci_P_chm13_lo'] = chm13_coord(fai, c['P'], r.get('switch_ci_P_hi') if c['P'].endswith('q') else r.get('switch_ci_P_lo'))
            r['switch_ci_P_chm13_hi'] = chm13_coord(fai, c['P'], r.get('switch_ci_P_lo') if c['P'].endswith('q') else r.get('switch_ci_P_hi'))
            r['eligible'] = bool(r['terminal'] and c.get('switch_sub') is not None and
                                 r['read_len'] >= r['s0_read'] + c['switch_sub'] + th['anchor_min_bp'])
            r['molecule'] = ','.join(sorted(molecule_keys(r['read'])))
            rows.append(r)
    # ---- duplicate molecules across datasets (carrier reads): keep the first dataset in DATASETS order
    order = {d: i for i, d in enumerate(DATASETS)}
    seen = {}
    for r in sorted([r for r in rows if r['role'] == 'carrier'], key=lambda r: order[r['dataset']]):
        ks = molecule_keys(r['read']); key = (r['claim_id'],)
        dup = next((seen[(r['claim_id'], k)] for k in ks if (r['claim_id'], k) in seen), None)
        r['duplicate_of'] = dup or ''
        if not dup:
            for k in ks: seen[(r['claim_id'], k)] = f"{r['dataset']}:{r['read']}"
    # ---- claim tables
    E = H.set_index('end_id')
    out = []; byds = []
    for c in claims:
        cid = c['claim_id']; cr = [r for r in rows if r['claim_id'] == cid and r['role'] == 'carrier']
        prim = [r for r in cr if not r['derived'] and not r.get('duplicate_of')]
        v = RE.claim_verdict(prim, crit, same_hap_key='contig')
        vh = RE.claim_verdict([r for r in cr if r['derived']], crit, same_hap_key='contig')
        v3 = RE.claim_verdict([r for r in cr if r['dataset'] in ('ont_ul_dorado', 'hifi_revio_3cells', 'ont_duplex')
                               and r['capture'] == 'end'], crit, same_hap_key='contig')
        ctl = [r for r in rows if r['claim_id'] == cid and r['role'] == 'hg002_otherhap_native']
        e = E.loc[c['end_id']]; A = asm[cid] or {}
        conc_ds = collections.Counter(r['dataset'] for r in prim if r['verdict'] == 'conclusive' and r.get('switch_concordant'))
        plats = sorted({DATASETS[d][2] for d in conc_ds}); techs = sorted({tech(DATASETS[d][2]) for d in conc_ds})
        rec = dict(claim_id=cid, end_id=c['end_id'], contig=c['contig'], side=c['side'], haplotype=e['haplotype'],
                   transmitting_parent='HG003 (father)' if 'PATERNAL' in c['contig'] else 'HG004 (mother)',
                   allele_id=e['structural_allele_id'], nativity_class=e['nativity_class'], core_structure=e['core_structure'],
                   structure=e['structure'], pattern=c['pattern'], donor_D=c['D'], proximal_P=c['P'], recipient_R=c['R'],
                   assembly_switch_sub=c.get('switch_sub'), assembly_switch_sub_distal=c.get('switch_sub_distal'),
                   chm13_p_switch_tip=c.get('p_switch'), chm13_p_switch_coord=chm13_coord(fai, c['P'], c.get('p_switch')),
                   asm_projection_switch_pass=A.get('switch_pass'), asm_projection_sub_lo=A.get('switch_ci_sub_lo'),
                   asm_projection_sub_hi=A.get('switch_ci_sub_hi'), asm_projection_P_lo=A.get('switch_ci_P_lo'),
                   asm_projection_P_hi=A.get('switch_ci_P_hi'),
                   verdict=v['verdict'], reason=reason_for(v, prim, c, A),
                   n_reads=v['n_reads'], n_eligible=sum(r['eligible'] for r in prim),
                   n_conclusive=v['n_conclusive_concordant'] + v['n_conclusive_discordant'],
                   n_conclusive_concordant=v['n_conclusive_concordant'], n_conclusive_discordant=v['n_conclusive_discordant'],
                   n_consistent=v['n_consistent'], n_contradicts=v['n_contradicts'],
                   n_conclusive_proximal_capture=sum(r['verdict'] == 'conclusive' and r['capture'] == 'proximal' for r in prim),
                   n_duplicate_molecules_removed=sum(bool(r.get('duplicate_of')) for r in cr if not r['derived']),
                   datasets_conclusive_concordant=';'.join(f'{d}:{n}' for d, n in sorted(conc_ds.items(), key=lambda x: order[x[0]])),
                   n_datasets_conclusive=len(conc_ds), platforms_conclusive=','.join(plats), n_platforms_conclusive=len(plats),
                   technologies_conclusive=','.join(techs), n_technologies_conclusive=len(techs),
                   herro_verdict=vh['verdict'], herro_conclusive_concordant=vh['n_conclusive_concordant'], herro_contradicts=vh['n_contradicts'],
                   calibration_scope_verdict=v3['verdict'],
                   n_control_reads=len(ctl), n_control_conclusive=sum(r['verdict'] == 'conclusive' for r in ctl),
                   criteria_sha=RE.criteria_hash(crit), criteria_status=crit['status'])
        out.append(rec)
        for d in DATASETS:
            if d not in a.datasets: continue
            dr = [r for r in cr if r['dataset'] == d]
            vd = RE.claim_verdict(dr, crit, same_hap_key='contig')
            byds.append(dict(claim_id=cid, dataset=d, platform=DATASETS[d][2], derived=DATASETS[d][3], n_reads=len(dr),
                             n_eligible=sum(r['eligible'] for r in dr),
                             n_conclusive=vd['n_conclusive_concordant'] + vd['n_conclusive_discordant'],
                             n_conclusive_concordant=vd['n_conclusive_concordant'], n_consistent=vd['n_consistent'],
                             n_contradicts=vd['n_contradicts'], n_duplicates=sum(bool(r.get('duplicate_of')) for r in dr),
                             dataset_verdict=vd['verdict']))
    CL = pd.DataFrame(out); BD = pd.DataFrame(byds)
    # wide counts per dataset: eligible/conclusive/concordant/contradicting
    for d in [x for x in DATASETS if x in a.datasets]:
        s = BD[BD.dataset == d].set_index('claim_id')
        CL[f'counts_{d}'] = CL.claim_id.map(lambda k: '/'.join(str(int(s.loc[k, x])) for x in
                                                               ('n_eligible', 'n_conclusive', 'n_conclusive_concordant', 'n_contradicts')))
    os.makedirs(a.out, exist_ok=True)
    CL.to_csv(os.path.join(a.out, 'claims.tsv'), sep='\t', index=False)
    BD.to_csv(os.path.join(a.out, 'claims_by_dataset.tsv'), sep='\t', index=False)
    # controls by dataset
    ctl = pd.DataFrame([dict(claim_id=r['claim_id'], role=r['role'], dataset=r['dataset'], verdict=r['verdict'],
                             conc=bool(r.get('switch_concordant'))) for r in rows if r['role'] != 'carrier'])
    if len(ctl):
        cs = ctl.groupby(['role', 'dataset']).agg(n_reads=('verdict', 'size'), n_conclusive=('verdict', lambda x: int((x == 'conclusive').sum())),
                                                 n_contradicts=('verdict', lambda x: int((x == 'contradicts').sum()))).reset_index()
        cs.to_csv(os.path.join(a.out, 'controls_by_dataset.tsv'), sep='\t', index=False)
    # ---- reads.tsv: every conclusive or contradicting read (all roles)
    keep = [r for r in rows if r['verdict'] in ('conclusive', 'contradicts')]
    cols = ['claim_id', 'role', 'source_end', 'read', 'molecule', 'duplicate_of', 'dataset', 'platform', 'derived', 'source_url', 'region',
            'capture', 'run', 'read_len', 'orientation', 'input_strand', 'telo_offset', 'telo_len', 's0_read', 'n_units', 'segments',
            'seg_D_units', 'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'dR_units', 'dR_purity', 'log10_lr_vs_allP',
            'log10_lr_vs_allD', 'switch_read', 'switch_ci_read_lo', 'switch_ci_read_hi', 'switch_ci_sub_lo', 'switch_ci_sub_hi',
            'switch1_ci_sub_lo', 'switch1_ci_sub_hi', 'P', 'switch_ci_P_lo', 'switch_ci_P_hi', 'switch_ci_P_chm13_lo', 'switch_ci_P_chm13_hi',
            'switch_concordant', 'switch_concordant_catalogue', 'concordance_ref', 'donor_interval_units', 'donor_interval_P_purity',
            'anchor_chrom', 'anchor_tip_start', 'anchor_tip_end', 'anchor_mapq', 'anchor_len', 'anchor_unique_bp', 'anchor_tier',
            'internal_telomere_bp', 'adapter_hits', 'foldback', 'subtel_aligned_frac', 'input_contig', 'input_pos', 'input_end',
            'input_mapq', 'input_has_sa', 'input_qcov', 'assembly_concordant', 'eligible', 'verdict', 'D', 'R', 'claim_pattern']
    with open(os.path.join(a.out, 'reads.tsv'), 'w') as f: RE.write_tsv(keep, f, cols)
    with gzip.open(os.path.join(a.work, 'reads_all.tsv.gz'), 'wt') as f:
        RE.write_tsv(rows, f, cols + ['terminal', 'switch_pass', 'anchor_pass', 'integrity_pass', 'genome_alignments'])
    pd.DataFrame(stats).to_csv(os.path.join(a.work, 'job_stats.tsv'), sep='\t', index=False)
    write_bams(a, rows, keep)
    log(f'{len(rows)} read rows; {len(keep)} conclusive/contradicting; claims: ' +
        ', '.join(f'{k}={n}' for k, n in CL.verdict.value_counts().items()))


def clean_header(hd):
    """drop machine-specific directories from @PG command lines (published evidence BAMs carry no local paths)"""
    for pg in hd.get('PG', []):
        if 'CL' in pg: pg['CL'] = re.sub(r'(?:/[^ \t/]+)+/([^/ \t]+)', r'\1', pg['CL'])
    return hd


def write_bams(a, rows, keep):
    """Evidence BAMs of conclusive and contradicting reads: (a) CHM13 v2.0 alignments (whole read, tip frame, from the genome-wide
    mapping), (b) the original HG002 Q100 alignments (primary record of the slice). Read groups name the dataset and source URL;
    tags XC = claim, XV = verdict, XR = role, XK = switch concordance."""
    import pysam
    od = os.path.join(a.work, 'evidence_bam'); os.makedirs(od, exist_ok=True)
    rgs = {}
    for r in keep:
        rgs.setdefault(r['dataset'], dict(ID=r['dataset'], SM='HG002', PL='ONT' if r['platform'].startswith('ONT') else 'PACBIO',
                                          DS=f"{r['platform']}; source={r['source_url']}" + ('; derived (corrected) reads' if r['derived'] else '')))
    by_w = collections.defaultdict(list)
    for r in keep: by_w[r['_w']].append(r)
    # (a) CHM13
    sam = os.path.join(a.work, 'genome.sam')
    with pysam.AlignmentFile(sam) as fi:
        hd = fi.header.to_dict(); hd['RG'] = list(rgs.values())
        hd.setdefault('PG', []).append(dict(ID='readevidence', PN='subtel-atlas', CL='evidence/scripts/hg002_multiplatform.py'))
        hd = clean_header(hd)
        tmp = os.path.join(od, 'hg002_evidence.chm13v2.unsorted.bam')
        with pysam.AlignmentFile(tmp, 'wb', header=hd) as fo:
            for rec in fi:
                if rec.query_name not in by_w: continue
                rs = by_w[rec.query_name]; r0 = rs[0]
                rec.query_name = r0['read']; rec.set_tag('RG', r0['dataset'])
                rec.set_tag('XC', ','.join(sorted({r['claim_id'] for r in rs})))
                rec.set_tag('XV', ','.join(sorted({r['verdict'] for r in rs})))
                rec.set_tag('XR', ','.join(sorted({r['role'] for r in rs})))
                rec.set_tag('XK', ','.join(sorted({str(r.get('switch_concordant')) for r in rs})))
                fo.write(rec)
    out_a = os.path.join(od, 'hg002_evidence.chm13v2.bam')
    pysam.sort('-o', out_a, tmp); pysam.index(out_a); os.remove(tmp)
    # (b) HG002 Q100 contig alignments
    want = collections.defaultdict(dict)
    for r in keep: want[(r['dataset'], r['input_contig'])].setdefault(r['read'], []).append(r)
    hdr = None; recs = []
    for (ds, contig), names in want.items():
        for fn in [f for f in os.listdir(os.path.join(a.work, 'slices', ds)) if f.startswith(contig + '_') and f.endswith('.bam')]:
            with pysam.AlignmentFile(os.path.join(a.work, 'slices', ds, fn)) as fi:
                if hdr is None:
                    hdr = fi.header.to_dict(); hdr['RG'] = list(rgs.values()); hdr['PG'] = [dict(ID='readevidence', PN='subtel-atlas',
                                                                                                CL='evidence/scripts/hg002_multiplatform.py')]
                for rec in fi:
                    if rec.query_name in names and not rec.is_secondary and not rec.is_supplementary:
                        rs = names[rec.query_name]
                        recs.append((ds, rec.to_string(), sorted({r['claim_id'] for r in rs}), sorted({r['verdict'] for r in rs}),
                                     sorted({r['role'] for r in rs})))
    out_b = os.path.join(od, 'hg002_evidence.hg002q100.bam'); tmp = out_b + '.unsorted.bam'
    seen = set()
    with pysam.AlignmentFile(tmp, 'wb', header=clean_header(hdr)) as fo:
        H = fo.header
        for ds, s, cl, vv, rl in recs:
            rec = pysam.AlignedSegment.fromstring(s, H)
            k = (ds, rec.query_name)
            if k in seen: continue
            seen.add(k)
            for t in ('MM', 'ML', 'RG'):
                if rec.has_tag(t) and t != 'RG': rec.set_tag(t, None)
            rec.set_tag('RG', ds); rec.set_tag('XC', ','.join(cl)); rec.set_tag('XV', ','.join(vv)); rec.set_tag('XR', ','.join(rl))
            fo.write(rec)
    pysam.sort('-o', out_b, tmp); pysam.index(out_b); os.remove(tmp)
    log(f'evidence BAMs: {out_a}, {out_b}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['fetch', 'eval', 'report', 'all'])
    ap.add_argument('--datasets', default=','.join(DATASETS))
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--chm13-ends', help='chm13_ends500k.fa (subtel-atlas fetch-refs)')
    ap.add_argument('--chm13-fai', help='.fai of CHM13 v2.0 (chromosome lengths)')
    ap.add_argument('--genome', help='CHM13 v2.0 FASTA or minimap2 map-ont index')
    ap.add_argument('--windows', default=os.path.join(REPO, 'work', 'windows', 'ends.fa'))
    ap.add_argument('--index-dir', action='append', default=[], help='directories with local copies of the public .bai files')
    ap.add_argument('--work', default=os.path.join(REPO, 'work', 'evidence', 'hg002'))
    ap.add_argument('--out', default=os.path.join(REPO, 'evidence', 'hg002'))
    ap.add_argument('--min-len', type=int, default=5000)
    ap.add_argument('-t', '--threads', type=int, default=4)
    ap.add_argument('-j', '--jobs', type=int, default=4)
    a = ap.parse_args(); a.datasets = a.datasets.split(',')
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    cat = os.path.join(REPO, 'data', 'catalogue')
    H = load_ends(cat)
    ids = H[H.nativity_class != 'native'].end_id.tolist()
    claims = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=ids)
    log(f'{len(claims)} claims (all non-native HG002 Q100 ends); datasets: {",".join(a.datasets)}')
    if a.stage in ('fetch', 'all'): stage_fetch(a, claims, H)
    if a.stage in ('eval', 'all'): stage_eval(a, claims, H)
    if a.stage in ('report', 'all'): stage_report(a, claims, H)


if __name__ == '__main__':
    main()
