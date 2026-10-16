#!/usr/bin/env python3
"""Calibration runs of the read-evidence method (evidence/CRITERIA.md section 5).

Jobs = (claim, role, read set):
  carrier                : HG002 Q100 reads aligned to the claimed contig end (haplotype-assignable set of the claim)
  hg002_otherhap_native  : HG002 reads of the other haplotype at the same arm where that haplotype is strictly native (raw mosaic
                           structure only home/none; negative control)
  hg002_otherhap_corenative: other haplotype native only at the core level (non-home raw blocks < 10 kb; reported only)
  hg002_otherhap_other   : HG002 reads of the other haplotype where it carries a different non-native structure (reported only)
  hg002_otherhap_same    : the other haplotype carries the same donor block (a second carrier set; reported only)
  ont1kg_native_strict   : 1KG-ONT reads at the arm from individuals strictly native at the arm in both haplotypes (negative control)
  ont1kg_native_core     : 1KG-ONT reads from individuals core-native in both haplotypes, not strictly native (reported only)
HG002 claims are split into a training and a held-out half (sorted end_id, alternating). --phase selects which carrier claims run;
negative controls of the selected claims always run. Output: <out>/reads.tsv, <out>/claims.tsv, <out>/units.json.gz,
<out>/genome.sam (CHM13 alignments of the evaluated reads, git-ignored work copy).
"""
import argparse, os, sys, json, gzip, time, csv, multiprocessing as mp_
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import readevidence as RE  # noqa: E402

HG_SETS = {'ont_ul_dorado': ('map-ont', 'ONT_R10_UL'), 'ont_duplex': ('map-ont', 'ONT_R10_duplex'),
           'hifi_revio_3cells': ('map-hifi', 'HiFi_Revio')}
G = {}


def split(claims):
    ids = sorted(c['claim_id'] for c in claims)
    return {cid: ('train' if i % 2 == 0 else 'heldout') for i, cid in enumerate(ids)}


def same_block(a, b):
    return a is not None and b is not None and a['D'] == b['D'] and a['pattern'] == b['pattern'] and \
        a.get('switch_sub') is not None and b.get('switch_sub') is not None and abs(a['switch_sub'] - b['switch_sub']) <= 5000


def build_jobs(claims, hg_ends, k1g, work, phase, halves):
    import pandas as pd
    byid = {c['claim_id']: c for c in claims}
    jobs = []
    H = pd.read_csv(hg_ends, sep='\t')
    for c in claims:
        if phase != 'all' and halves[c['claim_id']] != phase: continue
        other = H[(H.home_arm == c['recipient_arm']) & (H.end_id != c['end_id'])]
        for ds, (preset, plat) in HG_SETS.items():
            f = os.path.join(work, 'hg002', ds, f"{c['contig']}_{c['side']}.bam")
            jobs.append(dict(claim=c['claim_id'], role='carrier', set=ds, bam=f, side=c['side'], preset=preset, platform=plat, sample='HG002'))
            for o in other.itertuples():
                oc = byid.get(o.end_id)
                strict = isinstance(o.structure, str) and all(t in ('home', 'none') for t in o.structure.split('|'))
                role = ('hg002_otherhap_native' if strict else 'hg002_otherhap_corenative') if o.nativity_class == 'native' else \
                    ('hg002_otherhap_same' if same_block(c, oc) else 'hg002_otherhap_other')
                f = os.path.join(work, 'hg002', ds, f"{o.contig}_{o.side}.bam")
                jobs.append(dict(claim=c['claim_id'], role=role, set=ds, bam=f, side=o.side, preset=preset, platform=plat, sample='HG002',
                                 source_end=o.end_id))
        for r in k1g:
            if r['home_arm'] != c['recipient_arm']: continue
            f = os.path.join(work, 'ont1kg', r['sample'], f"{r['chm13_end']}.bam")
            if not os.path.exists(f): continue
            jobs.append(dict(claim=c['claim_id'], role=r['role'], set='ont1kg_' + r['sample'], bam=f,
                             side='L' if r['chm13_end'].endswith('p') else 'R', preset='map-ont', platform='ONT_R9.4.1', sample=r['sample']))
    for i, j in enumerate(jobs): j['job'] = i
    return jobs


def _ctx(cid, preset):
    k = (cid, preset)
    if k not in G['ctx']:
        G['ctx'][k] = RE.ClaimContext(G['claims'][cid], G['ends'], G['crit'], preset=preset, mask_cache=G['mask_cache'])
    return G['ctx'][k]


def work1(j):
    ctx = _ctx(j['claim'], j['preset']); out = []; Q = []
    for name, seq, rec in RE.iter_reads(j['bam'], min_len=G['min_len']):
        meta = dict(dataset=j['set'], platform=j['platform'], role=j['role'], sample=j['sample'], job=j['job'],
                    source_end=j.get('source_end', j['claim'] if j['role'] == 'carrier' else ''))
        if rec is not None and not j['role'].startswith('ont1kg'): meta.update(RE.assembly_fields(rec, ctx.c))
        r, u, t = RE.evaluate_read(ctx, name, seq, side=j['side'], meta=meta)
        r['_w'] = f"W|{j['set']}|{name}"; r['_a'] = f"A|{j['job']}|{name}"
        Q.append((r['_w'], t))
        if len(t) - r['anchor_from'] >= ctx.th['anchor_min_bp']: Q.append((r['_a'], t[r['anchor_from']:]))
        out.append((r, u))
    return j, out, Q


def init(claims, ends, crit, mask_cache, min_len):
    G.update(claims=claims, ends=ends, crit=crit, mask_cache=mask_cache, min_len=min_len, ctx={})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phase', choices=['train', 'heldout', 'all'], required=True)
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--chm13-ends', required=True, help='chm13_ends500k.fa (subtel-atlas fetch-refs)')
    ap.add_argument('--genome', required=True, help='CHM13 v2.0 FASTA or minimap2 map-ont index')
    ap.add_argument('--work', default=os.path.join(REPO, 'work', 'readevidence'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--windows', default=os.path.join(REPO, 'work', 'windows', 'ends.fa'), help='tip-frame end windows (faidx-indexed)')
    ap.add_argument('--min-len', type=int, default=5000)
    ap.add_argument('-t', '--threads', type=int, default=8)
    a = ap.parse_args()
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    ends = RE.Chm13Ends(a.chm13_ends)
    cat = os.path.join(REPO, 'data', 'catalogue')
    hg_ends = os.path.join(a.work, 'hg002_ends.tsv')
    import pandas as pd
    ids = pd.read_csv(hg_ends, sep='\t').end_id.tolist()
    claims = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=ids)
    halves = split(claims)
    k1g = []
    for fn, role in (('ont1kg_strict_native_controls.tsv', 'ont1kg_native_strict'), ('ont1kg_native_controls.tsv', 'ont1kg_native_core')):
        seen = {(r['sample'], r['home_arm']) for r in k1g}
        for r in csv.DictReader(open(os.path.join(a.work, fn)), delimiter='\t'):
            if (r['sample'], r['home_arm']) not in seen: r['role'] = role; k1g.append(r)
    jobs = build_jobs(claims, hg_ends, k1g, a.work, a.phase, halves)
    os.makedirs(a.out, exist_ok=True)
    print(f'{len(jobs)} jobs ({a.phase}); criteria {crit["status"]} {RE.criteria_hash(crit)}', file=sys.stderr)
    t0 = time.time(); cmap = {c['claim_id']: c for c in claims}
    for c in claims: RE.shared_mask_cached(ends, c['R'], th, os.path.join(a.work, 'masks'))   # warm cache before forking
    ctxm = mp_.get_context('fork')
    with ctxm.Pool(a.threads, initializer=init, initargs=(cmap, ends, crit, os.path.join(a.work, 'masks'), a.min_len)) as pool:
        res = pool.map(work1, jobs, chunksize=1)
    print(f'stage 1: {time.time() - t0:.0f} s', file=sys.stderr)
    Q = {}
    for _, out, q in res:
        for k, s in q: Q.setdefault(k, s)
    Gm = RE.map_genome(list(Q.items()), a.genome, a.threads, sam_out=os.path.join(a.out, 'genome.sam'))
    print(f'genome mapping of {len(Q)} queries: {time.time() - t0:.0f} s', file=sys.stderr)
    init(cmap, ends, crit, os.path.join(a.work, 'masks'), a.min_len)
    import pysam
    wf = pysam.FastaFile(a.windows); asm = {}
    for c in claims:
        ctx = _ctx(c['claim_id'], 'map-hifi')
        asm[c['claim_id']] = RE.assembly_projection(ctx, wf.fetch(c['end_id']).upper()) if c['end_id'] in wf.references else None
    rows = []; units = {}
    for j, out, _ in res:
        ctx = _ctx(j['claim'], j['preset'])
        for r, u in out:
            RE.finalize(ctx, r, Gm.get(r['_w']), Gm.get(r['_a'])); RE.concordance(r, asm[j['claim']], th)
            r['half'] = halves[j['claim']]; r['claim_pattern'] = ctx.c['pattern']; r['D'] = ctx.c['D']; r['P'] = ctx.c['P']
            r['assembly_switch_sub'] = ctx.c.get('switch_sub'); r['assembly_p_switch'] = ctx.c.get('p_switch')
            rows.append(r); units[f"{j['job']}|{r['read']}"] = [(x['read_pos'], x['allele'][0], x['kind'][0], x['d_pos'], x['p_pos']) for x in u]
    cols = ['half', 'role', 'sample', 'source_end', 'claim_pattern', 'D', 'P', 'assembly_switch_sub', 'assembly_p_switch', 'job'] + RE.READ_COLS
    with open(os.path.join(a.out, 'reads.tsv'), 'w') as f: RE.write_tsv(rows, f, cols)
    with gzip.open(os.path.join(a.out, 'units.json.gz'), 'wt') as f: json.dump(units, f)
    # claim level: carrier set per dataset group (haplotype-assignable: the claimed contig)
    out = []
    for c in claims:
        rr = [r for r in rows if r['claim_id'] == c['claim_id'] and r['role'] == 'carrier']
        if not rr: continue
        v = RE.claim_verdict(rr, crit, same_hap_key='contig')
        ctl = [r for r in rows if r['claim_id'] == c['claim_id'] and r['role'] in ('hg002_otherhap_native', 'ont1kg_native_strict')]
        v.update(claim_id=c['claim_id'], half=halves[c['claim_id']], pattern=c['pattern'], D=c['D'], P=c['P'], R=c['R'],
                 switch_sub=c.get('switch_sub'), p_switch=c.get('p_switch'), core_structure=c['core_structure'],
                 n_control_reads=len(ctl), n_control_conclusive=sum(r['verdict'] == 'conclusive' for r in ctl),
                 n_control_contradicts=sum(r['verdict'] == 'contradicts' for r in ctl),
                 criteria_sha=RE.criteria_hash(crit), criteria_status=crit['status'],
                 **{'asm_' + k: v for k, v in (asm[c['claim_id']] or {}).items()})
        out.append(v)
    pd.DataFrame(out).to_csv(os.path.join(a.out, 'claims.tsv'), sep='\t', index=False)
    print(f'done: {len(rows)} read rows, {len(out)} claims in {time.time() - t0:.0f} s', file=sys.stderr)


if __name__ == '__main__':
    main()
