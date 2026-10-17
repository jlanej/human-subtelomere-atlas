#!/usr/bin/env python3
"""Trio transmission test: do the HG002 parents carry the HG002 non-native terminal blocks? Public, range-readable long-read
alignments of HG003 (father) and HG004 (mother), GIAB Ashkenazi trio (NIST RM8392), tested with the frozen criteria v1
(evidence/CRITERIA.md; module atlas/subtel_atlas/readevidence.py, decision logic unchanged).

Sources (SOURCES below): PacBio HiFi Revio 2023 and HiFi Sequel II (15/20 kb, chemistry 2) aligned by GIAB to CHM13 v2.0; UCSC
ultra-long ONT R9.4.1 (2019) aligned to GRCh38; ONT R10.4.1 sup per flow cell (giab_2023.05, ONT open data) aligned to GRCh38
(CRAM; decoded with the GRCh38 no-alt analysis set). The reference of the public alignment is used only to retrieve reads; every read
is evaluated assembly-free against CHM13 v2.0 (diagnostic units and the MAPQ-60 genome-wide anchor), as in the calibration.

Read capture per (source, claim): primary alignments in
  * the recipient arm R: outer 3 kb + min(max(switch, CHM13 projection), 100 kb) + 30 kb of the reference chromosome end;
  * the outer 5 kb of the donor end D and of the proximal comparator end P (reads whose tip block aligns to D or P).
Only reads >= 5 kb with a terminal telomere (module caller, platform motifs) are evaluated: reads without one can be neither conclusive
nor contradicting. Parents are unphased: claim rule of CRITERIA.md section 4 for unphased samples (>= 2 conclusive switch-concordant
reads confirm; contradicting reads, expected from the parent's other haplotype, are reported, not blocking). Switch concordance is
tested against the HG002 assembly projection of the claimed (transmitted) haplotype. Every claim is tested in both parents: the
transmitting parent (HG003 for PATERNAL, HG004 for MATERNAL claims) is the test; the other parent is reported only.

Stages: fetch | eval | report | all. Outputs: evidence/hg002/trio_claims.tsv, trio_reads.tsv, trio_fetch_manifest.tsv;
work/evidence/trio/ (slices, genome SAM, evidence BAM; git-ignored).
"""
import argparse, os, sys, csv, gzip, time, collections, pickle
import multiprocessing as mp_
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hg002_multiplatform as HM  # noqa: E402
from hg002_multiplatform import RE, REPO, log, fetch, count, sha, index_for, chm13_coord  # noqa: E402

GIAB = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/'
ONT = 'https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/'
F, M = 'HG003_NA24149_father/', 'HG004_NA24143_mother/'
# name: (sample, url, reference, preset, platform)
SOURCES = {
    'HG003_hifi_revio': ('HG003', GIAB + F + 'PacBio_HiFi-Revio_20231031/HG003_PacBio-HiFi-Revio_20231031_46x_CHM13v2.0.bam', 'CHM13', 'map-hifi', 'HiFi_Revio'),
    'HG004_hifi_revio': ('HG004', GIAB + M + 'PacBio_HiFi-Revio_20231031/HG004_PacBio-HiFi-Revio_20231031_36x_CHM13v2.0.bam', 'CHM13', 'map-hifi', 'HiFi_Revio'),
    'HG003_hifi_sequel2': ('HG003', GIAB + F + 'PacBio_CCS_15kb_20kb_chemistry2/CHM13v2.0/GIAB_5mC_CpG/HG003.chm13.deepvariant.haplotagged.bam', 'CHM13', 'map-hifi', 'HiFi_SequelII'),
    'HG004_hifi_sequel2': ('HG004', GIAB + M + 'PacBio_CCS_15kb_20kb_chemistry2/CHM13v2.0/GIAB_5mC_CpG/HG004.chm13.deepvariant.haplotagged.bam', 'CHM13', 'map-hifi', 'HiFi_SequelII'),
    'HG003_ont_ul_r941': ('HG003', GIAB + F + 'UCSC_Ultralong_OxfordNanopore_Promethion/HG003_GRCh38_ONT-UL_UCSC_20200508.bam', 'GRCh38', 'map-ont', 'ONT_R9.4.1'),
    'HG004_ont_ul_r941': ('HG004', GIAB + M + 'UCSC_Ultralong_OxfordNanopore_Promethion/HG004_GRCh38_ONT-UL_UCSC_20200508.bam', 'GRCh38', 'map-ont', 'ONT_R9.4.1'),
}
for s, fcs in (('HG003', ['PAO98397_289a47d8', 'PAO98397_4e702cd4', 'PAQ00476_291af615']),
               ('HG004', ['PAO98462_78a9a38d', 'PAO98462_947941d1', 'PAO98462_c3dbcd08', 'PAO98633_064894cd', 'PAO98633_162f2f43',
                          'PAO98633_9ddb25e9'])):
    for fc in fcs:
        SOURCES[f'{s}_ont_r10_{fc}'] = (s, ONT + f'{s.lower()}/sup/{fc}.pass.cram', 'GRCh38', 'map-ont', 'ONT_R10')
D_FLANK = 5000
G = {}


def chrom_lengths(url, ref_fa=None):
    import pysam
    kw = dict(reference_filename=ref_fa) if (ref_fa and url.endswith('.cram')) else {}
    with pysam.AlignmentFile(url, 'rc' if url.endswith('.cram') else 'rb', **kw) as f:
        return dict(zip(f.references, f.lengths))


def terminal_n(fasta, chrom, side, step=100000, cap=25000000):
    """length of the terminal N run of a reference chromosome (GRCh38 telomere gaps; 0 for CHM13)"""
    import pysam
    with pysam.FastaFile(fasta) as f:
        L = f.get_reference_length(chrom); n = 0
        while n < min(cap, L):
            s = (f.fetch(chrom, n, min(n + step, L)) if side == 'L' else f.fetch(chrom, max(L - n - step, 0), L - n)).upper()
            k = len(s) - len(s.lstrip('N')) if side == 'L' else len(s) - len(s.rstrip('N'))
            n += k
            if k < len(s): break
    return n


def arm_region(arm_end, lens, n, off=0):
    """outer n bp of the reference chromosome end, beyond a terminal N run of length off"""
    chrom = arm_end[:-1]; L = lens[chrom]
    return (f'{chrom}:{off + 1}-{off + n}', 'L') if arm_end.endswith('p') else (f'{chrom}:{L - off - n + 1}-{L - off}', 'R')


def plan(claims, a):
    """slices to fetch: per source, the union over claims of R-arm regions (largest flank) and D/P end flanks"""
    need = {}
    for c in claims:
        sw = max([x for x in (c.get('switch_sub'), c.get('p_switch')) if x is not None], default=HM.CAP)
        nR = 3000 + min(sw, HM.CAP) + HM.EXTRA
        need[c['R']] = max(need.get(c['R'], 0), nR)
        for e in (c['D'], c['P']):
            need[e] = max(need.get(e, 0), D_FLANK)
    return need


def slice_path(work, src, end, region): return os.path.join(work, 'slices', src, f"{end}_{region.split(':')[1]}.bam")


def stage_fetch(a, claims):
    need = plan(claims, a); jobs = []; lens = {}
    for src in a.sources:
        s, url, ref, preset, plat = SOURCES[src]
        idx = index_for(url, os.path.join(a.work, 'index'))
        lens[src] = chrom_lengths(url, a.grch38 if ref == 'GRCh38' else None)
        for end, n in sorted(need.items()):
            off = terminal_n(a.grch38, end[:-1], 'L' if end.endswith('p') else 'R') if ref == 'GRCh38' else 0
            reg, side = arm_region(end, lens[src], n, off)
            jobs.append((src, url, idx, reg, slice_path(a.work, src, end, reg), end, a.grch38 if url.endswith('.cram') else None))

    def go(j):
        src, url, idx, reg, out, end, ref = j
        t0 = time.time(); fetch(url, idx, reg, out, ref)
        return dict(source=src, sample=SOURCES[src][0], reference=SOURCES[src][2], end=end, url=url, region=reg,
                    file=os.path.relpath(out, REPO), n_records=count(out), sha256=sha(out), seconds=round(time.time() - t0, 1))
    res = []
    with ThreadPoolExecutor(a.jobs) as ex:
        for i, r in enumerate(ex.map(go, jobs)):
            res.append(r)
            if i % 20 == 0: log(f'fetched {i + 1}/{len(jobs)} {r["source"]} {r["region"]} {r["n_records"]} rec')
    man = os.path.join(a.out, 'trio_fetch_manifest.tsv')
    with open(man, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]), delimiter='\t'); w.writeheader(); w.writerows(res)
    log(f'{len(res)} slices, {sum(r["n_records"] for r in res)} records -> {man}')


def _ctx(cid, preset):
    k = (cid, preset)
    if k not in G['ctx']:
        G['ctx'][k] = RE.ClaimContext(G['claims'][cid], G['ends'], G['crit'], preset=preset, mask_cache=G['mask_cache'])
    return G['ctx'][k]


def work1(j):
    ctx = _ctx(j['claim'], j['preset']); th = ctx.th; out = []; Q = []; n = n_tel = 0
    motifs = RE.telomere_motifs(th, j['platform'])
    for name, seq, rec in RE.iter_reads(j['bam'], min_len=G['min_len'], reference=G['grch38'] if j['bam'].endswith('.cram') else None):
        n += 1
        t, o, tb = RE.orient(seq, j['side'], th, motifs)
        if not (tb and tb['junk_len'] <= th['telomere_max_offset'] and tb['telo_len'] >= th['telomere_min_bp']): continue
        n_tel += 1
        rg = rec.get_tag('RG') if rec.has_tag('RG') else (name.split('/')[0] if name.startswith('m') and '/' in name else '')
        meta = dict(dataset=j['src'], platform=j['platform'], role=j['role'], sample=j['sample'], job=j['job'], slice_end=j['end'],
                    run=rg, source_url=j['url'], region=j['region'], input_contig=rec.reference_name, input_pos=rec.reference_start,
                    input_end=rec.reference_end, input_mapq=rec.mapping_quality, input_strand='-' if rec.is_reverse else '+',
                    input_has_sa=rec.has_tag('SA'), reference=j['reference'])
        r, u, t = RE.evaluate_read(ctx, name, seq, side=j['side'], meta=meta)
        r['_w'] = f"W|{j['src']}|{name}"; r['_a'] = f"A|{j['job']}|{name}"
        Q.append((r['_w'], t))
        if len(t) - r['anchor_from'] >= th['anchor_min_bp']: Q.append((r['_a'], t[r['anchor_from']:]))
        out.append((r, u))
    return j, out, Q, dict(job=j['job'], n_reads=n, n_telomeric=n_tel)


def init(claims, ends, crit, mask_cache, min_len, grch38):
    G.update(claims=claims, ends=ends, crit=crit, mask_cache=mask_cache, min_len=min_len, grch38=grch38, ctx={})


def build_jobs(claims, a, man):
    reg = {(r['source'], r['end']): r['region'] for r in man}
    fil = {(r['source'], r['end']): os.path.join(REPO, r['file']) for r in man}
    jobs = []
    for c in claims:
        tp = 'HG003' if 'PATERNAL' in c['contig'] else 'HG004'
        for src in a.sources:
            s, url, ref, preset, plat = SOURCES[src]
            role = 'transmitting_parent' if s == tp else 'other_parent'
            for end in sorted({c['R'], c['D'], c['P']}):
                if (src, end) not in reg: continue
                jobs.append(dict(claim=c['claim_id'], src=src, sample=s, url=url, reference=ref, preset=preset, platform=plat, role=role,
                                 end=end, side='L' if end.endswith('p') else 'R', bam=fil[(src, end)], region=reg[(src, end)]))
    for i, j in enumerate(jobs): j['job'] = i
    return jobs


def stage_eval(a, claims):
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    ends = RE.Chm13Ends(a.chm13_ends)
    man = list(csv.DictReader(open(os.path.join(a.out, 'trio_fetch_manifest.tsv')), delimiter='\t'))
    jobs = build_jobs(claims, a, man)
    cmap = {c['claim_id']: c for c in claims}
    mask_dir = os.path.join(REPO, 'work', 'readevidence', 'masks')
    for c in claims: RE.shared_mask_cached(ends, c['R'], th, mask_dir)
    log(f'{len(jobs)} jobs; criteria {crit["status"]} {RE.criteria_hash(crit)}')
    t0 = time.time(); jobs.sort(key=lambda j: -os.path.getsize(j['bam'])); res = []
    with mp_.get_context('fork').Pool(a.threads, initializer=init, initargs=(cmap, ends, crit, mask_dir, a.min_len, a.grch38)) as pool:
        for i, x in enumerate(pool.imap_unordered(work1, jobs, chunksize=1)):
            res.append(x)
            if i % 50 == 0: log(f'stage 1: {i + 1}/{len(jobs)} jobs, {time.time() - t0:.0f} s')
    Q = {}
    for _, out, q, _ in res:
        for k, s in q: Q.setdefault(k, s)
    log(f'genome mapping of {len(Q)} queries ({sum(len(s) for s in Q.values()) / 1e6:.0f} Mb)')
    Gm = RE.map_genome(list(Q.items()), a.genome, a.threads, sam_out=os.path.join(a.work, 'genome.sam'))
    with open(os.path.join(a.work, 'stage.pkl'), 'wb') as f:
        pickle.dump(dict(res=[(j, out, st) for j, out, _, st in res], Gm=Gm), f)
    log(f'stage done: {time.time() - t0:.0f} s')


def stage_report(a, claims):
    import pandas as pd, pysam
    crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    ends = RE.Chm13Ends(a.chm13_ends)
    cmap = {c['claim_id']: c for c in claims}
    init(cmap, ends, crit, os.path.join(REPO, 'work', 'readevidence', 'masks'), a.min_len, a.grch38)
    S = pickle.load(open(os.path.join(a.work, 'stage.pkl'), 'rb')); Gm = S['Gm']
    wf = pysam.FastaFile(a.windows); asm = {}
    for c in claims:
        asm[c['claim_id']] = RE.assembly_projection(_ctx(c['claim_id'], 'map-hifi'), wf.fetch(c['end_id']).upper()) \
            if c['end_id'] in wf.references else None
    fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.chm13_fai)}
    rows = []; stats = []
    for j, out, st in S['res']:
        ctx = _ctx(j['claim'], j['preset']); stats.append(dict(st, claim=j['claim'], source=j['src'], end=j['end'], role=j['role']))
        for r, u in out:
            RE.finalize(ctx, r, Gm.get(r['_w']), Gm.get(r['_a'])); RE.concordance(r, asm[j['claim']], th)
            c = ctx.c; r.update(D=c['D'], P=c['P'], R=c['R'], claim_pattern=c['pattern'])
            r['switch_ci_P_chm13_lo'] = chm13_coord(fai, c['P'], r.get('switch_ci_P_hi') if c['P'].endswith('q') else r.get('switch_ci_P_lo'))
            r['switch_ci_P_chm13_hi'] = chm13_coord(fai, c['P'], r.get('switch_ci_P_lo') if c['P'].endswith('q') else r.get('switch_ci_P_hi'))
            r['eligible'] = bool(r['terminal'] and c.get('switch_sub') is not None and
                                 r['read_len'] >= r['s0_read'] + c['switch_sub'] + th['anchor_min_bp'])
            rows.append(r)
    H = HM.load_ends(os.path.join(REPO, 'data', 'catalogue')).set_index('end_id')
    hg = pd.read_csv(os.path.join(a.out, 'claims.tsv'), sep='\t').set_index('claim_id') if os.path.exists(os.path.join(a.out, 'claims.tsv')) else None
    out = []
    for c in claims:
        for s in ('HG003', 'HG004'):
            rr = [r for r in rows if r['claim_id'] == c['claim_id'] and r['sample'] == s]
            # one read evaluated from several slices of the same source (R, D, P) cannot occur for primary records; guard anyway
            seen = set(); uq = []
            for r in rr:
                k = (r['dataset'], r['read'])
                if k in seen: continue
                seen.add(k); uq.append(r)
            v = RE.claim_verdict(uq, crit, same_hap_key=None)
            conc = collections.Counter(r['dataset'] for r in uq if r['verdict'] == 'conclusive' and r.get('switch_concordant'))
            ncon = sum(r['verdict'] == 'conclusive' for r in uq); nctr = sum(r['verdict'] == 'contradicts' for r in uq)
            role = 'transmitting_parent' if s == ('HG003' if 'PATERNAL' in c['contig'] else 'HG004') else 'other_parent'
            e = H.loc[c['end_id']]
            rec = dict(claim_id=c['claim_id'], end_id=c['end_id'], hg002_haplotype=e['haplotype'], parent=s, parent_role=role,
                       allele_id=e['structural_allele_id'], core_structure=e['core_structure'], donor_D=c['D'], proximal_P=c['P'],
                       recipient_R=c['R'], assembly_switch_sub=c.get('switch_sub'),
                       hg002_verdict=(hg.loc[c['claim_id'], 'verdict'] if hg is not None else None),
                       verdict=v['verdict'], n_telomeric_reads=len(uq), n_eligible=sum(r['eligible'] for r in uq),
                       n_switch_pass=sum(bool(r.get('switch_pass')) for r in uq), n_conclusive=ncon,
                       n_conclusive_concordant=v['n_conclusive_concordant'], n_conclusive_discordant=v['n_conclusive_discordant'],
                       n_consistent=v['n_consistent'], n_contradicts=nctr,
                       conclusive_fraction=round(ncon / (ncon + nctr), 3) if ncon + nctr else None,
                       datasets_conclusive_concordant=';'.join(f'{d}:{n}' for d, n in sorted(conc.items())),
                       reads_by_source=';'.join(f'{d}:{n}' for d, n in sorted(collections.Counter(r['dataset'] for r in uq).items())),
                       reason=HM.reason_for(v, uq, c, asm[c['claim_id']]) if v['verdict'] != 'confirmed' else '',
                       criteria_sha=RE.criteria_hash(crit))
            out.append(rec)
    T = pd.DataFrame(out); T.to_csv(os.path.join(a.out, 'trio_claims.tsv'), sep='\t', index=False)
    keep = [r for r in rows if r['verdict'] in ('conclusive', 'contradicts')]
    cols = ['claim_id', 'sample', 'role', 'read', 'dataset', 'platform', 'source_url', 'region', 'reference', 'slice_end', 'run', 'read_len',
            'orientation', 'input_contig', 'input_pos', 'input_end', 'input_mapq', 'input_strand', 'input_has_sa', 'telo_offset',
            'telo_len', 's0_read', 'n_units', 'segments', 'seg_D_units', 'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'dR_units',
            'dR_purity', 'log10_lr_vs_allP', 'log10_lr_vs_allD', 'switch_read', 'switch_ci_read_lo', 'switch_ci_read_hi',
            'switch_ci_sub_lo', 'switch_ci_sub_hi', 'P', 'switch_ci_P_lo', 'switch_ci_P_hi', 'switch_ci_P_chm13_lo', 'switch_ci_P_chm13_hi',
            'switch_concordant', 'concordance_ref', 'donor_interval_units', 'donor_interval_P_purity', 'anchor_chrom', 'anchor_tip_start',
            'anchor_tip_end', 'anchor_mapq', 'anchor_len', 'anchor_unique_bp', 'anchor_tier', 'internal_telomere_bp', 'adapter_hits',
            'foldback', 'subtel_aligned_frac', 'eligible', 'verdict', 'D', 'R']
    with open(os.path.join(a.out, 'trio_reads.tsv'), 'w') as f: RE.write_tsv(keep, f, cols)
    with gzip.open(os.path.join(a.work, 'trio_reads_all.tsv.gz'), 'wt') as f:
        RE.write_tsv(rows, f, cols + ['terminal', 'switch_pass', 'anchor_pass', 'integrity_pass', 'genome_alignments'])
    pd.DataFrame(stats).to_csv(os.path.join(a.work, 'job_stats.tsv'), sep='\t', index=False)
    # evidence BAM (CHM13 v2.0 alignments of conclusive and contradicting parent reads)
    od = os.path.join(a.work, 'evidence_bam'); os.makedirs(od, exist_ok=True)
    by_w = collections.defaultdict(list)
    for r in keep: by_w[r['_w']].append(r)
    rgs = {r['dataset']: dict(ID=r['dataset'], SM=r['sample'], PL='ONT' if r['platform'].startswith('ONT') else 'PACBIO',
                              DS=f"{r['platform']}; source={r['source_url']}") for r in keep}
    if keep:
        with pysam.AlignmentFile(os.path.join(a.work, 'genome.sam')) as fi:
            hd = fi.header.to_dict(); hd['RG'] = list(rgs.values())
            hd.setdefault('PG', []).append(dict(ID='readevidence', PN='subtel-atlas', CL='evidence/scripts/trio_reads.py'))
            hd = HM.clean_header(hd)
            tmp = os.path.join(od, 'trio_evidence.chm13v2.unsorted.bam')
            with pysam.AlignmentFile(tmp, 'wb', header=hd) as fo:
                for rec in fi:
                    if rec.query_name not in by_w: continue
                    rs = by_w[rec.query_name]; rec.query_name = rs[0]['read']; rec.set_tag('RG', rs[0]['dataset'])
                    rec.set_tag('XC', ','.join(sorted({r['claim_id'] for r in rs}))); rec.set_tag('XV', ','.join(sorted({r['verdict'] for r in rs})))
                    fo.write(rec)
        o = os.path.join(od, 'trio_evidence.chm13v2.bam'); pysam.sort('-o', o, tmp); pysam.index(o); os.remove(tmp)
    log(f'{len(rows)} telomeric parent reads evaluated; {len(keep)} conclusive/contradicting; ' +
        ', '.join(f'{k}={n}' for k, n in T[T.parent_role == "transmitting_parent"].verdict.value_counts().items()) + ' (transmitting parent)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['fetch', 'eval', 'report', 'all'])
    ap.add_argument('--sources', default=','.join(SOURCES))
    ap.add_argument('--criteria', default=os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    ap.add_argument('--chm13-ends'); ap.add_argument('--chm13-fai'); ap.add_argument('--genome')
    ap.add_argument('--grch38', help='GRCh38 no-alt analysis set FASTA (faidx) for decoding the ONT R10 CRAMs')
    ap.add_argument('--windows', default=os.path.join(REPO, 'work', 'windows', 'ends.fa'))
    ap.add_argument('--work', default=os.path.join(REPO, 'work', 'evidence', 'trio'))
    ap.add_argument('--out', default=os.path.join(REPO, 'evidence', 'hg002'))
    ap.add_argument('--min-len', type=int, default=5000)
    ap.add_argument('-t', '--threads', type=int, default=4)
    ap.add_argument('-j', '--jobs', type=int, default=4)
    a = ap.parse_args(); a.sources = a.sources.split(',')
    if a.grch38: a.grch38 = os.path.abspath(a.grch38)
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    cat = os.path.join(REPO, 'data', 'catalogue')
    H = HM.load_ends(cat); ids = H[H.nativity_class != 'native'].end_id.tolist()
    claims = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=ids)
    log(f'{len(claims)} claims; sources: {",".join(a.sources)}')
    if a.stage in ('fetch', 'all'): stage_fetch(a, claims)
    if a.stage in ('eval', 'all'): stage_eval(a, claims)
    if a.stage in ('report', 'all'): stage_report(a, claims)


if __name__ == '__main__':
    main()
