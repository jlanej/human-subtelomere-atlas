#!/usr/bin/env python3
"""Build the reviewer-facing read-evidence resource from the frozen-criteria runs (criteria v1, hash a89c5eaf023ed032).

Inputs (all already produced; nothing here re-decides a v1 verdict):
  evidence/hg002/claims.tsv, evidence/hg002/trio_claims.tsv, evidence/population/claims.tsv, controls.tsv
  work/evidence/hg002/reads_all.tsv.gz, work/evidence/trio/trio_reads_all.tsv.gz, work/evidence/population/reads_full.tsv
  the read slices fetched from the public alignments (work/evidence/*/slices) and the end windows work/windows/ends.fa
  (identical to windows cut from the public assemblies; SHA-256 in data/catalogue/ends.tsv).

Stages:
  collect   conclusive reads -> FASTQ (work/evidence/release_build/fq) and own-assembly placement against every end window of the
            same person (HG002 / HG003 / HG004: the 92 HG002 Q100 windows; 1KG-ONT individuals: their Release 2 windows)
  tables    evidence/claims.tsv, evidence/reads.tsv.gz, evidence/controls.tsv (unified; definitive tier, evidence_lib.py)
  bams      per-claim CRAM + BAM on CHM13 v2.0 and BAM on the person's own end window(s) (+ window FASTA) for every claim with a
            v1-conclusive read; showcase subset under evidence/showcase/; release tarballs + SHA256SUMS in work/evidence/release/
  igv       IGV sessions for the showcase claims
Usage: python evidence/scripts/build_evidence_set.py collect|tables|bams|igv|all [-t 6]
"""
import argparse, os, sys, csv, gzip, json, collections, subprocess, tempfile, shutil, time
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_lib as L  # noqa: E402

W = os.path.join(REPO, 'work', 'evidence')
B = os.path.join(W, 'release_build')
CRIT_SHA = 'a89c5eaf023ed032'


def log(*a): print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


def ends_table():
    import pandas as pd
    return pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False)


def load_read_tables():
    """carrier/claim read rows of the three cohorts as pandas frames (all verdicts)"""
    import pandas as pd
    H = pd.read_csv(os.path.join(W, 'hg002', 'reads_all.tsv.gz'), sep='\t', low_memory=False)
    T = pd.read_csv(os.path.join(W, 'trio', 'trio_reads_all.tsv.gz'), sep='\t', low_memory=False)
    P = pd.read_csv(os.path.join(W, 'population', 'reads_full.tsv'), sep='\t', low_memory=False)
    return H, T, P


# ------------------------------------------------------------------------------------------------------------------ collect
def _slice_files_hg002(row):
    contig, side = row['source_end'].split('#')[2], row['source_end'].split('#')[3]
    return [os.path.join(W, 'hg002', 'slices', row['dataset'], f'{contig}_{side}.bam')]


def _slice_files_trio(row):
    d = os.path.join(W, 'trio', 'slices', row['dataset'])
    return [os.path.join(d, x) for x in sorted(os.listdir(d)) if x.endswith(('.bam', '.cram')) and x.split('.')[0].split('_')[0] == row['slice_end']]


def _pop_slices():
    P = json.load(open(os.path.join(W, 'population', 'plan.json')))
    return {(x['sample'], x['region']): os.path.join(W, 'population', x['file']) for x in P['slices']}


def _read_records(files, names, reference=None):
    """{name: (seq, qual_string, rec_strand)} for the primary records of `names` in the slice files (stored orientation)"""
    import pysam
    out = {}
    for fn in files:
        if not os.path.exists(fn): continue
        kw = dict(reference_filename=reference) if reference and fn.endswith('.cram') else {}
        with pysam.AlignmentFile(fn, **kw) as fh:
            for r in fh.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_name not in names or r.query_name in out: continue
                q = ''.join(chr(min(v, 93) + 33) for v in r.query_qualities) if r.query_qualities is not None else 'I' * r.query_length
                out[r.query_name] = (r.query_sequence, q)
    return out


def collect(a):
    import pandas as pd, pysam
    os.makedirs(os.path.join(B, 'fq'), exist_ok=True)
    H, T, P = load_read_tables(); E = ends_table()
    wf = pysam.FastaFile(os.path.join(REPO, 'work', 'windows', 'ends.fa'))
    want = []   # (person, cohort, claim_id, dataset, read, files)
    for r in H[(H.verdict == 'conclusive') & (H.role == 'carrier')].to_dict('records'):
        want.append(('HG002', 'HG002_Q100', r['claim_id'], r['dataset'], r['read'], _slice_files_hg002(r)))
    for r in T[T.verdict == 'conclusive'].to_dict('records'):
        want.append((r['sample'], 'HG002_trio', r['claim_id'], r['dataset'], r['read'], _slice_files_trio(r)))
    ps = _pop_slices()
    for r in P[(P.verdict == 'conclusive') & (P.role == 'carrier')].to_dict('records'):
        want.append((r['sample'], '1KG_ONT', r['claim_id'], r['dataset'], r['read'], [ps[(r['sample'], r['slice_region'])]]))
    log(f'{len(want)} conclusive read rows to collect')
    # sequences, per (person, dataset, files)
    groups = collections.defaultdict(set)
    for p, c, cid, ds, rd, files in want: groups[(p, ds, tuple(files))].add(rd)
    seqs = collections.defaultdict(dict)   # person -> read -> (seq, qual)
    for (p, ds, files), names in groups.items():
        got = _read_records(list(files), names, reference=a.chm13)
        miss = names - set(got)
        if miss: raise SystemExit(f'{p} {ds}: {len(miss)} reads not found in {files}')
        for n, v in got.items(): seqs[p].setdefault(n, v)
    for p, d in seqs.items():
        with gzip.open(os.path.join(B, 'fq', f'{p}.fq.gz'), 'wt') as fh:
            for n, (s, q) in sorted(d.items()): fh.write(f'@{n}\n{s}\n+\n{q}\n')
    # placement windows per person
    hg = E[E.set == 'HG002_Q100'].end_id.tolist()
    win_ids = {p: (hg if p in ('HG002', 'HG003', 'HG004') else E[(E.set == 'HPRC_R2') & (E['sample'] == p)].end_id.tolist())
               for p in seqs}
    rows = []
    for p in sorted(seqs):
        ids = [e for e in win_ids[p] if e in wf.references]
        pl = L.placement({n: s for n, (s, q) in seqs[p].items()}, {e: wf.fetch(e) for e in ids})
        for n, v in pl.items():
            rows.append(dict(person=p, read=n, n_windows=len(ids), placement=json.dumps(v)))
        log(f'{p}: {len(seqs[p])} reads placed on {len(ids)} windows')
    pd.DataFrame(rows).to_csv(os.path.join(B, 'placement.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(person=p, cohort=c, claim_id=cid, dataset=ds, read=rd) for p, c, cid, ds, rd, _ in want]).to_csv(
        os.path.join(B, 'conclusive_rows.tsv'), sep='\t', index=False)


# ------------------------------------------------------------------------------------------------------------------ tables
CLAIM_COLS = ['claim_id', 'slug', 'cohort', 'individual', 'haplotype', 'arm', 'allele_id', 'nativity_class', 'core_structure',
              'donor_D', 'proximal_P', 'recipient_R', 'pattern', 'claimed_end_ids', 'model_end', 'other_haplotype',
              'assembly_switch_sub', 'assembly_switch_chm13_P', 'asm_projection_switch_ci_sub', 'n_reads_tested',
              'n_conclusive', 'n_conclusive_concordant', 'n_conclusive_discordant', 'n_consistent', 'n_contradicts',
              'n_concordant_unique_anchor', 'n_concordant_own_on_claimed_end', 'n_definitive', 'n_definitive_molecules',
              'n_definitive_datasets', 'definitive_datasets', 'reads_by_tier', 'conclusive_concordant_by_dataset', 'v1_verdict',
              'v1_class', 'definitive_verdict', 'diagnostics_note', 'source', 'criteria_sha']
READ_COLS = ['read_id', 'claim_id', 'cohort', 'role', 'individual', 'dataset', 'platform', 'derived', 'counted', 'duplicate_of',
             'molecule', 'source_url', 'index_url', 'source_reference', 'region', 'orient_side', 'read_group', 'read_len', 'input_contig',
             'input_pos', 'input_strand', 'input_mapq', 'input_has_sa', 'telo_offset', 'telo_len', 's0_read', 'n_units',
             'seg_D_units', 'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'dR_units', 'dR_purity', 'log10_lr_vs_allP',
             'log10_lr_vs_allD', 'switch_read', 'switch_ci_read_lo', 'switch_ci_read_hi', 'switch_ci_sub_lo', 'switch_ci_sub_hi',
             'switch_P', 'switch_ci_P_tip_lo', 'switch_ci_P_tip_hi', 'switch_ci_chm13', 'switch_concordant', 'concordance_ref',
             'anchor_chrom', 'anchor_chm13', 'anchor_mapq', 'anchor_len', 'anchor_unique_bp', 'anchor_tier', 'internal_telomere_bp',
             'adapter_hits', 'foldback', 'subtel_aligned_frac', 'v1_verdict', 'own_assembly', 'own_best_end', 'own_best_identity',
             'own_claimed_identity', 'own_claimed_matches', 'own_runner_up_end', 'own_runner_up_matches', 'own_margin_matches',
             'own_on_claimed_end', 'tier', 'samtools_command']


def _slug(cid):
    return cid.replace('#', '_').replace('|', '__').replace('>', '_')


def _hap_of(end_id):
    c = end_id.split('#')[2]
    return 'maternal' if 'MATERNAL' in c else 'paternal' if 'PATERNAL' in c else end_id.split('#')[1]


def _fai(a):
    return {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.chm13 + '.fai')}


def _read_row(r, cohort, claim_id, claimed, fai, place, own_label, counted=True):
    """unified per-read row from a pipeline read dict"""
    g = lambda k: (None if (r.get(k) is None or (isinstance(r.get(k), float) and r.get(k) != r.get(k))) else r.get(k))  # noqa: E731
    url = g('source_url') or g('cram_url'); region = g('region') or g('slice_region')
    cram = isinstance(url, str) and url.endswith('.cram')
    ref = g('reference') or ('CHM13v2.0' if cohort != 'HG002_Q100' else 'HG002 Q100 v1.1')
    P = g('P'); R = g('R')
    sw = L.chm13_interval(fai, P, g('switch_ci_P_lo'), g('switch_ci_P_hi')) if P and g('switch_ci_P_lo') is not None else None
    an = L.chm13_interval(fai, R, g('anchor_tip_start'), g('anchor_tip_end')) if R and g('anchor_tip_start') is not None else None
    pf = L.placement_fields(place.get((r.get('_person'), r['read'])), claimed) if (r.get('_person'), r['read']) in place else {}
    se = g('slice_end'); side = (se[-1] == 'p' and 'L' or 'R') if isinstance(se, str) else \
        (g('source_end').split('#')[3] if isinstance(g('source_end'), str) else None)
    row = dict(orient_side=side, read_id=r['read'], claim_id=claim_id, cohort=cohort, role=r.get('role'), individual=r.get('_person'),
               dataset=g('dataset'), platform=g('platform'), derived=bool(L.b(g('derived'))), counted=counted,
               duplicate_of=g('duplicate_of') or '', molecule=g('molecule') or r['read'], source_url=url,
               index_url=L.index_url(url) if url else None, source_reference=ref, region=region,
               read_group=g('read_group') or g('run') or '', read_len=g('read_len'), input_contig=g('input_contig'),
               input_pos=g('input_pos'), input_strand=g('input_strand'), input_mapq=g('input_mapq'), input_has_sa=g('input_has_sa'),
               switch_P=P, switch_ci_P_tip_lo=g('switch_ci_P_lo'), switch_ci_P_tip_hi=g('switch_ci_P_hi'),
               switch_ci_chm13=f'{sw[0]}:{sw[1]}-{sw[2]}' if sw else '', anchor_chm13=f'{an[0]}:{an[1]}-{an[2]}' if an else '',
               v1_verdict=g('verdict'), own_assembly=own_label if pf else '',
               samtools_command=L.samtools_cmd(url, region, r['read'], cram=cram,
                                               ref='chm13v2.0.fa.gz' if ref == 'CHM13v2.0' else 'GRCh38_no_alt_analysis_set.fa')
               if url else '')
    for k in ('telo_offset', 'telo_len', 's0_read', 'n_units', 'seg_D_units', 'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'dR_units',
              'dR_purity', 'log10_lr_vs_allP', 'log10_lr_vs_allD', 'switch_read', 'switch_ci_read_lo', 'switch_ci_read_hi',
              'switch_ci_sub_lo', 'switch_ci_sub_hi', 'switch_concordant', 'concordance_ref', 'anchor_chrom', 'anchor_mapq',
              'anchor_len', 'anchor_unique_bp', 'anchor_tier', 'internal_telomere_bp', 'adapter_hits', 'foldback', 'subtel_aligned_frac'):
        row[k] = g(k)
    row.update(pf)
    row['tier'] = L.read_tier(row) if counted or row['v1_verdict'] != 'conclusive' else 'not_counted_' + L.read_tier(row)
    return row


def _tier_counts(rows):
    c = collections.Counter(r['tier'] for r in rows if r['counted'])
    return ';'.join(f'{t}:{c[t]}' for t in L.TIER_ORDER if c.get(t))


def tables(a):
    import pandas as pd
    H, T, P = load_read_tables(); E = ends_table().set_index('end_id'); fai = _fai(a)
    PL = pd.read_csv(os.path.join(B, 'placement.tsv'), sep='\t')
    place = {(p, r): json.loads(v) for p, r, v in zip(PL.person, PL.read, PL.placement)}
    claims, reads, controls = [], [], []
    keep_v = ('conclusive', 'consistent', 'contradicts')
    # ---------------- HG002 Q100
    HC = pd.read_csv(os.path.join(REPO, 'evidence', 'hg002', 'claims.tsv'), sep='\t')
    H['_person'] = 'HG002'
    for c in HC.to_dict('records'):
        cid = c['claim_id']; cr = H[(H.claim_id == cid) & (H.role == 'carrier')]
        rows = []
        for r in cr.to_dict('records'):
            counted = (not L.b(r['derived'])) and not (isinstance(r['duplicate_of'], str) and r['duplicate_of'])
            if r['verdict'] in keep_v:
                rows.append(_read_row(r, 'HG002_Q100', cid, [c['end_id']], fai, place, 'HG002 Q100 v1.1 end windows (92)', counted))
        # HG002 placement: the Q100 primary alignment is on the claimed contig by construction of the read set; definitive also
        # needs the window placement on the claimed end
        reads += rows
        prim = [r for r in rows if r['counted']]
        conc = [r for r in prim if r['v1_verdict'] == 'conclusive' and L.b(r['switch_concordant'])]
        dfn = [r for r in prim if r['tier'] == 'definitive']
        dv, nm, ds = L.definitive_claim(c['verdict'], dfn, molecule_key=lambda r: r['molecule'])
        assert len(conc) == c['n_conclusive_concordant'], (cid, len(conc), c['n_conclusive_concordant'])
        e = E.loc[c['end_id']]
        claims.append(dict(
            claim_id=cid, slug=_slug(cid), cohort='HG002_Q100', individual='HG002', haplotype=_hap_of(cid), arm=e['home_arm'],
            allele_id=c['allele_id'], nativity_class=c['nativity_class'], core_structure=c['core_structure'], donor_D=c['donor_D'],
            proximal_P=c['proximal_P'], recipient_R=c['recipient_R'], pattern=c['pattern'], claimed_end_ids=c['end_id'],
            model_end=c['end_id'], other_haplotype='', assembly_switch_sub=c['assembly_switch_sub'],
            assembly_switch_chm13_P=f"{c['proximal_P'][:-1]}:{int(c['chm13_p_switch_coord'])}" if L.f(c['chm13_p_switch_coord']) else '',
            asm_projection_switch_ci_sub=f"{int(c['asm_projection_sub_lo'])}-{int(c['asm_projection_sub_hi'])}" if L.f(c['asm_projection_sub_lo']) is not None else '',
            n_reads_tested=c['n_reads'], n_conclusive=c['n_conclusive'], n_conclusive_concordant=c['n_conclusive_concordant'],
            n_conclusive_discordant=c['n_conclusive_discordant'], n_consistent=c['n_consistent'], n_contradicts=c['n_contradicts'],
            n_concordant_unique_anchor=sum((L.f(r['anchor_unique_bp']) or 0) >= L.UNIQUE_TIER_BP for r in conc),
            n_concordant_own_on_claimed_end=sum(L.b(r.get('own_on_claimed_end')) for r in conc), n_definitive=len(dfn),
            n_definitive_molecules=nm, n_definitive_datasets=len(ds), definitive_datasets=','.join(ds), reads_by_tier=_tier_counts(rows),
            conclusive_concordant_by_dataset=c['datasets_conclusive_concordant'] if isinstance(c['datasets_conclusive_concordant'], str) else '',
            v1_verdict=c['verdict'], v1_class=c['verdict'], definitive_verdict=dv,
            diagnostics_note=c['reason'] if isinstance(c['reason'], str) else '',
            source='HG002 Q100 v1.1 public read alignments (16 files; 12 independent datasets after molecule de-duplication)',
            criteria_sha=c['criteria_sha']))
        # controls: other-haplotype strictly native reads at the same arm
        ctl = H[(H.claim_id == cid) & (H.role == 'hg002_otherhap_native')]
        for r in ctl[ctl.verdict.isin(keep_v)].to_dict('records'):
            reads.append(_read_row(r, 'HG002_Q100', f'control|HG002_otherhap|{cid}', [c['end_id']], fai, {}, '', True))
        if len(ctl):
            el = ctl[(ctl.terminal == True) & (ctl.anchor_pass == True) & (ctl.integrity_pass == True)]  # noqa: E712
            controls.append(dict(control_id=f'control|HG002_otherhap|{cid}', cohort='HG002_Q100', individual='HG002',
                                 arm=e['home_arm'], allele_model=cid, control_type='HG002 other haplotype, strictly native',
                                 n_reads=len(ctl), n_eligible=len(el), n_conclusive=int((ctl.verdict == 'conclusive').sum()),
                                 n_contradicts=int((ctl.verdict == 'contradicts').sum()), n_consistent=int((ctl.verdict == 'consistent').sum())))
    # ---------------- trio (HG003 / HG004), placement against the HG002 Q100 windows (the transmitted haplotype is in HG002)
    TCl = pd.read_csv(os.path.join(REPO, 'evidence', 'hg002', 'trio_claims.tsv'), sep='\t')
    T['_person'] = T['sample']
    for c in TCl.to_dict('records'):
        cid = f"{c['parent']}|{c['claim_id']}"
        tr = T[(T.claim_id == c['claim_id']) & (T['sample'] == c['parent'])].drop_duplicates(['dataset', 'read'])
        rows = [_read_row(r, 'HG002_trio', cid, [c['end_id']], fai, place,
                          'HG002 Q100 v1.1 end windows (transmitted-haplotype proxy)', True)
                for r in tr.to_dict('records') if r['verdict'] in keep_v]
        for r in rows: r['role'] = c['parent_role']
        reads += rows
        conc = [r for r in rows if r['v1_verdict'] == 'conclusive' and L.b(r['switch_concordant'])]
        assert len(conc) == c['n_conclusive_concordant'], cid
        dfn = [r for r in rows if r['tier'] == 'definitive']
        dv, nm, ds = L.definitive_claim(c['verdict'], dfn)
        e = E.loc[c['end_id']]; hc = HC.set_index('claim_id').loc[c['claim_id']]
        claims.append(dict(
            claim_id=cid, slug=_slug(cid), cohort='HG002_trio', individual=c['parent'],
            haplotype=f"{c['parent_role']} (HG002 {_hap_of(c['end_id'])})", arm=e['home_arm'], allele_id=c['allele_id'],
            nativity_class=hc['nativity_class'], core_structure=c['core_structure'], donor_D=c['donor_D'], proximal_P=c['proximal_P'],
            recipient_R=c['recipient_R'], pattern=hc['pattern'], claimed_end_ids=c['end_id'], model_end=c['end_id'],
            other_haplotype='', assembly_switch_sub=c['assembly_switch_sub'],
            assembly_switch_chm13_P=f"{hc['proximal_P'][:-1]}:{int(hc['chm13_p_switch_coord'])}" if L.f(hc['chm13_p_switch_coord']) else '',
            asm_projection_switch_ci_sub=f"{int(hc['asm_projection_sub_lo'])}-{int(hc['asm_projection_sub_hi'])}" if L.f(hc['asm_projection_sub_lo']) is not None else '',
            n_reads_tested=c['n_telomeric_reads'], n_conclusive=c['n_conclusive'], n_conclusive_concordant=c['n_conclusive_concordant'],
            n_conclusive_discordant=c['n_conclusive_discordant'], n_consistent=c['n_consistent'], n_contradicts=c['n_contradicts'],
            n_concordant_unique_anchor=sum((L.f(r['anchor_unique_bp']) or 0) >= L.UNIQUE_TIER_BP for r in conc),
            n_concordant_own_on_claimed_end=sum(L.b(r.get('own_on_claimed_end')) for r in conc), n_definitive=len(dfn),
            n_definitive_molecules=nm, n_definitive_datasets=len(ds), definitive_datasets=','.join(ds), reads_by_tier=_tier_counts(rows),
            conclusive_concordant_by_dataset=c['datasets_conclusive_concordant'] if isinstance(c['datasets_conclusive_concordant'], str) else '',
            v1_verdict=c['verdict'], v1_class=c['parent_role'], definitive_verdict=dv,
            diagnostics_note=(c['reason'] if isinstance(c['reason'], str) else '') + ('' if c['parent_role'] == 'transmitting_parent' else
                                                                                        ' (other parent: reported only)'),
            source='GIAB HG003/HG004 public long-read alignments (HiFi on CHM13 v2.0; ONT on GRCh38)', criteria_sha=c['criteria_sha']))
    # ---------------- 1KG-ONT
    PC = pd.read_csv(os.path.join(REPO, 'evidence', 'population', 'claims.tsv'), sep='\t')
    AN = pd.read_csv(os.path.join(REPO, 'evidence', 'population', 'diagnostics', 'claim_annotations.tsv'), sep='\t').set_index('claim_id')
    P['_person'] = P['sample']
    for c in PC.to_dict('records'):
        cid = c['claim_id']; R = c['R'] if isinstance(c['R'], str) else None
        cr = P[(P.claim_id == cid) & (P.role == 'carrier')]
        cr = cr[(cr.slice_end == R) | (cr.verdict == 'conclusive')] if R else cr.iloc[0:0]
        claimed = c['end_ids'].split(',')
        rows = [_read_row(r, '1KG_ONT', cid, claimed, fai, place, f"{c['sample']} HPRC Release 2 end windows", True)
                for r in cr.to_dict('records') if r['verdict'] in keep_v]
        reads += rows
        conc = [r for r in rows if r['v1_verdict'] == 'conclusive' and L.b(r['switch_concordant'])]
        assert len(conc) == c['n_conclusive_concordant'], cid
        dfn = [r for r in rows if r['tier'] == 'definitive']
        dv, nm, ds = L.definitive_claim(c['verdict'], dfn)
        haps = sorted({x.split('#')[1] for x in claimed})
        note = ''
        if cid in AN.index:
            x = AN.loc[cid]; note = '; '.join(str(v) for v in (x['diagnostic_notes'], x['recommendation']) if isinstance(v, str) and v)
        psw = L.chm13_coord(fai, c['P'], c['assembly_p_switch']) if isinstance(c['P'], str) and L.f(c['assembly_p_switch']) is not None else None
        claims.append(dict(
            claim_id=cid, slug=_slug(cid), cohort='1KG_ONT', individual=c['sample'],
            haplotype='both (homozygous)' if len(haps) == 2 else f'hap{haps[0]}', arm=c['arm'], allele_id=c['allele'],
            nativity_class=c['nativity_class'], core_structure=E.loc[c['model_end'], 'core_structure'] if isinstance(c['model_end'], str) else '',
            donor_D=c['D'], proximal_P=c['P'], recipient_R=c['R'], pattern=c['pattern'], claimed_end_ids=c['end_ids'],
            model_end=c['model_end'], other_haplotype=c['other_haplotype'], assembly_switch_sub=c['assembly_switch_sub'],
            assembly_switch_chm13_P=f"{c['P'][:-1]}:{psw}" if psw else '',
            asm_projection_switch_ci_sub=c['asm_projection_switch_ci_sub'] if isinstance(c['asm_projection_switch_ci_sub'], str) else '',
            n_reads_tested=c['n_reads'], n_conclusive=c['n_conclusive_concordant'] + c['n_conclusive_discordant'],
            n_conclusive_concordant=c['n_conclusive_concordant'], n_conclusive_discordant=c['n_conclusive_discordant'],
            n_consistent=c['n_consistent'], n_contradicts=c['n_contradicts'],
            n_concordant_unique_anchor=sum((L.f(r['anchor_unique_bp']) or 0) >= L.UNIQUE_TIER_BP for r in conc),
            n_concordant_own_on_claimed_end=sum(L.b(r.get('own_on_claimed_end')) for r in conc), n_definitive=len(dfn),
            n_definitive_molecules=nm, n_definitive_datasets=len(ds), definitive_datasets=','.join(ds), reads_by_tier=_tier_counts(rows),
            conclusive_concordant_by_dataset=f"1KG_ONT_{c['sample']}:{c['n_conclusive_concordant']}" if c['n_conclusive_concordant'] else '',
            v1_verdict=c['verdict'], v1_class=c['verdict_class'] if c['verdict'] in ('confirmed', 'supported') else
            f"{c['verdict_class']}:{c['assessability']}", definitive_verdict=dv, diagnostics_note=note,
            source=L.ONT1KG.format(s=c['sample']), criteria_sha=c['criteria_sha']))
    K = pd.read_csv(os.path.join(REPO, 'evidence', 'population', 'controls.tsv'), sep='\t')
    pc = P[P.role == 'control_native_strict']
    for k in K.to_dict('records'):
        kr = pc[pc.claim_id == k['control_id']]
        el = kr[(kr.terminal == True) & (kr.anchor_pass == True) & (kr.integrity_pass == True)]  # noqa: E712
        for r in kr[kr.verdict.isin(keep_v) | kr.index.isin(el.index)].to_dict('records'):
            reads.append(_read_row(r, '1KG_ONT', 'control|' + k['control_id'], [], fai, {}, '', True))
        controls.append(dict(control_id='control|' + k['control_id'], cohort='1KG_ONT', individual=k['sample'], arm=k['arm'],
                             allele_model=k['allele_model'], control_type='1KG-ONT individual strictly native in both haplotypes',
                             n_reads=k['n_reads'], n_eligible=k['n_eligible'], n_conclusive=k['n_conclusive'],
                             n_contradicts=k['n_contradicts'], n_consistent=k['n_consistent']))
    C = pd.DataFrame(claims)[CLAIM_COLS]; Rr = pd.DataFrame(reads)[READ_COLS]; K2 = pd.DataFrame(controls)
    for col in ('assembly_switch_sub',):
        C[col] = C[col].map(lambda x: '' if L.f(x) is None else int(float(x)))
    C.to_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', index=False)
    Rr.to_csv(os.path.join(REPO, 'evidence', 'reads.tsv.gz'), sep='\t', index=False, compression='gzip')
    K2.to_csv(os.path.join(REPO, 'evidence', 'controls.tsv'), sep='\t', index=False)
    log(f'claims {len(C)}, reads {len(Rr)}, controls {len(K2)}')
    log(C.groupby(['cohort', 'definitive_verdict']).size().to_string())


# ------------------------------------------------------------------------------------------------------------------ bams
SHOWCASE = ['HG002#2#chr3_MATERNAL#R', 'HG002#2#chr9_MATERNAL#L', 'HG002#2#chr19_MATERNAL#R', 'HG002#1#chr17_PATERNAL#L',
            'HG004|HG002#2#chr9_MATERNAL#L', 'HG01258|20q|SA-20q-01', 'HG01258|17p|SA-17p-02', 'NA20346|2q|SA-2q-04']
SHOWCASE_MAX_READS = 6
DS_PL = lambda p: 'PACBIO' if str(p).startswith('HiFi') or 'hifi' in str(p).lower() else 'ONT'  # noqa: E731


def _mm2(ref, fq, out_sam, threads, preset='map-ont', soft_supp=True):
    """soft_supp: -Y (full SEQ in supplementary records); off for the small own-window references, where reads longer than the
    window give many supplementary records"""
    with open(out_sam, 'w') as so:
        subprocess.run(['minimap2', '-ax', preset] + (['-Y'] if soft_supp else []) + ['--MD', '--secondary=no', '-t', str(threads), ref, fq], stdout=so,
                       stderr=subprocess.DEVNULL, check=True)


def _mm2_version():
    return subprocess.run(['minimap2', '--version'], capture_output=True, text=True).stdout.strip()


def _header(sq, rgs, kind, claim_id, ref_label, y=True):
    return dict(HD=dict(VN='1.6', SO='coordinate'), SQ=sq, RG=rgs,
                PG=[dict(ID='minimap2', PN='minimap2', VN=_mm2_version(),
                         CL=f"minimap2 -ax map-ont {'-Y ' if y else ''}--MD --secondary=no {ref_label} <conclusive_reads.fq>")],
                CO=[f'subtel-atlas read-level evidence, criteria v1 ({CRIT_SHA}); claim {claim_id}; {kind}',
                    'records: v1-conclusive reads of the claim; tags XC=claim, XT=tier (definitive etc.), XV=v1 verdict, '
                    'XW=read switch 95% interval (tip-frame read coordinates), XU=anchor unique bp, XO=own-assembly best end, '
                    'XN=counted towards the claim (0 = derived or duplicate molecule)',
                    'read sequences as in the public source alignment (RG DS); base modification tags removed'])


def _write_sorted(hdr, recs_sam, out_bam):
    import pysam
    tmp = out_bam + '.unsorted.bam'
    with pysam.AlignmentFile(tmp, 'wb', header=hdr) as fo:
        for s in recs_sam:
            fo.write(pysam.AlignedSegment.fromstring(s, fo.header))
    pysam.sort('--no-PG', '-o', out_bam, tmp); pysam.index(out_bam); os.remove(tmp)


def _to_cram(bam, cram, ref):
    """CRAM (reference not embedded; @SQ M5 + public UR). The CRAM writer replaces @SQ UR by the local -T path, so the BAM header
    is restored afterwards (samtools reheader) - published files carry no local paths."""
    subprocess.run(['samtools', 'view', '--no-PG', '-C', '-T', ref, '--output-fmt-option', 'embed_ref=0', '-o', cram + '.tmp', bam], check=True)
    hdr = subprocess.run(['samtools', 'view', '-H', '--no-PG', bam], capture_output=True, text=True, check=True).stdout
    with open(cram + '.hdr.sam', 'w') as fh: fh.write(hdr)
    with open(cram, 'wb') as fo:
        subprocess.run(['samtools', 'reheader', '--no-PG', cram + '.hdr.sam', cram + '.tmp'], stdout=fo, check=True)
    os.remove(cram + '.tmp'); os.remove(cram + '.hdr.sam')
    subprocess.run(['samtools', 'index', cram], check=True)


def fixcram(a):
    """re-apply the public header to existing CRAMs (claim directories and showcase), then rebuild sums and release"""
    n = 0
    for root in (os.path.join(B, 'claims'), os.path.join(REPO, 'evidence', 'showcase')):
        for dp, _, fs in os.walk(root):
            for f in fs:
                if f.endswith('.chm13.cram'):
                    cram = os.path.join(dp, f); _to_cram(cram[:-5] + '.bam', cram, a.chm13); n += 1
    log(f'{n} CRAMs rewritten')
    _sums(os.path.join(REPO, 'evidence', 'showcase'), os.path.join(REPO, 'evidence', 'showcase', 'SHA256SUMS'))
    release(a)


def _md5(s):
    import hashlib
    return hashlib.md5(s.upper().encode()).hexdigest()


def bams(a):
    import pandas as pd, pysam
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False)
    Rr = pd.read_csv(os.path.join(REPO, 'evidence', 'reads.tsv.gz'), sep='\t', low_memory=False)
    E = ends_table(); Ei = E.set_index('end_id')
    conc = Rr[(Rr.v1_verdict == 'conclusive') & ~Rr.claim_id.str.startswith('control|')].copy()
    log(f'{len(conc)} conclusive rows in {conc.claim_id.nunique()} claims')
    seqs = {}
    for fn in sorted(os.listdir(os.path.join(B, 'fq'))):
        p = fn.split('.')[0]
        with gzip.open(os.path.join(B, 'fq', fn), 'rt') as fh:
            ln = fh.read().split('\n')
        for i in range(0, len(ln) - 3, 4): seqs[(p, ln[i][1:])] = (ln[i + 1], ln[i + 3])
    # ---- one CHM13 mapping of all conclusive reads
    od = os.path.join(B, 'map'); os.makedirs(od, exist_ok=True)
    fq = os.path.join(od, 'all_conclusive.fq'); sam = os.path.join(od, 'all_conclusive.chm13.sam')
    keys = sorted({(p, r) for p, r in zip(conc.individual, conc.read_id)})
    if not os.path.exists(sam):
        with open(fq, 'w') as fh:
            for p, r in keys:
                s, q = seqs[(p, r)]; fh.write(f'@{p}~{r}\n{s}\n+\n{q}\n')
        log('mapping to CHM13 v2.0'); _mm2(a.genome, fq, sam + '.part', a.threads); os.replace(sam + '.part', sam)
    recs = collections.defaultdict(list)
    with pysam.AlignmentFile(sam) as fi:
        for r in fi:
            if r.is_secondary: continue
            for t in ('MM', 'ML', 'MN'):
                if r.has_tag(t): r.set_tag(t, None)
            recs[r.query_name].append(r.to_string())
    sq = [l for l in open(os.path.join(B, 'chm13v2.0.dict')) if l.startswith('@SQ')]
    sq = [dict(x.split(':', 1) for x in l.rstrip('\n').split('\t')[1:]) for l in sq]
    for d in sq: d['LN'] = int(d['LN'])
    outroot = os.path.join(B, 'claims'); man = []
    for cid, g in conc.groupby('claim_id', sort=True):
        c = C.set_index('claim_id').loc[cid]; slug = c['slug']; cdir = os.path.join(outroot, c['cohort'], slug)
        os.makedirs(cdir, exist_ok=True)
        rgs = {}
        for r in g.itertuples():
            rid = r.dataset if isinstance(r.dataset, str) else r.individual
            rgs.setdefault(rid, dict(ID=rid, SM=r.individual, PL=DS_PL(r.platform), PM=str(r.platform),
                                     DS=f'source={r.source_url}; index={r.index_url}; region={r.region}; '
                                        f'source_reference={r.source_reference}'))
        def tagged(reflabel_rows, name_map):
            out = []
            for r in g.itertuples():
                k = f'{r.individual}~{r.read_id}'
                for s in name_map.get(k, []):
                    f = s.split('\t'); f[0] = r.read_id
                    rid = r.dataset if isinstance(r.dataset, str) else r.individual
                    f = [x for x in f if not x.startswith('RG:Z:')]
                    lo, hi = L.f(r.switch_ci_read_lo), L.f(r.switch_ci_read_hi)
                    f += [f'RG:Z:{rid}', f'XC:Z:{cid}', f'XT:Z:{r.tier}', f'XV:Z:{r.v1_verdict}',
                          f'XW:Z:{int(lo) if lo is not None else ""}-{int(hi) if hi is not None else ""}',
                          f'XU:i:{int(L.f(r.anchor_unique_bp) or 0)}', f'XO:Z:{r.own_best_end if isinstance(r.own_best_end, str) else ""}',
                          f'XN:i:{int(L.b(r.counted))}']
                    out.append('\t'.join(f))
            return out
        # dedupe identical read ids (same molecule in two datasets keeps the counted row)
        g = g.sort_values('counted', ascending=False).drop_duplicates('read_id')
        hdr = _header(sq, list(rgs.values()), 'reads aligned to CHM13 v2.0 (minimap2 map-ont)', cid, 'chm13v2.0.fa.gz')
        bam = os.path.join(cdir, f'{slug}.chm13.bam')
        _write_sorted(hdr, tagged(None, recs), bam); _to_cram(bam, os.path.join(cdir, f'{slug}.chm13.cram'), a.chm13)
        # own-assembly end windows: claimed end(s) + every other end of the person at the same arm
        ce = c['claimed_end_ids'].split(',')
        e0 = Ei.loc[ce[0]]
        same = E[(E['sample'] == e0['sample']) & (E['set'] == e0['set']) & (E.home_arm == e0['home_arm'])].end_id.tolist()
        wins = list(dict.fromkeys(ce + same))
        wf = pysam.FastaFile(os.path.join(REPO, 'work', 'windows', 'ends.fa'))
        wfa = os.path.join(cdir, f'{slug}.own_windows.fa'); wsq = []
        with open(wfa, 'w') as fh:
            for e in wins:
                s = wf.fetch(e); x = Ei.loc[e]
                nm = e.replace('#', '_')
                fh.write(f">{nm} end_id={e} claimed={'yes' if e in ce else 'no'} allele={x['structural_allele_id']} tip_frame=1 "
                         f"source={L.https(x['assembly_url'])} region={x['contig']}:{int(x['win_start'])}-{int(x['win_end'])} "
                         f"strand={x['strand']} sha256={x['window_sha256']}\n")
                for i in range(0, len(s), 80): fh.write(s[i:i + 80] + '\n')
                wsq.append(dict(SN=nm, LN=len(s), M5=_md5(s), AS=x['assembly_name'] if isinstance(x['assembly_name'], str) else '',
                                DS=f"end_id={e}; claimed={'yes' if e in ce else 'no'}; tip frame (position 1 = telomeric tip)"))
        pysam.faidx(wfa)
        tfq = os.path.join(od, f'{slug}.fq')
        with open(tfq, 'w') as fh:
            for r in g.itertuples():
                s, q = seqs[(r.individual, r.read_id)]; fh.write(f'@{r.individual}~{r.read_id}\n{s}\n+\n{q}\n')
        osam = os.path.join(od, f'{slug}.own.sam'); _mm2(wfa, tfq, osam, 2, soft_supp=False)
        orec = collections.defaultdict(list)
        with pysam.AlignmentFile(osam) as fi:
            for r in fi:
                if r.is_secondary: continue
                for t in ('MM', 'ML', 'MN'):
                    if r.has_tag(t): r.set_tag(t, None)
                orec[r.query_name].append(r.to_string())
        hdr = _header(wsq, list(rgs.values()), "reads aligned to the person's own assembly end windows (tip frame)", cid,
                      f'{slug}.own_windows.fa', y=False)
        _write_sorted(hdr, tagged(None, orec), os.path.join(cdir, f'{slug}.own.bam'))
        os.remove(tfq); os.remove(osam)
        g.to_csv(os.path.join(cdir, f'{slug}.reads.tsv'), sep='\t', index=False)
        man.append(dict(claim_id=cid, cohort=c['cohort'], slug=slug, n_reads=len(g), n_definitive=int((g.tier == 'definitive').sum()),
                        definitive_verdict=c['definitive_verdict'], dir=os.path.relpath(cdir, B)))
    pd.DataFrame(man).to_csv(os.path.join(B, 'claims_manifest.tsv'), sep='\t', index=False)
    log(f'{len(man)} claim directories')
    showcase(a); release(a)


def showcase(a):
    """repository subset: up to SHOWCASE_MAX_READS definitive reads per showcase claim (most diagnostic units first)"""
    import pandas as pd, pysam
    sd = os.path.join(REPO, 'evidence', 'showcase'); os.makedirs(sd, exist_ok=True)
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False).set_index('claim_id')
    rows = []
    for cid in SHOWCASE:
        c = C.loc[cid]; src = os.path.join(B, 'claims', c['cohort'], c['slug']); dst = os.path.join(sd, c['slug'])
        os.makedirs(dst, exist_ok=True)
        g = pd.read_csv(os.path.join(src, f"{c['slug']}.reads.tsv"), sep='\t')
        g = g[(g.tier == 'definitive') & (g.counted == True)]  # noqa: E712
        g = g.assign(_u=g.seg_D_units.fillna(0) + g.seg_P_units.fillna(0)).sort_values(['_u', 'anchor_unique_bp'], ascending=False)
        pick = g.drop_duplicates('dataset').head(SHOWCASE_MAX_READS)       # one per dataset first (independence), then fill
        pick = pd.concat([pick, g[~g.read_id.isin(pick.read_id)]]).head(SHOWCASE_MAX_READS).drop(columns='_u')
        names = set(pick.read_id)
        for kind in ('chm13', 'own'):
            with pysam.AlignmentFile(os.path.join(src, f"{c['slug']}.{kind}.bam")) as fi:
                out = os.path.join(dst, f"{c['slug']}.{kind}.bam")
                with pysam.AlignmentFile(out, 'wb', header=fi.header) as fo:
                    for r in fi:
                        if r.query_name in names: fo.write(r)
            pysam.index(out)
        _to_cram(os.path.join(dst, f"{c['slug']}.chm13.bam"), os.path.join(dst, f"{c['slug']}.chm13.cram"), a.chm13)
        shutil.copy(os.path.join(src, f"{c['slug']}.own_windows.fa"), dst); pysam.faidx(os.path.join(dst, f"{c['slug']}.own_windows.fa"))
        pick.to_csv(os.path.join(dst, f"{c['slug']}.reads.tsv"), sep='\t', index=False)
        rows.append(dict(claim_id=cid, slug=c['slug'], cohort=c['cohort'], n_definitive_reads_in_claim=int(c['n_definitive']),
                         n_reads_in_showcase=len(pick), datasets=','.join(sorted(pick.dataset.astype(str).unique()))))
    pd.DataFrame(rows).to_csv(os.path.join(sd, 'showcase_claims.tsv'), sep='\t', index=False)
    _sums(sd, os.path.join(sd, 'SHA256SUMS'))
    log(f'showcase: {sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(sd) for f in fs) / 1e6:.1f} MB')


def _sums(root, out, exclude=('SHA256SUMS',)):
    lines = []
    for dp, _, fs in os.walk(root):
        for f in sorted(fs):
            if f in exclude: continue
            p = os.path.join(dp, f); lines.append(f'{L.sha256(p)}  {os.path.relpath(p, root)}')
    with open(out, 'w') as fh: fh.write('\n'.join(sorted(lines, key=lambda x: x.split('  ')[1])) + '\n')


TABLES_CONTENT = ('claims.tsv, reads.tsv.gz, controls.tsv, criteria, population_reads_all.tsv.gz '
                  '(all 77,333 evaluated 1KG-ONT read rows)')


def release_tables(rd, filt):
    """tables bundle. evidence/population/reads.tsv.gz (every evaluated population read row; written by
    population_report.py, git-ignored because of its size) ships here as read_evidence_v1/population_reads_all.tsv.gz."""
    import tarfile
    out = os.path.join(rd, 'subtel_read_evidence_v1_tables.tar.gz')
    with tarfile.open(out, 'w:gz') as tf:
        for f in ('claims.tsv', 'reads.tsv.gz', 'controls.tsv', 'criteria_v1.json'):
            tf.add(os.path.join(REPO, 'evidence', f), arcname=f'read_evidence_v1/{f}', filter=filt)
        tf.add(os.path.join(B, 'claims_manifest.tsv'), arcname='read_evidence_v1/claims_manifest.tsv', filter=filt)
        tf.add(os.path.join(REPO, 'evidence', 'population', 'reads.tsv.gz'),
               arcname='read_evidence_v1/population_reads_all.tsv.gz', filter=filt)
    return dict(file=os.path.basename(out), bytes=os.path.getsize(out), content=TABLES_CONTENT)


def release(a):
    """work/evidence/release/: one tar.gz per cohort (claim directories) + tables, SHA256SUMS, MANIFEST.tsv"""
    import tarfile, pandas as pd
    rd = os.path.join(W, 'release'); os.makedirs(rd, exist_ok=True)
    def filt(ti):
        ti.uid = ti.gid = 0; ti.uname = ti.gname = ''; ti.mtime = 0
        return ti
    man = []
    for cohort in sorted(os.listdir(os.path.join(B, 'claims'))):
        out = os.path.join(rd, f'subtel_read_evidence_v1_{cohort}.tar.gz')
        with tarfile.open(out, 'w:gz') as tf:
            tf.add(os.path.join(B, 'claims', cohort), arcname=f'read_evidence_v1/{cohort}', filter=filt)
        man.append(dict(file=os.path.basename(out), bytes=os.path.getsize(out), content=f'per-claim CRAM/BAM/window FASTA, {cohort}'))
    man.append(release_tables(rd, filt))
    for m in man: m['sha256'] = L.sha256(os.path.join(rd, m['file']))
    M = pd.DataFrame(man)[['file', 'bytes', 'sha256', 'content']]
    M.to_csv(os.path.join(rd, 'MANIFEST.tsv'), sep='\t', index=False)
    _sums(rd, os.path.join(rd, 'SHA256SUMS'), exclude=('SHA256SUMS', 'MANIFEST.tsv'))
    M.to_csv(os.path.join(REPO, 'evidence', 'release_manifest.tsv'), sep='\t', index=False)   # committed copy (GitHub release assets)
    log('release: ' + ', '.join(f"{m['file']} {m['bytes'] / 1e6:.1f} MB" for m in man))


# ------------------------------------------------------------------------------------------------------------------ segments
def segments(a):
    """Showcase only: each read split at its switch interval, in tip frame - distal segment (telomere boundary to the lower
    switch bound) and proximal segment (upper bound to the read end) - aligned genome-wide to CHM13 v2.0 (minimap2 map-ont).
    Records are named <read_id>:distal / <read_id>:proximal (<slug>.segments.chm13.bam, segments.tsv)."""
    import pandas as pd, pysam
    from subtel_atlas import readevidence as RE
    th = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json'))['thresholds']
    S = pd.read_csv(os.path.join(REPO, 'evidence', 'showcase', 'showcase_claims.tsv'), sep='\t')
    seqs = {}
    for fn in sorted(os.listdir(os.path.join(B, 'fq'))):
        p = fn.split('.')[0]
        with gzip.open(os.path.join(B, 'fq', fn), 'rt') as fh: ln = fh.read().split('\n')
        for i in range(0, len(ln) - 3, 4): seqs[(p, ln[i][1:])] = (ln[i + 1], ln[i + 3])
    od = os.path.join(B, 'segments'); os.makedirs(od, exist_ok=True); fq = os.path.join(od, 'showcase_segments.fq')
    meta = {}
    with open(fq, 'w') as fh:
        for s in S.itertuples():
            rr = pd.read_csv(os.path.join(REPO, 'evidence', 'showcase', s.slug, f'{s.slug}.reads.tsv'), sep='\t')
            for r in rr.itertuples():
                seq, q = seqs[(r.individual, r.read_id)]
                t, o, tb = RE.orient(seq, r.orient_side, th, RE.telomere_motifs(th, r.platform))
                qq = q if o == '+' else q[::-1]
                assert len(t) == int(r.read_len), (r.read_id, len(t), r.read_len)
                s0 = int(r.s0_read); lo, hi = int(r.switch_ci_read_lo), int(r.switch_ci_read_hi)
                for part, (x0, x1) in (('distal', (s0, lo)), ('proximal', (hi, len(t)))):
                    nm = f'{s.slug}~{r.read_id}:{part}'; fh.write(f'@{nm}\n{t[x0:x1]}\n+\n{qq[x0:x1]}\n')
                    meta[nm] = (s.slug, r.read_id, part, x0, x1, r.dataset if isinstance(r.dataset, str) else r.individual)
    sam = os.path.join(od, 'showcase_segments.sam'); _mm2(a.genome, fq, sam, a.threads)
    sq = [l for l in open(os.path.join(B, 'chm13v2.0.dict')) if l.startswith('@SQ')]
    sq = [dict(x.split(':', 1) for x in l.rstrip('\n').split('\t')[1:]) for l in sq]
    for d in sq: d['LN'] = int(d['LN'])
    by = collections.defaultdict(list); rows = []
    with pysam.AlignmentFile(sam) as fi:
        for r in fi:
            if r.is_secondary: continue
            slug, rid, part, x0, x1, ds = meta[r.query_name]
            f = r.to_string().split('\t'); f[0] = f'{rid}:{part}'
            f += [f'RG:Z:{ds}', f'XP:Z:{part}', f'XR:Z:{x0}-{x1}']
            by[slug].append('\t'.join(f))
            if not r.is_supplementary:
                rows.append(dict(slug=slug, read_id=rid, segment=part, tip_frame_read_interval=f'{x0}-{x1}', length=x1 - x0,
                                 chrom=r.reference_name if not r.is_unmapped else '', pos=(r.reference_start + 1) if not r.is_unmapped else '',
                                 end=r.reference_end if not r.is_unmapped else '', mapq=r.mapping_quality if not r.is_unmapped else '',
                                 strand='-' if r.is_reverse else '+', has_supplementary=r.has_tag('SA')))
    T = pd.DataFrame(rows)
    for s in S.itertuples():
        d = os.path.join(REPO, 'evidence', 'showcase', s.slug)
        with pysam.AlignmentFile(os.path.join(d, f'{s.slug}.chm13.bam')) as fi: rg = fi.header.to_dict().get('RG', [])
        hdr = _header(sq, rg, 'read segments split at the read switch interval (distal: telomere boundary to switch; proximal: '
                      'switch to read end), tip frame, aligned to CHM13 v2.0', s.claim_id, 'chm13v2.0.fa.gz')
        out = os.path.join(d, f'{s.slug}.segments.chm13.bam'); _write_sorted(hdr, by[s.slug], out)
        T[T.slug == s.slug].drop(columns='slug').to_csv(os.path.join(d, f'{s.slug}.segments.tsv'), sep='\t', index=False)
    log(f'segment alignments for {T.read_id.nunique()} showcase reads')


# ------------------------------------------------------------------------------------------------------------------ igv
IGV_GENOME = L.CHM13_PLAIN_URL     # public uncompressed CHM13 v2.0 FASTA with .fai (IGV also offers the hosted genome 'hs1')


def _intervals(mask, offset):
    """runs of True in a window-relative boolean mask -> [(tip_lo, tip_hi)]"""
    import numpy as np
    m = np.asarray(mask, bool); out = []
    if not m.any(): return out
    d = np.diff(np.r_[0, m.astype(int), 0]); st = np.where(d == 1)[0]; en = np.where(d == -1)[0]
    return [(int(a) + offset, int(b) + offset) for a, b in zip(st, en) if b - a >= 200]


def _bed_model(ctx, fai):
    """BED9 of the claim model on CHM13: SNV diagnostic sites (on D and on P), D-private and P-private sequence, homologous blocks"""
    out = []; D, P = ctx.c['D'], ctx.c['P']
    for (dp, pp), (da, pa) in zip(ctx.pm.sites.tolist(), ctx.pm.site_alleles):
        for end, x, lab, col in ((D, dp, 'donor_site', '213,94,0'), (P, pp, 'proximal_site', '0,114,178')):
            pos = L.chm13_coord(fai, end, x)
            out.append((end[:-1], pos - 1, pos, f'{lab}:D={da}/P={pa}', 0, '.', pos - 1, pos, col))
    for end, mask, off, lab, col in ((D, ctx.pm.d_private, ctx.d_win[0], 'donor_private', '213,94,0'),
                                     (P, ctx.pm.p_private, ctx.p_win[0], 'proximal_private', '0,114,178')):
        for a, b in _intervals(mask, off):
            c, s, e = L.chm13_interval(fai, end, a, b - 1)
            out.append((c, s - 1, e, lab, 0, '.', s - 1, e, col))
    for bl in ctx.pm.blocks:
        for end, a, b in ((D, bl['d_start'], bl['d_end']), (P, bl['p_start'], bl['p_end'])):
            c, s, e = L.chm13_interval(fai, end, a, b - 1)
            out.append((c, s - 1, e, f"homologous_block:{bl['ident']:.3f}", 0, '.', s - 1, s - 1, '150,150,150'))
    return sorted(out, key=lambda x: (x[0], x[1]))


def _locus(fai, end, lo, hi, pad=2000):
    c, s, e = L.chm13_interval(fai, end, max(0, lo - pad), hi + pad)
    return f'{c}:{max(1, s)}-{min(e, fai[c])}'


def igv(a):
    import pandas as pd
    from subtel_atlas import readevidence as RE
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False).set_index('claim_id')
    S = pd.read_csv(os.path.join(REPO, 'evidence', 'showcase', 'showcase_claims.tsv'), sep='\t')
    crit = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json')); fai = _fai(a)
    ends = RE.Chm13Ends(a.chm13_ends); cat = os.path.join(REPO, 'data', 'catalogue')
    mask_cache = os.path.join(REPO, 'work', 'readevidence', 'masks')
    idx = []
    for s in S.itertuples():
        c = C.loc[s.claim_id]; d = os.path.join(REPO, 'evidence', 'showcase', s.slug)
        model = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=[c['model_end']])[0]
        ctx = RE.ClaimContext(model, ends, crit, preset='map-ont', mask_cache=mask_cache)
        D, P, R = model['D'], model['P'], model['R']
        bed = _bed_model(ctx, fai)
        bfn = f'{s.slug}.diagnostic_sites.chm13.bed'
        with open(os.path.join(d, bfn), 'w') as fh:
            fh.write(f'track name="claim model {s.claim_id}" itemRgb="On" description="CHM13 {D} (orange) vs {P} (blue): diagnostic '
                     f'SNV sites, private sequence, homologous blocks (grey), criteria v1"\n')
            for b in bed: fh.write('\t'.join(str(x) for x in b) + '\n')
        rr = pd.read_csv(os.path.join(d, f'{s.slug}.reads.tsv'), sep='\t')
        sg = pd.read_csv(os.path.join(d, f'{s.slug}.segments.tsv'), sep='\t')
        # loci: donor block on D (tip-frame span of the D units of the showcase reads, from the window) and the P/R end from the
        # telomere to the furthest anchor end
        an_hi = max([int(L.f(x) or 0) for x in rr.anchor_len] + [0])
        p_hi = max(60000, int(max(L.f(x) or 0 for x in rr.switch_ci_P_tip_hi)) + an_hi)
        locP = _locus(fai, R, 0, p_hi, 0)
        locD = _locus(fai, D, ctx.d_win[0], ctx.d_win[1], 0)
        loci = f'{locD} {locP}' if D[:-1] != R[:-1] or True else locP
        sess = os.path.join(d, f'{s.slug}.chm13.igv_session.xml')
        with open(sess, 'w') as fh:
            fh.write(f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{IGV_GENOME}" locus="{loci}" version="8">
    <Resources>
        <Resource name="read segments split at the switch (distal / proximal)" path="{s.slug}.segments.chm13.bam" index="{s.slug}.segments.chm13.bam.bai" type="bam"/>
        <Resource name="whole reads on CHM13 v2.0" path="{s.slug}.chm13.bam" index="{s.slug}.chm13.bam.bai" type="bam"/>
        <Resource name="claim model ({D} vs {P})" path="{bfn}" type="bed"/>
    </Resources>
</Session>
''')
        w0 = c['claimed_end_ids'].split(',')[0].replace('#', '_')
        with open(os.path.join(d, f'{s.slug}.own.igv_session.xml'), 'w') as fh:
            fh.write(f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{s.slug}.own_windows.fa" locus="{w0}:1-100000" version="8">
    <Resources>
        <Resource name="reads on the person's own assembly end windows" path="{s.slug}.own.bam" index="{s.slug}.own.bam.bai" type="bam"/>
    </Resources>
</Session>
''')
        # data-driven description
        dist = sg[sg.segment == 'distal']; prox = sg[sg.segment == 'proximal']
        n = len(rr); nD = int((dist.chrom == D[:-1]).sum()); nP = int((prox.chrom == R[:-1]).sum())
        nP60 = int(((prox.chrom == R[:-1]) & (prox.mapq.astype(float) >= 60)).sum())
        dq = sorted(dist[dist.chrom == D[:-1]].mapq.astype(int).tolist())
        oth = ', '.join(f'{k} {v}' for k, v in dist[dist.chrom != D[:-1]].chrom.value_counts().items())
        n_sites = int(len(ctx.pm.sites)); n_blk = len(ctx.pm.blocks)
        own_ok = int(rr.own_on_claimed_end.astype(str).eq('True').sum())
        sw = rr.switch_ci_chm13.dropna().tolist()
        units_txt = (f'{n_sites} SNV sites in {n_blk} homologous block(s) plus private 200 bp bins' if n_sites else
                     f'no SNV sites (no homologous block of >= 2 kb at >= 90% between the two windows), so all units are private '
                     f'200 bp bins: sequence present in CHM13 {D[3:]} but absent from CHM13 {P[3:]}, and the reverse')
        txt = f"""# {s.claim_id}

**Claim.** {c['individual']} {c['arm']} ({c['haplotype']}), structural allele {c['allele_id']} ({c['nativity_class'].replace('_', ' ')}):
the telomere-adjacent block is closer to CHM13 {D[3:]} (donor D) than to CHM13 {R[3:]}; the assembly places the switch back
to CHM13 {P[3:]} {int(c['assembly_switch_sub']):,} bp from the telomere boundary (CHM13 {c['assembly_switch_chm13_P'] or 'projection undefined'}).

**Verdicts.** v1 (pre-registered): **{c['v1_verdict']}**, {int(c['n_conclusive_concordant'])} conclusive switch-concordant reads.
Definitive tier (post hoc, stricter): **{c['definitive_verdict'].replace('_', ' ')}**, {int(c['n_definitive'])} definitive reads
(datasets: {c['definitive_datasets']}). This directory holds {n} of them; the release bundle has all reads of the claim.

**Unit model** (CHM13 only): {units_txt}.

## Files

| file | content |
|---|---|
| `{s.slug}.reads.tsv` | per-read evidence (units, purity, switch interval, anchor, own-assembly placement, public source, samtools command) |
| `{s.slug}.chm13.cram`, `.chm13.bam` | whole reads aligned to CHM13 v2.0 (minimap2 map-ont, -Y) |
| `{s.slug}.segments.chm13.bam`, `.segments.tsv` | each read split at its switch interval; distal and proximal segments aligned separately |
| `{s.slug}.own.bam`, `{s.slug}.own_windows.fa` | reads aligned to the person's own assembly end windows (tip frame: position 1 = telomeric tip) |
| `{bfn}` | the claim model on CHM13: diagnostic SNV sites, donor- and proximal-private sequence, homologous blocks |
| `{s.slug}.chm13.igv_session.xml` | IGV: two loci, donor block {locD} and recipient end {locP} |
| `{s.slug}.own.igv_session.xml` | IGV on the own-assembly windows |

## What you should see

1. **CHM13 session, track "read segments".** The distal segment of each read (telomere boundary to the switch) aligns to
   the donor end CHM13 {D[3:]} for {nD} of {n} reads{f' (MAPQ {dq[0]}-{dq[-1]})' if dq else ''}{f'; elsewhere: {oth}' if oth else ''}.
   The proximal segment aligns to CHM13 {R[3:]} for {nP} of {n} reads, {nP60} of them at MAPQ 60. The reads therefore join
   {D[3:]}-like sequence at the telomere to {R[3:]} sequence, switching inside CHM13 {sw[0] if sw else 'the interval in reads.tsv'}
   (switch interval of the first read; all reads in `reads.tsv`).
2. **Track "claim model".** Orange features are {D[3:]} sequence and alleles, blue features {P[3:]}. With mismatch colouring
   on, the distal segments match the orange features and the proximal segments the blue ones.
3. **Track "whole reads".** The whole-read alignment of a carrier read often stays on one chromosome. Donor and recipient
   copies are near-identical in shared blocks, so the aligner need not split the read. This is why the evidence uses
   diagnostic units, not split alignments.
4. **Own-assembly session.** {own_ok} of {n} reads place best on the claimed end window `{w0}` (`claimed=yes` in the FASTA
   header). The other windows are the person's other ends at this arm.
"""
        open(os.path.join(d, 'README.md'), 'w').write(txt)
        idx.append(dict(claim_id=s.claim_id, slug=s.slug, D=D, P=P, R=R, locus_donor=locD, locus_recipient=locP, n_reads=n,
                        distal_on_D=nD, proximal_on_R=nP, proximal_on_R_mapq60=nP60, own_on_claimed=own_ok, snv_sites=n_sites))
    pd.DataFrame(idx).to_csv(os.path.join(REPO, 'evidence', 'showcase', 'showcase_loci.tsv'), sep='\t', index=False)
    _sums(os.path.join(REPO, 'evidence', 'showcase'), os.path.join(REPO, 'evidence', 'showcase', 'SHA256SUMS'))
    log(f'IGV sessions for {len(idx)} showcase claims')
    print(pd.DataFrame(idx).to_string())


def main():
    inp = os.environ.get('ATLAS_INPUT_DIR', '')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('stage', choices=['collect', 'tables', 'bams', 'segments', 'igv', 'fixcram', 'release', 'all'])
    ap.add_argument('--chm13', default=os.path.join(inp, 'ref', 'chm13v2.0.fa.gz'))
    ap.add_argument('--genome', default=os.path.join(REPO, 'work', 'readevidence', 'chm13v2.map-ont.mmi'))
    ap.add_argument('--chm13-ends', default=os.path.join(inp, 'subtel', 'foreign', 'ref', 'chm13_ends500k.fa'))
    ap.add_argument('-t', '--threads', type=int, default=6)
    a = ap.parse_args()
    for s in (['collect', 'tables', 'bams', 'segments', 'igv'] if a.stage == 'all' else [a.stage]):
        globals()[s](a)


if __name__ == '__main__':
    main()
