"""Claim-level verdicts, controls, per-allele aggregates and evidence BAMs for run_population.py (steps `report` and `bams`).

Claim rule: subtel_atlas.readevidence.claim_verdict, unphased (1000 Genomes ONT reads carry no haplotype assignment). An
'assessability' reason is attached to every claim that is not confirmed; reasons are tested in this order (first match wins):
  model_window_undefined  the donor span of the claimed block is not collinear on CHM13 D (d_lo > d_hi), so the D window of
                     readevidence.claim_windows is empty; not assessable under v1
  discordant         conclusive reads exist but none is switch-concordant with the assembly switch
  one_conclusive     exactly one conclusive concordant read (verdict supported or not_confirmed)
                     (these two come first: a conclusive read makes the claim assessable; with C2 a read can be conclusive and
                     catalogue-concordant even when the assembly projection fails)
  diagnostic_units   the claim's unit model applied to the person's own assembly end window (an error-free read) does not pass
                     the switch test, and no read is conclusive
  no_reads           no read >= 5 kb in the claim read set
  no_telomere_read   no read with a terminal telomere (section 3a)
  read_length        no telomere-reaching read long enough to span the switch plus a 10 kb anchor
                     (read subtelomere length < assembly switch + 10 kb)
  anchor             long enough telomere-reaching reads exist but none anchors at MAPQ >= 60 on the recipient chromosome
  integrity          anchored, switch-passing telomere-reaching reads fail integrity (internal telomere, adapter, foldback or
                     aligned fraction)
  contradicting_only anchored telomere-reaching reads exist; some show the native arm (contradicts) and none is conclusive
  switch_test        anchored, telomere-reaching, long-enough reads fail the switch test without contradicting
Not assessable = model_window_undefined, diagnostic_units, no_reads, no_telomere_read, read_length, anchor; the others are
'assessed, not confirmed' (verdict_class).
"""
import os, sys, json, csv, gzip, math, collections, subprocess, tempfile, shutil
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import readevidence as RE  # noqa: E402

READ_OUT = ['claim_id', 'role', 'sample', 'read', 'cram_url', 'slice_region', 'slice_end', 'R', 'D', 'P', 'claim_pattern', 'read_group', 'input_contig', 'input_pos',
            'input_strand', 'input_mapq', 'input_has_sa', 'read_len', 'orientation', 'terminal', 'telo_offset', 'telo_len', 's0_read',
            'n_sites', 'snv_D', 'snv_P', 'snv_other', 'snv_uncovered', 'bin_D', 'bin_P', 'n_units', 'segments', 'seg_D_units',
            'seg_D_purity', 'seg_P_units', 'seg_P_purity', 'dR_units', 'dR_purity', 'log10_lr_vs_allP', 'log10_lr_vs_allD',
            'switch_pass', 'switch_ci_sub_lo', 'switch_ci_sub_hi', 'switch1_ci_sub_lo', 'switch1_ci_sub_hi', 'switch_ci_P_lo',
            'switch_ci_P_hi', 'assembly_switch_sub', 'switch_concordant', 'concordance_ref', 'switch_concordant_catalogue',
            'switch_concordant_hap2', 'donor_interval_units', 'donor_interval_P_purity', 'anchor_pass', 'anchor_chrom', 'anchor_mapq',
            'anchor_len', 'anchor_tip_start', 'anchor_tip_end', 'anchor_unique_bp', 'anchor_tier', 'internal_telomere_bp',
            'adapter_hits', 'foldback', 'subtel_aligned_frac', 'integrity_pass', 'verdict', 'genome_alignments']


def _b(x):
    return x in (True, 'True', 'true', '1')


def _f(x):
    try:
        v = float(x); return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def load_rows(work):
    rows = list(csv.DictReader(open(os.path.join(work, 'reads_full.tsv')), delimiter='\t'))
    for r in rows:
        for k in ('switch_concordant', 'switch_pass', 'terminal', 'anchor_pass', 'integrity_pass'): r[k] = _b(r.get(k))
    return rows


def binom_upper(k, n, conf=0.95):
    """exact (Clopper-Pearson) one-sided upper bound"""
    from scipy.stats import beta
    return 1.0 if n == 0 else (1.0 if k >= n else float(beta.ppf(conf, k + 1, n - k)))


def reason(c, rr, asm, th):
    if not c.get('model_valid'): return 'model_window_undefined'
    conc = [r for r in rr if r['verdict'] == 'conclusive']      # conclusive reads make the claim assessable
    nconc = sum(r['switch_concordant'] for r in conc)
    if conc and not nconc: return 'discordant'
    if nconc == 1: return 'one_conclusive'
    if asm is None or not _b((asm or {}).get('switch_pass')): return 'diagnostic_units'
    if not rr: return 'no_reads'
    term = [r for r in rr if r['terminal']]
    if not term: return 'no_telomere_read'
    sw = (_f(asm.get('switch_ci_sub_hi')) or c['model'].get('switch_sub') or 0)
    need = sw + th['anchor_min_bp']
    longr = [r for r in term if int(r['read_len']) - int(r['s0_read']) >= need]
    if not longr: return 'read_length'
    anch = [r for r in longr if r['anchor_pass']]
    if not anch: return 'anchor'
    if any(r['switch_pass'] and r['anchor_pass'] and not r['integrity_pass'] for r in term): return 'integrity'
    if any(r['verdict'] == 'contradicts' for r in rr): return 'contradicting_only'
    return 'switch_test'


def unit_detail(pj, pattern, th):
    """why the assembly projection fails the switch test (first failing condition of section 3f)"""
    if not pj or not pj.get('segments'): return 'no_projection'
    if _b(pj.get('switch_pass')): return ''
    ss = [(x.split(':')[0], int(x.split(':')[1]), float(x.split(':')[2])) for x in pj['segments'].split(';')]
    if [x[0] for x in ss] != (['D', 'P'] if pattern == 'DP' else ['P', 'D', 'P']): return 'segment_order'
    for lab in ('D', 'P'):
        if any(l == lab and n < th['min_units_segment'] for l, n, _ in ss): return f'{lab}_units<{th["min_units_segment"]}'
    for lab in ('D', 'P'):
        if any(l == lab and p < th['min_segment_purity'] for l, _, p in ss): return f'{lab}_purity<{th["min_segment_purity"]}'
    if min(_f(pj.get('log10_lr_vs_allP')) or 0, _f(pj.get('log10_lr_vs_allD')) or 0) < th['min_log10_lr']: return 'log10_LR<5'
    return 'D_vs_R_check'


def report(a):
    import pandas as pd
    P = json.load(open(os.path.join(a.work, 'plan.json'))); crit = RE.load_criteria(a.criteria); th = crit['thresholds']
    rows = load_rows(a.work)
    asm = json.load(open(os.path.join(a.work, 'assembly_projection.json')))
    man = {r['file']: r for r in csv.DictReader(open(os.path.join(a.work, 'fetch_manifest.tsv')), delimiter='\t')}
    by = collections.defaultdict(list)
    for r in rows: by[(r['role'], r['claim_id'])].append(r)
    out = []
    for c in P['claims']:
        m = c['model'] or {}
        rall = by.get(('carrier', c['claim']), [])
        # claim read set: reads whose primary alignment is in the recipient-end slice, plus donor-slice reads only when conclusive
        # (anchored at MAPQ >= 60 on the recipient chromosome). A donor-slice read that is not anchored on R has no evidence of
        # recipient origin and is in nearly all cases a read of the donor arm itself; it must not count as 'consistent'.
        rr = [r for r in rall if r['slice_end'] == m.get('R') or r['verdict'] == 'conclusive']
        v = RE.claim_verdict(rr, crit, same_hap_key=None)
        v_all = RE.claim_verdict(rall, crit, same_hap_key=None)
        pj = asm.get(f"{c['claim']}||{c['model_end']}") if c['model'] else None
        conc = [r for r in rr if r['verdict'] == 'conclusive' and r['switch_concordant']]
        contra = [r for r in rr if r['verdict'] == 'contradicts']
        term = [r for r in rr if r['terminal']]
        frac = len(conc) / (len(conc) + len(contra)) if conc or contra else None
        ci = [(_f(r['switch_ci_sub_lo']), _f(r['switch_ci_sub_hi'])) for r in conc]
        sl = [x for x in P['slices'] if x['sample'] == c['sample'] and x['end'] in (m.get('R'), m.get('D'))]
        out.append(dict(
            claim_id=c['claim'], sample=c['sample'], arm=c['arm'], allele=c['allele'], tier=c['tier'],
            nativity_class=c['nativity_class'], n_haplotypes_with_allele=c['n_haplotypes'], end_ids=','.join(c['end_ids']),
            model_end=c['model_end'], other_haplotype=c['other_haplotype'], other_end=c['other_end'] or '',
            other_allele=c['other_allele'] or '', pattern=m.get('pattern'), D=m.get('D'), P=m.get('P'), R=m.get('R'),
            assembly_switch_sub=m.get('switch_sub'), assembly_p_switch=m.get('p_switch'),
            asm_projection_switch_pass=(pj or {}).get('switch_pass'), asm_projection_segments=(pj or {}).get('segments'),
            asm_projection_fail=unit_detail(pj, m.get('pattern'), th) if c.get('model_valid') else 'model_window_undefined', asm_projection_switch_ci_sub=
            f"{(pj or {}).get('switch_ci_sub_lo')}-{(pj or {}).get('switch_ci_sub_hi')}" if pj and pj.get('switch_pass') else '',
            cram_url=RE_url(c['sample']), regions=';'.join(x['region'] for x in sl),
            n_slice_records=sum(int(man.get(x['file'], {}).get('n_reads', 0)) for x in sl),
            n_reads=len(rr), n_recipient_slice_reads=sum(r['slice_end'] == m.get('R') for r in rall), n_telomere_reaching=len(term), n_telomere_anchored=sum(r['terminal'] and r['anchor_pass'] for r in rr),
            n_conclusive_concordant=v['n_conclusive_concordant'], n_conclusive_discordant=v['n_conclusive_discordant'],
            n_consistent=v['n_consistent'], n_contradicts=v['n_contradicts'],
            n_conclusive_via_donor_slice=sum(r['slice_end'] != m.get('R') for r in conc),
            conclusive_fraction=round(frac, 3) if frac is not None else None, expected_fraction=c['expected_fraction'],
            read_switch_ci_sub_median=f"{int(pd.Series([x[0] for x in ci]).median())}-{int(pd.Series([x[1] for x in ci]).median())}"
            if ci else '', verdict=v['verdict'], verdict_incl_donor_slice_consistent=v_all['verdict'],
            n_donor_slice_reads=len(rall) - sum(r['slice_end'] == m.get('R') for r in rall),
            n_donor_slice_consistent_excluded=sum(r['slice_end'] != m.get('R') and r['verdict'] == 'consistent' for r in rall),
            assessability='assessed_confirmed' if v['verdict'] == 'confirmed' else reason(c, rr, pj, th),
            criteria_sha=RE.criteria_hash(crit), platform='ONT_R9.4.1'))
    C = pd.DataFrame(out)
    # unphased reads cannot separate two non-native alleles of one individual at one arm when both models accept the same read
    cr = collections.defaultdict(set)
    for r in rows:
        if r['role'] == 'carrier' and r['verdict'] == 'conclusive' and r['switch_concordant']: cr[r['claim_id']].add(r['read'])
    sh, spec = [], []
    for c in C.itertuples():
        others = set().union(*[cr[o] for o in C[(C['sample'] == c.sample) & (C.arm == c.arm) & (C.claim_id != c.claim_id)].claim_id])
        mine = cr[c.claim_id]; sh.append(len(mine & others)); spec.append(len(mine - others))
    C['n_conclusive_shared_with_other_claim'] = sh; C['n_conclusive_allele_specific'] = spec
    C['allele_specific_confirmed'] = (C.verdict == 'confirmed') & (C.n_conclusive_allele_specific >= th['claim_min_conclusive'])
    # secondary own-assembly line (CRITERIA.md section 3, assembly_concordant): primary alignment of each conclusive read to the
    # individual's own Release 2 end windows at the arm (both haplotypes; from the `bams` step). Not used by the claim rule.
    own = own_placement(a)
    if own is not None:
        oc, oo, ok, oe = [], [], [], []
        for c in C.itertuples():
            mine = set(c.end_ids.split(',')); pl = [own.get((c.sample, n)) for n in cr[c.claim_id]]; pl = [p for p in pl if p]
            oc.append(sum(p['contig'] in mine for p in pl)); oo.append(sum(p['contig'] not in mine for p in pl))
            ok.append(sum(p['contig'] in mine and p['strict'] for p in pl))
            oe.append(','.join(sorted({p['contig'] for p in pl if p['contig'] not in mine})))
        C['n_conclusive_own_on_claimed_end'] = oc; C['n_conclusive_own_assembly_concordant'] = ok
        C['n_conclusive_own_on_other_end'] = oo; C['own_other_ends'] = oe
    NA = {'model_window_undefined', 'diagnostic_units', 'no_reads', 'no_telomere_read', 'read_length', 'anchor'}
    C['verdict_class'] = [v if v in ('confirmed', 'supported') else ('not_assessable' if r in NA else 'assessed_not_confirmed')
                          for v, r in zip(C.verdict, C.assessability)]
    os.makedirs(a.out, exist_ok=True)
    C.to_csv(os.path.join(a.out, 'claims.tsv'), sep='\t', index=False)
    # controls
    ctl = []
    for (role, cid), rr in sorted(by.items()):
        if role != 'control_native_strict': continue
        elig = [r for r in rr if r['terminal'] and r['anchor_pass'] and r['integrity_pass']]
        ctl.append(dict(control_id=cid, sample=cid.split('|')[0], arm=cid.split('|')[1], allele_model=rr[0]['model_claim'],
                        n_reads=len(rr), n_recipient_slice_reads=sum(r['slice_end'] == m.get('R') for r in rall), n_telomere_reaching=sum(r['terminal'] for r in rr), n_eligible=len(elig),
                        n_conclusive=sum(r['verdict'] == 'conclusive' for r in rr),
                        n_conclusive_concordant=sum(r['verdict'] == 'conclusive' and r['switch_concordant'] for r in rr),
                        n_consistent=sum(r['verdict'] == 'consistent' for r in rr), n_contradicts=sum(r['verdict'] == 'contradicts' for r in rr),
                        claim_rule_verdict=RE.claim_verdict(rr, crit)['verdict']))
    K = pd.DataFrame(ctl); K.to_csv(os.path.join(a.out, 'controls.tsv'), sep='\t', index=False)
    # alleles
    al = []
    for (arm, sa), g in C.groupby(['arm', 'allele']):
        assessable = g[~g.assessability.isin(['model_window_undefined', 'diagnostic_units', 'no_reads', 'no_telomere_read', 'read_length', 'anchor'])]
        kk = K[K.control_id.str.split('|').str[1:3].str.join('|') == f'{arm}|{sa}'] if len(K) else K
        al.append(dict(allele=sa, arm=arm, tier=g.tier.iloc[0], nativity_class=g.nativity_class.iloc[0], n_carriers_tested=len(g),
                       n_homozygous=int((g.n_haplotypes_with_allele == 2).sum()), n_confirmed=int((g.verdict == 'confirmed').sum()),
                       n_supported=int((g.verdict == 'supported').sum()), n_not_confirmed=int((g.verdict == 'not_confirmed').sum()),
                       n_allele_specific_confirmed=int(g.allele_specific_confirmed.sum()),
                       n_not_assessable=int((g.verdict_class == 'not_assessable').sum()),
                       n_assessed_not_confirmed=int((g.verdict_class == 'assessed_not_confirmed').sum()),
                       n_reached_anchor_stage=len(assessable),
                       not_confirmed_reasons=';'.join(f'{k}:{v}' for k, v in g[g.verdict != 'confirmed'].assessability.value_counts().items()),
                       n_conclusive_reads=int(g.n_conclusive_concordant.sum()), n_contradicting_reads=int(g.n_contradicts.sum()),
                       n_control_samples=int(kk['sample'].nunique()) if len(kk) else 0,
                       n_control_eligible_reads=int(kk.n_eligible.sum()) if len(kk) else 0,
                       n_control_conclusive=int(kk.n_conclusive.sum()) if len(kk) else 0))
    pd.DataFrame(al).sort_values(['arm', 'allele']).to_csv(os.path.join(a.out, 'alleles.tsv'), sep='\t', index=False)
    # reads: every carrier and control read row (committed copy, compact columns)
    with gzip.open(os.path.join(a.out, 'reads.tsv.gz'), 'wt') as f:
        f.write('\t'.join(READ_OUT) + '\n')
        for r in rows: f.write('\t'.join('' if r.get(k) in (None, 'None') else str(r.get(k)) for k in READ_OUT) + '\n')
    # reads.tsv: reads in the claim read sets (and control reads) that reach the telomere or are not non-informative
    cols = [k for k in READ_OUT if k != 'genome_alignments']
    inf = [r for r in rows if (r['verdict'] != 'non-informative' or r['terminal'])
           and (r['role'] != 'carrier' or r['slice_end'] == r['R'] or r['verdict'] == 'conclusive')]
    with open(os.path.join(a.out, 'reads.tsv'), 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in inf: f.write('\t'.join('' if r.get(k) in (None, 'None') else str(r.get(k)) for k in cols) + '\n')
    print(C.verdict.value_counts().to_dict(), C.assessability.value_counts().to_dict(), file=sys.stderr)


def own_placement(a):
    """{(sample, read): {contig, mapq, strict}} from <work>/bam/<sample>.own.bam; strict = primary on the end window, no
    supplementary alignment, at least 95% of the read aligned (the HG002 assembly_concordant definition)"""
    man = os.path.join(a.work, 'bam', 'MANIFEST.tsv')
    if not os.path.exists(man): return None
    import pysam
    out = {}
    for m in csv.DictReader(open(man), delimiter='\t'):
        if m['kind'] != 'own': continue
        with pysam.AlignmentFile(os.path.join(REPO, m['file'])) as f:
            for r in f:
                if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
                q = (r.query_alignment_end - r.query_alignment_start) / (r.infer_read_length() or 1)
                out[(m['sample'], r.query_name)] = dict(contig=r.reference_name, mapq=r.mapping_quality,
                                                        strict=bool(not r.has_tag('SA') and q >= 0.95))
    return out


def RE_url(s):
    return 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/{s}.t2t.cram'.format(s=s)


def bams(a):
    """Conclusive reads (carrier and control) re-aligned with minimap2 -ax map-ont to CHM13 v2.0 and to the individual's own
    Release 2 end windows at the claimed arm (both haplotypes; tip-frame windows, telomere at position 1). One read group per
    sample: ID/SM = sample, DS = source CRAM URL. Each record carries XC:Z (claim ids) and XS:Z (source slice region)."""
    import pysam
    P = json.load(open(os.path.join(a.work, 'plan.json'))); rows = load_rows(a.work)
    import pandas as pd
    E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False, usecols=['end_id', 'sample', 'set', 'home_arm'])
    E = E[E['set'] == 'HPRC_R2']
    conc = [r for r in rows if r['verdict'] == 'conclusive']
    bys = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in conc: bys[r['sample']][(r['read'], r['slice_region'], r['cram_url'])].add(r['claim_id'])
    od = os.path.join(a.work, 'bam'); os.makedirs(os.path.join(od, 'ref'), exist_ok=True)
    wf = pysam.FastaFile(a.windows); mm = RE.find_minimap2(); man = []
    slices = {x['region'] + '|' + x['sample']: x for x in P['slices']}
    for s, reads in sorted(bys.items()):
        want = collections.defaultdict(set)
        for (n, reg, url), cids in reads.items(): want[reg].add(n)
        td = tempfile.mkdtemp(prefix='popev_')   # scratch outside the repository
        fq = os.path.join(td, f'{s}.fq')
        with open(fq, 'w') as f:
            for reg, names in want.items():
                x = slices[reg + '|' + s]
                for name, seq, rec in RE.iter_reads(os.path.join(a.work, x['file'])):
                    if name in names:
                        q = ''.join(chr(min(v, 93) + 33) for v in rec.query_qualities) if rec.query_qualities is not None else 'I' * len(seq)
                        f.write(f'@{name}\n{seq}\n+\n{q}\n')
        arms = sorted({cid.split('|')[1] for v in reads.values() for cid in v})
        own = os.path.join(od, 'ref', f'{s}.own_windows.fa')
        with open(own, 'w') as f:
            for e in E[(E['sample'] == s) & E.home_arm.isin(arms)].end_id:
                if e in wf.references: f.write(f'>{e}\n{wf.fetch(e)}\n')
        url = RE_url(s); rg = f'@RG\\tID:{s}\\tSM:{s}\\tPL:ONT\\tDS:{url}'
        tag = {n: (','.join(sorted(c)), reg) for (n, reg, u), c in reads.items()}
        for ref, kind in ((a.genome, 'chm13'), (own, 'own')):
            out = os.path.join(od, f'{s}.{kind}.bam'); sam = os.path.join(td, f'{kind}.sam')
            with open(sam, 'w') as so:
                subprocess.run([mm, '-ax', 'map-ont', '-Y', '--MD', '-t', str(a.threads), '-R', rg, ref, fq], stdout=so,
                               stderr=subprocess.DEVNULL, check=True)
            with pysam.AlignmentFile(sam) as fi:
                hdr = fi.header.to_dict()
                # program lines without local paths (published evidence set)
                hdr['PG'] = [dict(ID='minimap2', PN='minimap2', VN=p.get('VN', ''), CL=f"minimap2 -ax map-ont -Y --MD -R <RG> "
                                  f"{'chm13v2.0' if kind == 'chm13' else s + '.own_windows.fa'} <conclusive_reads.fq>")
                             for p in hdr.get('PG', []) if p.get('ID') == 'minimap2'][:1]
                hdr.setdefault('CO', []).append(f'conclusive 1000 Genomes ONT reads (criteria v1 a89c5eaf023ed032), source {url}; '
                                                f'reference {"CHM13 v2.0" if kind == "chm13" else "own Release 2 end windows (tip frame)"}')
                tmp = os.path.join(td, f'{kind}.unsorted.bam')
                with pysam.AlignmentFile(tmp, 'wb', header=hdr) as fo:
                    for r in fi:
                        c, reg = tag.get(r.query_name, ('', ''))
                        r.set_tag('XC', c, 'Z'); r.set_tag('XS', reg, 'Z'); fo.write(r)
            pysam.sort('--no-PG', '-o', out, tmp); pysam.index(out)
            man.append(dict(sample=s, kind=kind, file=os.path.relpath(out, REPO), n_reads=len(reads), cram_url=url))
        shutil.rmtree(td, ignore_errors=True)
    with open(os.path.join(od, 'MANIFEST.tsv'), 'w') as f:
        w = csv.DictWriter(f, fieldnames=['sample', 'kind', 'file', 'n_reads', 'cram_url'], delimiter='\t'); w.writeheader(); w.writerows(man)
    print(f'{len(man)} evidence BAMs in {od}', file=sys.stderr)
