"""subtel-atlas command line: fetch-refs, query, reads, reads-evidence, verify-refs."""
import argparse, sys, os, time
from . import __version__, CATALOGUE_VERSION

DESC = f"""Assign human telomeric chromosome ends to the subtelomere allele catalogue {CATALOGUE_VERSION}
(arm, CHM13 mosaic labelled by paralog family, S5/S10/S20 terminal states, structural allele SA-<arm>-NN with its nativity
class, sequence haplogroup HG-<arm>-NN, boundary stability). Labels such as 'foreign' or a 'switch' describe which CHM13 v2.0
end (or paralog family of near-identical ends) a stretch of sequence most resembles; inside shared, paralogous subtelomeric
sequence they are not necessarily structural breakpoints."""

EPI = """examples:
  subtel-atlas fetch-refs                                  # CHM13 ends + catalogue tables into ~/.cache/subtel-atlas/v0.1
  subtel-atlas query assembly.fa.gz > ends.tsv             # all telomeric contig ends of an assembly
  subtel-atlas query --windows tips.fa --format json       # pre-cut telomere-first windows (header tokens arm=20q side=R)
  subtel-atlas reads telomeric_reads.fq.gz > reads.tsv     # experimental: telomeric long reads
  subtel-atlas reads-evidence --claim 'HG002#1#chr4_PATERNAL#R' --genome chm13v2.0.fa slice.bam -o 4q   # read-level evidence
"""


def _log(quiet):
    return (lambda *a: None) if quiet else (lambda *a: print(*a, file=sys.stderr, flush=True))


def cmd_fetch(a):
    from . import refs
    man = refs.fetch_refs(a.dir, chm13=a.chm13, catalogue=a.catalogue, chm13_url=a.chm13_url, raw_url=a.catalogue_url, force=a.force,
                          log=_log(a.quiet))
    print(f'reference directory ready: {a.dir or refs.default_dir()} (catalogue {man["catalogue_version"]})', file=sys.stderr)


def cmd_verify(a):
    from . import refs
    ok = refs.verify(a.dir or refs.default_dir())
    sys.exit(0 if all(ok.values()) else 1)


def _common_query(p):
    p.add_argument('--refs', metavar='DIR', help='reference directory from fetch-refs (default $SUBTEL_ATLAS_REFS or ~/.cache/subtel-atlas/v0.1)')
    p.add_argument('--no-verify-refs', action='store_true', help='skip the SHA-256 check of the reference directory (custom or test '
                   'references; results are then not versioned against the catalogue)')
    p.add_argument('--format', choices=['tsv', 'json'], default='tsv', help='output format (default tsv; json adds blocks, core blocks, '
                   'top-3 haplogroups, per-superpopulation frequencies and run metadata)')
    p.add_argument('-o', '--output', metavar='FILE', help='write to FILE instead of standard output')
    p.add_argument('--min-telomere', type=int, default=500, metavar='BP', help='minimum clean (CCCTAA)n block for a telomeric end (default 500, '
                   'the Release 2 tip definition)')
    p.add_argument('--aligner', choices=['auto', 'minimap2', 'mappy'], default='auto', help='aligner for the identity tracks and arm bins: '
                   'minimap2 executable (the catalogue computation, exact) or the mappy binding (approximation, see README); auto = '
                   'minimap2 if on PATH')
    p.add_argument('--minimap2', metavar='EXE', help='minimap2 executable (default: from PATH or $SUBTEL_ATLAS_MINIMAP2)')
    p.add_argument('-t', '--threads', type=int, default=1, help='threads/processes (default 1)')
    p.add_argument('--quiet', action='store_true', help='no progress messages on standard error')


def _arm_opts(p):
    p.add_argument('--no-name-arms', action='store_true', help='do not take the arm from chromosome-like contig names (chr1, chr1_MATERNAL, '
                   'HG002#1#chr1); every arm is then inferred')
    p.add_argument('--arm-map', metavar='TSV', help='two columns: end ID (<contig>#L|R or window name) and arm (e.g. 20q)')


def _run(a, ends, n_in, mode):
    from . import query as Q, bins as B
    log = _log(a.quiet)
    backend = a.aligner if a.aligner != 'auto' else ('minimap2' if B.find_minimap2(a.minimap2) else 'mappy')
    if backend == 'mappy':
        print('warning: mappy backend in use (minimap2 executable not found or not selected). It approximates the catalogue '
              'alignment: in validation ~92% of structural-allele calls agreed with the minimap2 computation. Install minimap2 '
              '(conda install -c bioconda minimap2) for catalogue-exact results.', file=sys.stderr)
    ref = Q.Reference(a.refs, verify=not a.no_verify_refs)
    genome = getattr(a, 'chm13_genome', None)
    if genome: log(f'  arm probes against {genome} ({"minimap2" if backend == "minimap2" else "mappy"})')
    t0 = time.time()
    recs = Q.analyse(ends, ref, backend, a.threads, genome, a.minimap2, log)
    log(f'  {len(recs)} ends analysed in {time.time() - t0:.1f} s')
    m = Q.meta(backend, dict(min_telomere=a.min_telomere, input=os.path.basename(a.input)), n_in, len(recs)); m['mode'] = mode
    fh = open(a.output, 'w') if a.output else sys.stdout
    (Q.write_json if a.format == 'json' else Q.write_tsv)(recs, fh, m)
    if a.output: fh.close()


def cmd_query(a):
    from . import query as Q
    log = _log(a.quiet); amap = Q.load_arm_map(a.arm_map); ends = []; n = 0
    for name, comment, seq in Q.iter_records(a.input):
        n += 1
        if a.windows: ends.append(Q.window_end(name, comment, seq, not a.no_name_arms, amap))
        else: ends += Q.contig_ends(name, seq, a.min_telomere, name_arms=not a.no_name_arms, arm_map=amap, want_probes=bool(a.chm13_genome))
    log(f'{a.input}: {n} records, {len(ends)} {"windows" if a.windows else "telomeric ends"}')
    _run(a, ends, n, 'windows' if a.windows else 'contig')


def cmd_reads(a):
    from . import query as Q
    log = _log(a.quiet); ends = []; n = 0
    it = Q.iter_bam(a.input) if a.input.endswith(('.bam', '.cram')) else Q.iter_records(a.input)
    for name, comment, seq in it:
        n += 1
        if len(seq) < a.min_length: continue
        e = Q.read_end(name, seq, a.min_telomere)
        if e: ends.append(e)
        if a.max_reads and len(ends) >= a.max_reads: break
    log(f'{a.input}: {n} reads scanned, {len(ends)} telomeric (experimental mode)')
    a.chm13_genome = None
    _run(a, ends, n, 'reads')


def cmd_evidence(a):
    import gzip, json
    from . import readevidence as RE, refs
    log = _log(a.quiet)
    crit = RE.load_criteria(a.criteria)
    th = crit['thresholds']
    ends_fa = a.chm13_ends or os.path.join(a.refs or refs.default_dir(), refs.CHM13_FA)
    if not os.path.exists(ends_fa): raise SystemExit(f'{ends_fa} not found: run subtel-atlas fetch-refs or give --chm13-ends')
    ends = RE.Chm13Ends(ends_fa)
    if a.claim:
        cat = a.catalogue or refs.repo_catalogue_dir()
        if not cat: raise SystemExit('--claim needs the catalogue tables (ends.tsv, mosaic_blocks.tsv): give --catalogue DIR')
        cl = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=[a.claim])
        if not cl: raise SystemExit(f'{a.claim}: not a non-native end of the catalogue')
        claim = cl[0]
    else:
        if not (a.donor and a.recipient and a.d_range):
            raise SystemExit('give --claim END_ID, or --donor, --recipient and --d-range (plus --switch-sub/--p-switch)')
        d_lo, d_hi = (int(x) for x in a.d_range.split('-'))
        claim = dict(claim_id=a.claim_id or f'{a.recipient}>{a.donor}', end_id=None, R=a.recipient, D=a.donor,
                     P=a.proximal or a.recipient, recipient_arm=a.recipient[3:], pattern=a.pattern, switch_sub=a.switch_sub,
                     switch_sub_distal=a.switch_sub_distal, p_switch=a.p_switch, d_lo=d_lo, d_hi=d_hi,
                     window_sub_len=a.window_sub_len, core_structure='', contig=a.contig)
    side = a.side or ('L' if claim['R'].endswith('p') else 'R')
    preset = a.preset
    log(f"claim {claim['claim_id']}: D={claim['D']} P={claim['P']} R={claim['R']} pattern={claim['pattern']} "
        f"switch_sub={claim.get('switch_sub')} | criteria {crit['criteria_version']} ({crit['status']}, {RE.criteria_hash(crit)})")
    reads = RE.iter_reads(a.input, region=a.region, reference=a.reference, index=a.index, min_len=a.min_length)
    sam = (a.output + '.chm13.sam') if a.evidence_bam else None
    win = None
    if a.assembly_window:
        import pysam
        with pysam.FastxFile(a.assembly_window) as fx:
            recs = {r.name: r.sequence for r in fx}
        win = recs.get(claim.get('end_id') or '') or next(iter(recs.values()))
    rows, doc = RE.run_claim(claim, reads, ends, crit, a.genome, preset=preset, threads=a.threads, side=side, dataset=a.dataset,
                             platform=a.platform, mask_cache=a.mask_cache, minimap2=a.minimap2, sam_out=sam, assembly_window=win, log=log)
    with open(a.output + '.reads.tsv', 'w') as f: RE.write_tsv(rows, f)
    v = RE.claim_verdict(rows, crit, same_hap_key='contig' if a.phased else None)
    v.update(claim_id=claim['claim_id'], criteria_sha=RE.criteria_hash(crit)); doc['claim_verdict'] = v
    with gzip.open(a.output + '.json.gz', 'wt') as f: json.dump(doc, f, indent=None, sort_keys=True, default=str)
    if sam:
        import pysam
        keep = {r['_key']: r for r in rows}
        with pysam.AlignmentFile(sam) as fi, pysam.AlignmentFile(a.output + '.chm13.unsorted.bam', 'wb', template=fi) as fo:
            for rec in fi:
                k, kind = rec.query_name.rsplit('|', 1)
                if kind != 'w' or k not in keep: continue
                r = keep[k]; rec.query_name = r['read']; rec.set_tag('XV', r['verdict']); rec.set_tag('XC', claim['claim_id'])
                fo.write(rec)
        pysam.sort('-o', a.output + '.chm13.bam', a.output + '.chm13.unsorted.bam'); pysam.index(a.output + '.chm13.bam')
        os.remove(a.output + '.chm13.unsorted.bam'); os.remove(sam)
    log(f"{len(rows)} reads: " + ', '.join(f"{k}={sum(r['verdict'] == k for r in rows)}" for k in
                                         ('conclusive', 'consistent', 'contradicts', 'non-informative')) + f" -> claim {v['verdict']}")


def main(argv=None):
    p = argparse.ArgumentParser(prog='subtel-atlas', description=DESC, epilog=EPI, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--version', action='version', version=f'subtel-atlas {__version__} (catalogue {CATALOGUE_VERSION})')
    sp = p.add_subparsers(dest='cmd', required=True, metavar='command')
    f = sp.add_parser('fetch-refs', help='download and checksum-verify reference data',
                      description='Create a reference directory: the 48 CHM13 v2.0 500 kb chromosome ends (HTTP range reads from the public '
                      'CHM13 v2.0 FASTA, or a local copy with --chm13) and the catalogue tables and representative windows (repository '
                      'checkout or GitHub raw URL). Every file is SHA-256-checked against the checksums shipped with the package.')
    f.add_argument('--dir', help='reference directory (default $SUBTEL_ATLAS_REFS or ~/.cache/subtel-atlas/v0.1)')
    f.add_argument('--chm13', metavar='FASTA', help='local CHM13 v2.0 FASTA (plain or gzip) instead of range reads')
    f.add_argument('--catalogue', metavar='DIR', help='local data/catalogue directory')
    from .refs import CHM13_URL, RAW_URL
    f.add_argument('--chm13-url', default=CHM13_URL, help=f'uncompressed CHM13 v2.0 FASTA with .fai (default {CHM13_URL})')
    f.add_argument('--catalogue-url', default=RAW_URL, help=f'base URL of the catalogue files (default {RAW_URL})')
    f.add_argument('--force', action='store_true', help='re-fetch even if verified files exist')
    f.add_argument('--quiet', action='store_true')
    f.set_defaults(func=cmd_fetch)
    v = sp.add_parser('verify-refs', help='check a reference directory against the package checksums')
    v.add_argument('--dir'); v.set_defaults(func=cmd_verify)
    q = sp.add_parser('query', help='assign telomeric contig ends (or pre-cut windows) of an assembly',
                      description='Find telomeric contig ends (same telomere caller and 500 bp rule as the catalogue), cut the 100 kb '
                      'window (extended to keep 80 kb of subtelomere when the telomere boundary lies beyond 20 kb, as in the catalogue), '
                      'and report arm, CHM13 mosaic (per-base identity tracks, paralog-family labels), S5/S10/S20 states, the nearest '
                      'structural allele with its nativity class and evidence tier, the nearest haplogroup, identities, boundary '
                      'stability under a +/-500 bp boundary shift and novelty flags. One output row per end.',
                      epilog='Arms: chromosome-named contigs give the arm directly (L end = p, R end = q; assumes standard orientation); '
                      'otherwise --chm13-genome enables the catalogue probe rule, else the mosaic rule is used (flag arm_inferred; '
                      'catalogue accuracy 84%% at group level for unnamed ends without probes, see README).')
    q.add_argument('input', help='FASTA/FASTQ (optionally gzip): contigs, or with --windows tip-frame windows')
    q.add_argument('--windows', action='store_true', help='records are pre-cut telomere-first windows in tip frame (first base = contig tip); '
                   'the telomere is called on the window; header tokens arm=<arm> and side=L|R are used if present')
    q.add_argument('--chm13-genome', metavar='FASTA|MMI', help='full CHM13 v2.0 genome (FASTA or minimap2 asm5 .mmi) for the arm probe rule '
                   '(10 kb probes at 250 kb and 1 Mb from the tip; needs contigs >= 520 kb; ~8 GB RAM)')
    _common_query(q); _arm_opts(q); q.set_defaults(func=cmd_query)
    r = sp.add_parser('reads', help='EXPERIMENTAL: assign telomeric long reads',
                      description='EXPERIMENTAL (validated only on simulated reads; see tests/validation/VALIDATION.md). Reads with a '
                      'telomere at either end are oriented telomere-first and analysed like windows (100 kb, longer for telomeres beyond 20 kb). Arms are inferred by '
                      'the mosaic rule; read errors lower identities and can shift switch calls.')
    r.add_argument('input', help='FASTQ/FASTA (optionally gzip) or BAM/CRAM (needs pysam)')
    r.add_argument('--min-length', type=int, default=10000, help='skip reads shorter than this (default 10000)')
    r.add_argument('--max-reads', type=int, default=0, help='stop after this many telomeric reads (0 = all)')
    _common_query(r); r.set_defaults(func=cmd_reads)
    e = sp.add_parser('reads-evidence', help='read-level evidence for a non-native terminal block (assembly-free)',
                      description='Test single long reads against a claim that the telomere-adjacent block of a chromosome end '
                      'is closer to CHM13 donor end D than to the recipient arm: terminal telomere, donor alleles at paralog-diagnostic '
                      'units distal to a maximum-likelihood switch, proximal alleles beyond it, and a MAPQ 60 genome-wide anchor on the '
                      'recipient chromosome (criteria: evidence/CRITERIA.md; thresholds: criteria JSON). Writes <prefix>.reads.tsv, '
                      '<prefix>.json.gz (per-unit calls) and optionally <prefix>.chm13.bam.')
    e.add_argument('input', help='reads: BAM/CRAM (local path or URL; with --region) or FASTA/FASTQ')
    e.add_argument('--region', help='region of the BAM/CRAM (e.g. chr4_PATERNAL:192379018-192384017 or chr4:193569000-193574055)')
    e.add_argument('--reference', help='FASTA for CRAM decoding'); e.add_argument('--index', help='.bai/.crai (local copy) for a remote file')
    e.add_argument('--claim', metavar='END_ID', help='claim from the catalogue (non-native end, e.g. HG002#1#chr4_PATERNAL#R)')
    e.add_argument('--catalogue', metavar='DIR', help='data/catalogue directory (default: repository checkout)')
    e.add_argument('--donor', help='explicit claim: donor CHM13 end (e.g. chr3q)'); e.add_argument('--recipient', help='recipient CHM13 end (e.g. chr20q)')
    e.add_argument('--proximal', help='proximal comparator (default = recipient)'); e.add_argument('--d-range', help='donor block on D, tip frame (a-b)')
    e.add_argument('--switch-sub', type=int, help='assembly switch, subtelomere coordinate'); e.add_argument('--switch-sub-distal', type=int)
    e.add_argument('--p-switch', type=int, help='assembly switch on CHM13 P (tip frame)')
    e.add_argument('--pattern', choices=['DP', 'PDP'], default='DP'); e.add_argument('--window-sub-len', type=int, default=80000)
    e.add_argument('--claim-id'); e.add_argument('--contig', help='claimed contig (own-assembly alignment concordance)')
    e.add_argument('--side', choices=['L', 'R', 'p', 'q'], help='chromosome end of stored reads (default: from the recipient arm)')
    e.add_argument('--genome', required=True, help='CHM13 v2.0 FASTA or minimap2 map-ont index (anchor)')
    e.add_argument('--chm13-ends', help='chm13_ends500k.fa (default: reference directory from fetch-refs)')
    e.add_argument('--refs', help='reference directory (fetch-refs)'); e.add_argument('--mask-cache', help='cache directory for shared masks')
    e.add_argument('--criteria', help='criteria JSON (default: packaged criteria_v1.json)')
    e.add_argument('--preset', choices=['map-ont', 'map-hifi'], default='map-ont'); e.add_argument('--dataset'); e.add_argument('--platform')
    e.add_argument('--phased', action='store_true', help='reads are one haplotype-assignable set (claim-level contradiction rule applies)')
    e.add_argument('--assembly-window', metavar='FASTA', help="the person's tip-frame end window (record named by end_id, else the "
                   'first record): switch concordance against the assembly analysed with the same unit model (criteria change C2)')
    e.add_argument('--min-length', type=int, default=5000); e.add_argument('--evidence-bam', action='store_true')
    e.add_argument('--minimap2'); e.add_argument('-t', '--threads', type=int, default=1)
    e.add_argument('-o', '--output', required=True, help='output prefix'); e.add_argument('--quiet', action='store_true')
    e.set_defaults(func=cmd_evidence)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
