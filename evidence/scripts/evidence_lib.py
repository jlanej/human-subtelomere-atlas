"""Shared helpers for the read-evidence resource (evidence/claims.tsv, evidence/reads.tsv.gz, evidence BAM/CRAM set, replication).

Nothing here changes the frozen v1 decision logic (atlas/subtel_atlas/readevidence.py, criteria hash a89c5eaf023ed032). It adds:
  * own-assembly placement of a read (which of the person's own assembly end windows it aligns to best);
  * the reporting-only 'definitive' tier (EVIDENCE TIERS below);
  * public-source helpers: end windows cut from the public assemblies by range (verified against the catalogue SHA-256), the
    samtools command that retrieves one read from its public alignment file, CHM13 coordinate conversion.

EVIDENCE TIERS (reporting only; defined after the population diagnostics, evidence/population/DIAGNOSTICS.md; the tier can only
remove claims from the v1 'confirmed' set, never add one):
  definitive read  = v1 'conclusive' AND switch-concordant AND unique-anchor tier (>= 1 kb of anchor on CHM13 R positions that no
                     other CHM13 end covers at >= 95%) AND own-assembly placement on the claimed end (the read aligns with more
                     matching bases to the claimed end window than to any other end window of the same person).
  definitive claim = v1 'confirmed' AND >= 2 definitive reads from distinct molecules (HG002: non-derived reads after
                     molecule de-duplication); 'multi-dataset' when they come from >= 2 independent datasets.
"""
import os, re, hashlib, subprocess, collections

DEFINITIVE_MIN_READS = 2
UNIQUE_TIER_BP = 1000          # = criteria_v1 thresholds.anchor_unique_tier_bp (reported, not changed)
CHM13_FA_URL = 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz'
CHM13_PLAIN_URL = 'https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa'
Q100_BASE = 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/'
ONT1KG = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/{s}.t2t.cram'


def b(x):
    return x in (True, 'True', 'true', '1', 1)


def f(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def https(url):
    """s3://human-pangenomics/... -> public HTTPS URL"""
    return url.replace('s3://human-pangenomics/', 'https://human-pangenomics.s3.amazonaws.com/') if isinstance(url, str) else url


def chm13_coord(chrom_len, end, x):
    """1-based CHM13 v2.0 chromosome coordinate of tip-frame position x (0 = telomeric end) on CHM13 end `end` (e.g. chr20q)"""
    if x is None or x == '' or (isinstance(x, float) and x != x): return None
    x = int(float(x)); chrom = end[:-1]
    return x + 1 if end.endswith('p') else chrom_len[chrom] - x


def chm13_interval(chrom_len, end, lo, hi):
    """tip-frame interval [lo, hi] on `end` -> (chrom, start, end) 1-based inclusive in CHM13 coordinates"""
    a, c = chm13_coord(chrom_len, end, lo), chm13_coord(chrom_len, end, hi)
    if a is None or c is None: return None
    return end[:-1], min(a, c), max(a, c)


def rc(s):
    return s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


# ------------------------------------------------------------------------------------------------------------- own placement
def placement(seqs, windows, preset='map-ont', top=4):
    """Own-assembly placement. seqs: {read: sequence}; windows: {end_id: tip-frame window sequence} (all ends of one person).
    Score per window = sum of matching bases (mappy mlen) over all hits of the read to that window (best_n=10). Returns
    {read: [(end_id, matches, aligned_block_len, identity), ...]} sorted by matches (top `top`). Same rule as
    evidence/scripts/population_diagnostics.py (conclusive_read_placement.tsv)."""
    import mappy as mp, tempfile
    td = tempfile.mkdtemp(prefix='place_'); fa = os.path.join(td, 'w.fa')
    with open(fa, 'w') as fh:
        for e, s in windows.items(): fh.write(f'>{e}\n{s}\n')
    al = mp.Aligner(fa, preset=preset, best_n=10); out = {}
    for n, s in seqs.items():
        agg = collections.defaultdict(lambda: [0, 0])
        for h in al.map(s):
            agg[h.ctg][0] += h.mlen; agg[h.ctg][1] += h.blen
        out[n] = sorted([(c, m, bl, round(m / bl, 4)) for c, (m, bl) in agg.items()], key=lambda x: -x[1])[:top]
    os.remove(fa); os.rmdir(td)
    return out


def placement_fields(pl, claimed):
    """summary of a placement list for one read against the claimed end id(s)"""
    claimed = set(claimed)
    if not pl:
        return dict(own_best_end='', own_best_matches=0, own_best_identity=None, own_claimed_matches=0, own_claimed_identity=None,
                    own_runner_up_end='', own_runner_up_matches=0, own_margin_matches=None, own_on_claimed_end=False)
    best = pl[0]; cl = next((p for p in pl if p[0] in claimed), None)
    ru = next((p for p in pl[1:] if p[0] != best[0]), None)
    # best non-claimed competitor (the end the read would place on if it were not from the claimed end)
    comp = next((p for p in pl if p[0] not in claimed), None)
    on = best[0] in claimed and (comp is None or best[1] > comp[1])
    return dict(own_best_end=best[0], own_best_matches=best[1], own_best_identity=best[3],
                own_claimed_matches=cl[1] if cl else 0, own_claimed_identity=cl[3] if cl else None,
                own_runner_up_end=comp[0] if comp else '', own_runner_up_matches=comp[1] if comp else 0,
                own_margin_matches=(cl[1] - (comp[1] if comp else 0)) if cl else None, own_on_claimed_end=bool(on))


# ------------------------------------------------------------------------------------------------------------- tiers
def read_tier(r):
    """tier of one finalised read row (after placement fields are attached)"""
    v = r.get('v1_verdict', r.get('verdict'))
    if v != 'conclusive': return v or 'non-informative'
    if not b(r.get('switch_concordant')): return 'conclusive_discordant'
    uniq = (f(r.get('anchor_unique_bp')) or 0) >= UNIQUE_TIER_BP
    on = b(r.get('own_on_claimed_end'))
    if uniq and on: return 'definitive'
    if on: return 'conclusive_shared_anchor'
    if uniq: return 'conclusive_unique_anchor_off_end'
    return 'conclusive_shared_anchor_off_end'


TIER_ORDER = ['definitive', 'conclusive_shared_anchor', 'conclusive_unique_anchor_off_end', 'conclusive_shared_anchor_off_end',
              'conclusive_discordant', 'consistent', 'contradicts', 'non-informative']


def definitive_claim(v1_verdict, def_rows, molecule_key=lambda r: r.get('molecule') or r.get('read_id')):
    """(definitive_verdict, n_definitive_molecules, datasets)"""
    mols = {molecule_key(r) for r in def_rows}
    ds = sorted({r.get('dataset') for r in def_rows if r.get('dataset')})
    if v1_verdict == 'confirmed' and len(mols) >= DEFINITIVE_MIN_READS:
        return ('definitive_multi_dataset' if len(ds) >= 2 else 'definitive'), len(mols), ds
    if v1_verdict == 'confirmed':
        return 'confirmed_not_definitive', len(mols), ds
    return 'not_definitive', len(mols), ds


# ------------------------------------------------------------------------------------------------------------- public sources
def window_from_public(end_row, cwd=None):
    """Cut a catalogue end window from the public assembly FASTA by range (samtools faidx; bgzip + .fai/.gzi on the server),
    in tip frame, and verify it against ends.tsv window_sha256. end_row: mapping with assembly_url, contig, win_start, win_end,
    strand, window_sha256."""
    url = https(end_row['assembly_url'])
    reg = f"{end_row['contig']}:{int(end_row['win_start'])}-{int(end_row['win_end'])}"
    p = subprocess.run(['samtools', 'faidx', url, reg], capture_output=True, text=True, cwd=cwd)
    if p.returncode: raise RuntimeError(f'samtools faidx {url} {reg}: {p.stderr[-300:]}')
    s = ''.join(p.stdout.split('\n')[1:])
    s = s if end_row['strand'] == '+' else rc(s)
    h = hashlib.sha256(s.encode()).hexdigest()
    if h != end_row['window_sha256']: raise RuntimeError(f"{end_row['end_id']}: window SHA-256 mismatch ({h[:12]})")
    return s


def index_url(url):
    return url + ('.crai' if url.endswith('.cram') else '.bai')


def samtools_cmd(url, region, read, cram=False, ref='chm13v2.0.fa.gz'):
    """exact command that prints this read's primary record from its public alignment file (bash)"""
    t = f"-T {ref} " if cram else ''
    return (f"samtools view -F 0x900 {t}-X '{url}' '{index_url(url)}' '{region}' | awk -v r='{read}' '$1==r'")


def molecule_keys(name):
    """ONT read UUIDs (duplex 'a;b' -> both), PacBio movie/ZMW (as evidence/scripts/hg002_multiplatform.py)"""
    if name.startswith('m') and '/' in name: return {'/'.join(name.split('/')[:2])}
    return set(re.split(r'[;,]', name))


def sha256(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''): h.update(blk)
    return h.hexdigest()
