"""Shared paths and helpers for the cross-validation against Tan et al. 2026 (bioRxiv doi:10.64898/2026.09.03.748401).
Input  : the Zenodo archive subtelomere_500kb_2kbtelo.zip (doi:10.5281/zenodo.21807156; CC BY 4.0), path from $TAN_ZIP
         (default work/tan2026/subtelomere_500kb_2kbtelo.zip). One bgzipped FASTA per chromosome arm; each record is the
         500 kb nearest the telomere with the telomere trimmed to 2 kb. Record names:
         <sample>#<hap>#<cohort>-<assembler>_<contig>:<start>-<end>/<nc|rc>#chr<N>:<p|q>-arm
Outputs: intermediates under work/tan2026/ (git-ignored; the sequences are not redistributed), small tables under
         data/validation/, one figure under figures/validation/.
"""
import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import catalogue_common as C  # noqa: E402

REPO = C.REPO
TW = os.path.join(C.WORK, 'tan2026')
TIPS = os.path.join(TW, 'tips')
ZIP = os.environ.get('TAN_ZIP', os.path.join(TW, 'subtelomere_500kb_2kbtelo.zip'))
VAL = os.path.join(REPO, 'data', 'validation')
FIG = os.path.join(REPO, 'figures', 'validation')
THREADS = min(4, int(os.environ.get('ATLAS_THREADS', '4')))
SUB = 100_000          # subtelomere length classified per record (first 100 kb from the telomere boundary)
COHORTS_INDEP = ['HGSVC3', 'CPC', 'APC', 'JASA']
for d in (TW, TIPS, VAL, FIG): os.makedirs(d, exist_ok=True)

_NAME = re.compile(r'^(?P<sample>[^#]+)#(?P<hap>\d+)#(?P<cohort>[A-Za-z0-9]+)-(?P<asm>[a-z]+)_(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+)/'
                   r'(?P<ori>nc|rc)#chr(?P<chrom>[0-9XY]+):(?P<side>[pq])-arm$')

def ints(s):
    """'123,4567' -> [123, 4567]; empty / NaN -> []"""
    return [int(float(x)) for x in str(s).split(',') if x not in ('', 'nan')] if isinstance(s, str) else []

def load_alleles():
    """released structural alleles (data/catalogue/alleles_structural.tsv) grouped by catalogue arm, in the form used by
    subtel_atlas.alleles.match_allele (n = defining ends)"""
    import pandas as pd
    A = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'alleles_structural.tsv'), sep='\t', low_memory=False)
    by = {}
    for r in A.itertuples():
        by.setdefault(r.arm, []).append(dict(allele_id=r.allele_id, structure=r.structure, sw=ints(r.switches_median),
                                             n=int(r.n_defining), tier=r.tier))
    return A, by

def assign(arm, structure, sw, sub_len, by_arm, rule='nearest'):
    """Nearest structural allele with the catalogue's own matcher (subtel_atlas.alleles.match_allele, shared with the query
    tool): full length (subtelomere >= 80 kb) - same core structure, same number of switches, every switch within 10 kb of the
    allele median; partial - prefix rule. Several candidates -> nearest median switches, then most defining ends, then ID
    (rule='most_frequent' = the legacy order). Returns (allele_id or '', assignment, candidates)."""
    from subtel_atlas import alleles as AL
    full = sub_len >= AL.FULL_SUB
    best, cand = AL.match_allele(structure, sw, sub_len, full, by_arm.get(arm, []), legacy=(rule == 'most_frequent'))
    if not cand: return '', 'unassigned' if full else 'unassigned_partial', ''
    kind = 'full' if full else ('partial_unique' if len(cand) == 1 else 'partial_ambiguous')
    return best['allele_id'], kind, ','.join(c['allele_id'] for c in cand)

def core_of(blocks, clip=None):
    """core structure of one end from its mosaic blocks (DataFrame rows of one end with is_home, source, sub_start, sub_end,
    block) - subtel_atlas.alleles.core. clip: emulate a shorter window - blocks are cut at subtelomere coordinate `clip` before
    the core computation (blocks starting at or beyond it are dropped, the last block ends at it), so that blocks < 10 kb at
    the new distal end are merged exactly as in a catalogue window of that length."""
    from subtel_atlas import alleles as AL
    b = blocks.sort_values('block')
    if clip is not None:
        b = b[b.sub_start < clip].copy(); b['sub_end'] = b.sub_end.astype(float)
        if len(b): b.loc[b.index[-1], 'sub_end'] = min(float(b.sub_end.iloc[-1]), float(clip))
    bl, st, sw = AL.core(list(b.itertuples()))
    return st, sw, len(bl)

def truncate(structure, sw, L):
    """structure/switches restricted to the first L bp of subtelomere (switches closer than MIN_BLOCK to L dropped)"""
    from subtel_atlas import alleles as AL
    keep = [s for s in sw if s < L - AL.MIN_BLOCK]
    return '|'.join(structure.split('|')[:len(keep) + 1]), keep

def wilson(k, n, z=1.959964):
    if n == 0: return (float('nan'), float('nan'))
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (max(0.0, c - h), min(1.0, c + h))

def parse_name(n):
    m = _NAME.match(n)
    if not m: raise ValueError(n)
    d = m.groupdict(); d['arm'] = d['chrom'] + d['side']; d['hap'] = int(d['hap']); return d
