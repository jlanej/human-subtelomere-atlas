"""Constants and small helpers shared by the catalogue pipeline (pipeline/catalogue_common.py re-exports them) and the query tool."""
import gzip, hashlib

WINDOW = 100_000                             # telomere-anchored window length (bp)

_comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(_comp)[::-1]

def readfa(fn):
    """Yield (name, sequence) from a FASTA file (optionally gzip-compressed); name = first word of the header."""
    op = gzip.open if fn.endswith('.gz') else open
    name = None; buf = []
    with op(fn, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name is not None: yield name, ''.join(buf)
                name = line[1:].strip().split()[0]; buf = []
            else: buf.append(line.strip())
    if name is not None: yield name, ''.join(buf)

def sha256_str(s): return hashlib.sha256(s.encode()).hexdigest()
def sha256_file(fn, bs=1 << 20):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(bs), b''): h.update(b)
    return h.hexdigest()

def arm_of(chm13_end):
    """'chr20q' -> '20q'"""
    return chm13_end[3:] if isinstance(chm13_end, str) and chm13_end.startswith('chr') else chm13_end

# Catalogue arm groups: ends whose distal sequence is shared between chromosomes are catalogued as pooled groups.
#   acroP = 13p, 14p, 15p, 21p, 22p (acrocentric short arms); PAR1 = Xp, Yp; PAR2 = Xq, Yq. All other arms are their own group.
ARM_GROUPS = {'acroP': ['13p', '14p', '15p', '21p', '22p'], 'PAR1': ['Xp', 'Yp'], 'PAR2': ['Xq', 'Yq']}
_G = {a: g for g, arms in ARM_GROUPS.items() for a in arms}
def group_of(arm): return _G.get(arm, arm) if isinstance(arm, str) and arm else 'unassigned'
def group_ends(group):
    """CHM13 end names ('chr13p', ...) that count as 'home' for a catalogue arm group"""
    if group == 'unassigned': return []
    return ['chr' + a for a in ARM_GROUPS.get(group, [group])]
DEFINING_CONFIDENCE = {'named', 'high', 'medium', 'mosaic_group'}
CHM13_ENDS = [f'chr{c}{s}' for c in list(range(1, 23)) + ['X', 'Y'] for s in 'pq']   # the 48 CHM13 v2.0 ends, reference order
