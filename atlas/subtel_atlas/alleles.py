"""Structural alleles = recurring core mosaic structures (catalogue step 2f; used by pipeline/catalogue_08_alleles.py and the query tool).

Core structure: blocks shorter than MIN_BLOCK (10 kb, subtelomere coordinates) are dropped (their span split at its midpoint between
the neighbours) and adjacent blocks with the same source merged; structure = block sources telomere -> centromere ('home', CHM13 arm
of a foreign block, 'none'); switch vector = core block boundaries.
Definition: within (arm, structure), complete-linkage clustering of switch vectors on the Chebyshev distance cut at SW_TOL (10 kb).
Matching an end to defined alleles (match_allele): full-length ends - same structure, same number of switches, every switch within
SW_TOL of the allele's median switch; partial ends - the allele truncated at the observed subtelomere length must equal the end.
Several candidates -> the allele whose median switches are nearest (smallest maximum switch difference), then the allele with most
defining ends, then the ID. (Changed before the v0.1 publication: most-frequent-first absorbed ends of same-structure alleles < 10 kb
apart into the larger allele; TIEBREAK_LEGACY reproduces the old order for comparison.)
"""
import numpy as np
from .common import arm_of

MIN_BLOCK, SW_TOL, FULL_SUB = 10000, 10000, 80_000
DEF_SETS = {'HPRC_R2', 'CHM13', 'HG002_Q100'}


def core(bl, min_block=MIN_BLOCK):
    """bl: raw blocks (objects with is_home, source, sub_start, sub_end) in order. Returns (core blocks, structure, switches)."""
    bl = [dict(source='home' if b.is_home else b.source, hsrc=b.source, s=float(b.sub_start), e=float(b.sub_end)) for b in bl]
    while len(bl) > 1:
        L = [b['e'] - b['s'] for b in bl]; i = int(np.argmin(L))
        if L[i] >= min_block: break
        if i == 0: bl[1]['s'] = bl[0]['s']
        elif i == len(bl) - 1: bl[-2]['e'] = bl[-1]['e']
        else:
            m = (bl[i]['s'] + bl[i]['e']) / 2; bl[i - 1]['e'] = m; bl[i + 1]['s'] = m
        del bl[i]
        j = 1
        while j < len(bl):
            if bl[j]['source'] == bl[j - 1]['source']: bl[j - 1]['e'] = bl[j]['e']; del bl[j]
            else: j += 1
    toks = [b['source'] if b['source'] in ('home', 'none') else arm_of(b['source']) for b in bl]
    return bl, '|'.join(toks), [int(round(b['s'])) for b in bl[1:]]


def prefix_compatible(struct, sw, a_struct, a_sw, Lobs, min_block=MIN_BLOCK, sw_tol=SW_TOL):
    """end (struct, sw) observed up to Lobs; allele (a_struct, a_sw median). Allele truncated at Lobs must equal the end."""
    at = a_struct.split('|'); keep = [s for s in a_sw if s < Lobs - min_block]
    if struct.split('|') != at[:len(keep) + 1]: return False
    return len(sw) == len(keep) and all(abs(a - b) <= sw_tol for a, b in zip(sw, keep))


def cluster_structure(ids, switches, sw_tol=SW_TOL):
    """Alleles of one (arm, structure): ids and their switch vectors -> list of dict(members, sw, swmin, swmax, medoid)."""
    from scipy.cluster.hierarchy import linkage, fcluster      # catalogue build only (not needed for queries)
    X = np.array([s for s in switches], float).reshape(len(ids), -1); cl = []
    if len(ids) == 1 or X.shape[1] == 0: lab = np.ones(len(ids), int)
    else: lab = fcluster(linkage(X, 'complete', metric='chebyshev'), t=sw_tol, criterion='distance')
    for l in sorted(set(lab)):
        ix = np.where(lab == l)[0]; Xi = X[ix]
        if len(ix) > 1 and Xi.shape[1]:
            D = np.abs(Xi[:, None, :] - Xi[None, :, :]).max(2).sum(1); med = sorted(range(len(ix)), key=lambda j: (D[j], ids[ix[j]]))[0]
        else: med = sorted(range(len(ix)), key=lambda j: ids[ix[j]])[0]
        cl.append(dict(members=[ids[j] for j in ix], sw=np.median(Xi, 0).round().astype(int).tolist() if Xi.shape[1] else [],
                       swmin=Xi.min(0).astype(int).tolist() if Xi.shape[1] else [], swmax=Xi.max(0).astype(int).tolist() if Xi.shape[1] else [],
                       medoid=ids[ix[med]]))
    return cl


def max_switch_delta(a_sw, sw): return max((abs(a - b) for a, b in zip(a_sw, sw)), default=0)


def match_allele(structure, switches, sub_len, full, alleles, min_block=MIN_BLOCK, sw_tol=SW_TOL, legacy=False):
    """Candidates among defined alleles of the end's arm. alleles: dicts with allele_id, structure, sw (median switches), n (defining ends).
    Returns (best allele or None, candidate list)."""
    if full:
        cand = [c for c in alleles if c['structure'] == structure and len(c['sw']) == len(switches)
                and all(abs(a - b) <= sw_tol for a, b in zip(switches, c['sw']))]
    else:
        cand = [c for c in alleles if prefix_compatible(structure, switches, c['structure'], c['sw'], sub_len, min_block, sw_tol)]
    if not cand: return None, []
    best = sorted(cand, key=lambda c: tiebreak_key(c, switches, legacy))[0]
    return best, cand


def tiebreak_key(c, switches, legacy=False):
    """nearest median switch vector first, then most defining ends, then ID (legacy=True: most defining ends first)"""
    d = max_switch_delta(c['sw'], switches)
    return (-c['n'], d, c['allele_id']) if legacy else (d, -c['n'], c['allele_id'])
