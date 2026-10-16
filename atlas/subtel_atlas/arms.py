"""Home arm of an end (catalogue step 2b; used by pipeline/catalogue_05_arms.py and the query tool).

Rules in order: (1) contig name (chromosome-named contigs; L end = p arm, R end = q arm) -> 'named';
(2) 10 kb probes at 250 kb and 1 Mb from the tip mapped to the full CHM13 v2.0 genome (minimap2 asm5, primary alignments):
    informative if MAPQ >= 20, aligned >= 5 kb and orientation consistent with the arm; high = two informative probes agree (at
    group level), medium = one, conflict = disagree (highest MAPQ x aligned length wins);
(3) mosaic fallback: the CHM13 end with the largest summed alignment score over bins starting >= 50 kb from the boundary (all bins if
    the subtelomere is shorter); 'mosaic_group' inside a pooled group, else 'low'; 'none' without any aligned bin.
Catalogue arm = the arm or its pooled group (acroP, PAR1, PAR2).
"""
import re
import numpy as np
from .common import arm_of, group_of

PROBES = [(250_000, 10_000), (1_000_000, 10_000)]
_CHR = re.compile(r'(?:^|[^A-Za-z0-9])chr([1-9]|1[0-9]|2[0-2]|X|Y)(?=$|[^0-9A-Za-z]|_MATERNAL|_PATERNAL)')


def arm_from_name(contig, side):
    """'chr20' / 'chr20_MATERNAL' / 'HG002#1#chr20' + side R -> '20q'; None if the name carries no single chromosome."""
    m = {x for x in _CHR.findall(contig)}
    if len(m) != 1: return None
    return m.pop() + ('p' if side == 'L' else 'q')


def mosaic_arm(g, s0):
    """g: tiled-bin rows of one end (no fixed windows). Returns (arm, score fraction) or None."""
    prox = g[g.qstart - s0 >= 50000]; prox = prox if len(prox) else g; sc = {}
    for c in prox.cand.dropna():
        for x in c.split(';'):
            t, a, _ = x.split(':'); sc[t] = sc.get(t, 0) + int(a)
    if sc:
        b = max(sorted(sc), key=sc.get); return arm_of(b), round(sc[b] / sum(sc.values()), 3)
    return None


def probe_regions(contig_len, side):
    """Probe spans [a, b) on the contig (0-based) for one end; a probe is used only if it ends inside the tip-side half."""
    out = []
    for off, ln in PROBES:
        if off + ln > contig_len // 2: continue
        out.append((off, (off, off + ln) if side == 'L' else (contig_len - off - ln, contig_len - off)))
    return out


def best_probe_hit(best, q, aln, mq, chrom, ts, te, strand):
    """Keep, per probe q, the alignment with the largest aligned length x (MAPQ + 1) (catalogue rule)."""
    if q not in best or aln * (mq + 1) > best[q][0]: best[q] = (aln * (mq + 1), chrom, ts, te, strand, mq, aln)


def probe_row(e, off, side, hit, clen):
    sc, chrom, ts, te, st, mq, aln = hit
    L = clen[chrom]; mid = (ts + te) / 2; arm = chrom[3:] + ('p' if mid < L / 2 else 'q')
    tipst = st if side == 'L' else ('+' if st == '-' else '-')
    consistent = (arm.endswith('p') and tipst == '+') or (arm.endswith('q') and tipst == '-')
    dist = mid if arm.endswith('p') else L - mid
    return dict(end_id=e, offset=int(off), arm=arm, mapq=mq, aln=aln, consistent=consistent, dist_from_end=int(dist),
                informative=mq >= 20 and aln >= 5000 and consistent, score=mq * aln)


def probe_consensus(inf):
    """inf: DataFrame of informative probe rows of one end -> (probe arm, probe confidence)."""
    if not len(inf): return '', ''
    groups = set(inf.arm.map(group_of)); parm = inf.sort_values(['score', 'offset']).arm.iloc[-1]
    return parm, (('high' if len(inf) >= 2 else 'medium') if len(groups) == 1 else 'conflict')


def decide(contig_arm, parm, conf, ma):
    """(home arm, arm source, confidence) from the name, probe and mosaic calls."""
    if isinstance(contig_arm, str) and contig_arm:
        return contig_arm, 'contig_name', 'named'
    if parm:
        return parm, 'probe', conf
    return (ma, 'mosaic', ('mosaic_group' if group_of(ma) != ma else 'low')) if ma else ('unassigned', 'none', 'none')
