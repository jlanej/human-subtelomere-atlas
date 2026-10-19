"""Data behind the one-page pitch (report/PITCH.md): what a single showcase molecule looks like against the two CHM13
chromosome ends it is claimed to join.

Inputs (all committed or public):
  evidence/showcase/<claim>/<claim>.chm13.bam               whole reads (sequence) and CHM13 v2.0 @SQ lengths
  evidence/showcase/<claim>/<claim>.reads.tsv               per-read telomere length and switch interval (criteria v1)
  evidence/showcase/<claim>/<claim>.segments.chm13.bam      reads split at the switch, each part aligned to CHM13
  evidence/showcase/<claim>/<claim>.diagnostic_sites.chm13.bed   the claim model (sites where donor and recipient differ)
  CHM13 v2.0 (public, range reads with samtools faidx): outer 250 kb of each donor and recipient end

Outputs (report/pitch_data/, small and committed):
  profiles.tsv   per read and 2 kb window (1 kb step) of the read, starting at its telomeric end: differences per 1,000
                 aligned letters against the donor end and against the recipient end (minimap2 map-ont via mappy; best
                 primary hit; NaN = no alignment)
  sites.tsv      per read and diagnostic site covered by a split-read segment: position in the read (tip frame), donor
                 letter, recipient letter, read letter, call (D / P / other)
  segments.tsv   per read: where each part of the split read aligns on CHM13 (chromosome, position, MAPQ)

Letters in the claim model are in the tip frame (telomere outward); for q-arm ends this is the reverse complement of
CHM13's forward strand, so read letters on q-arm references are complemented before comparison.

Usage: python report/pitch_data.py [--claims SLUG,...]   (cache: work/pitch/)
"""
import os, sys, argparse, subprocess
import numpy as np, pandas as pd, pysam, mappy as mp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import readevidence as RE  # noqa: E402
TH = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json'))['thresholds']
SHOW = os.path.join(REPO, 'evidence', 'showcase')
OUT = os.path.join(REPO, 'report', 'pitch_data')
CACHE = os.path.join(REPO, 'work', 'pitch')
CHM13 = 'https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa'
WIN, STEP, SPAN, END_BP = 2000, 1000, 100_000, 250_000
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s): return s.translate(COMP)[::-1]


def end_seq(end, lens):
    """outer END_BP of a CHM13 end in tip frame (position 0 = telomeric tip); end like 'chr15q'"""
    chrom, arm = end[:-1], end[-1]
    os.makedirs(CACHE, exist_ok=True)
    fa = os.path.join(CACHE, f'{end}.fa')
    if not os.path.exists(fa):
        L = lens[chrom]
        reg = f'{chrom}:1-{END_BP}' if arm == 'p' else f'{chrom}:{L - END_BP + 1}-{L}'
        with open(fa + '.tmp', 'w') as fh:
            subprocess.run(['samtools', 'faidx', CHM13, reg], stdout=fh, check=True, cwd=CACHE)
        os.replace(fa + '.tmp', fa)
    s = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
    return s if arm == 'p' else rc(s)


def tip_frame(seq, platform, telo_len):
    """Tip-frame orientation (telomeric end first) exactly as the release (atlas readevidence.orient, criteria v1 telomere
    motifs incl. the ONT R9 basecalling artefacts): the orientation whose terminal telomere block reproduces the release's
    telo_len. A simple CCCTAA/TTAGGG count mis-orients R9 reads with artefact-rich telomeres."""
    hits = []
    for t in (seq, rc(seq)):
        b = RE.find_block(t[:70000], motifs=RE.telomere_motifs(TH, platform), maxjunk=max(TH['telomere_max_offset'] + 1, 1))
        if b and int(b['telo_len']) == int(telo_len): hits.append(t)
    assert len(hits) == 1, ('tip frame not reproduced', len(hits), telo_len)
    return hits[0]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--claims', default=None); a = ap.parse_args()
    loci = pd.read_csv(os.path.join(SHOW, 'showcase_loci.tsv'), sep='\t')
    if a.claims: loci = loci[loci.slug.isin(a.claims.split(','))]
    os.makedirs(OUT, exist_ok=True)
    prof, sites, segs = [], [], []
    for c in loci.itertuples():
        d = os.path.join(SHOW, c.slug)
        R = pd.read_csv(os.path.join(d, f'{c.slug}.reads.tsv'), sep='\t')
        with pysam.AlignmentFile(os.path.join(d, f'{c.slug}.chm13.bam')) as b:
            lens = dict(zip(b.references, b.lengths))
            reads = {}
            for r in b.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_name in reads: continue
                s = r.query_sequence
                reads[r.query_name] = rc(s) if r.is_reverse else s
        Dseq, Rseq = end_seq(c.D, lens), end_seq(c.R, lens)
        AD, AR = mp.Aligner(seq=Dseq, preset='map-ont'), mp.Aligner(seq=Rseq, preset='map-ont')
        for x in R.itertuples():
            t = tip_frame(reads[x.read_id], x.platform, x.telo_len)
            for st in range(0, min(len(t), SPAN) - WIN + 1, STEP):
                w = t[st:st + WIN]; v = []
                for A in (AD, AR):
                    h = [h for h in A.map(w) if h.is_primary]
                    v.append(1000 * (1 - max(hh.mlen / hh.blen for hh in h)) if h else np.nan)
                prof.append(dict(slug=c.slug, read_id=x.read_id, dataset=x.dataset, platform=x.platform, read_len=x.read_len,
                                 telo_len=x.telo_len, switch_lo=x.switch_ci_read_lo, switch_hi=x.switch_ci_read_hi,
                                 window_mid=st + WIN // 2, diff_per_kb_D=v[0], diff_per_kb_R=v[1]))
        # diagnostic sites read through the split-read segments
        model = {}
        for l in open(os.path.join(d, f'{c.slug}.diagnostic_sites.chm13.bed')):
            if l.startswith('track') or '_site' not in l: continue
            ch, s0, _, name = l.split('\t')[:4]
            model[(ch, int(s0))] = (name.split(':')[0], name.split('D=')[1][0], name.split('P=')[1][0])
        q_arm = {c.D[:-1]: c.D[-1] == 'q', c.R[:-1]: c.R[-1] == 'q'}
        T = pd.read_csv(os.path.join(d, f'{c.slug}.segments.tsv'), sep='\t')
        with pysam.AlignmentFile(os.path.join(d, f'{c.slug}.segments.chm13.bam')) as b:
            for r in b.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary: continue
                rid, seg = r.query_name.rsplit(':', 1)
                segs.append(dict(slug=c.slug, read_id=rid, segment=seg, chrom=r.reference_name, pos=r.reference_start + 1,
                                 end=r.reference_end, mapq=r.mapping_quality, aligned_bp=r.query_alignment_length))
                lo = int(T[(T.read_id == rid) & (T.segment == seg)].tip_frame_read_interval.iloc[0].split('-')[0])
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    k = (r.reference_name, rp)
                    if k not in model: continue
                    kind, Dl, Pl = model[k]
                    base = r.query_sequence[qp]
                    if q_arm.get(r.reference_name): base = base.translate(COMP)
                    # tip-frame position of this base in the read: segment start + offset along the segment
                    off = qp if not r.is_reverse else r.query_length - 1 - qp
                    sites.append(dict(slug=c.slug, read_id=rid, segment=seg, chrom=r.reference_name, chm13_pos=rp + 1,
                                      site_type=kind, read_pos=lo + off, donor_letter=Dl, recipient_letter=Pl,
                                      read_letter=base, call='D' if base == Dl else 'P' if base == Pl else 'other'))
        print(c.slug, 'reads', len(R), flush=True)
    pd.DataFrame(prof).to_csv(os.path.join(OUT, 'profiles.tsv'), sep='\t', index=False, float_format='%.2f')
    pd.DataFrame(sites).sort_values(['slug', 'read_id', 'read_pos']).to_csv(os.path.join(OUT, 'sites.tsv'), sep='\t', index=False)
    pd.DataFrame(segs).to_csv(os.path.join(OUT, 'segments.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
