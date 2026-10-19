"""Data for the guided read tour on the website (docs/evidence/tour.html, igv.js), written to docs/evidence/tour/data/.

Inputs: evidence/claims.tsv (the definitive claims), the per-claim read sets of the v0.1.0 release assets
(subtel_read_evidence_v1_<cohort>.tar.gz, unpacked to work/tour/release/read_evidence_v1/; `--release DIR` to point
elsewhere), data/catalogue/{ends,mosaic_blocks}.tsv, and the public CHM13 v2.0 FASTA (range reads).

Output (docs/evidence/tour/data/):
  claims/<slug>/combo.fa.gz (+.fai, .gzi) one reference per claim: the person's own end windows (hap1/hap2, tip frame) and
                                         the CHM13 v2.0 donor end D and recipient end R in tip frame (position 1 = the telomeric
                                         tip; q-arm ends reverse-complemented), long enough for the longest molecule
  claims/<slug>/asm.bam (+.bai)           catalogue blocks of the own windows placed on the CHM13 end their label names
  claims/<slug>/model.json, model.bed    the formal claim model in tip frame (website/tour_models.py)
  claims/<slug>/parts.bam (+.bai)         up to MAX_READS definitive molecules of the claim, whole, cut at their switch interval
                                         (criteria v1, reads.tsv): <read>/part1 = telomeric end to the switch, aligned to D only;
                                         <read>/part2 = switch to the read end, aligned to R only (minimap2 map-ont / map-hifi)
  claims/<slug>/own.bam (+.bai), own_blocks.bed
                                         the same molecules on the person's own assembly windows (from the release asset,
                                         filtered to the shown reads) and the catalogue mosaic blocks of those windows
  extra/...                              chapter-specific subsets (whole molecule vs each end; reads by dataset)
  tour_data.json                         everything the page needs (reads, alignments, lengths, chapters)

Usage (from the repository root; needs minimap2, samtools, bgzip and network; cache $TOUR_CACHE or work/tour/):
  python website/tour_build.py [--release DIR]; python website/tour_models.py; python website/tour_build.py --sites-only
"""
import os, json, shutil, subprocess, argparse, re
import numpy as np, pandas as pd, pysam

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import readevidence as RE  # noqa: E402
TH = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json'))['thresholds']
OUT = os.path.join(REPO, 'docs', 'evidence', 'tour', 'data')
CACHE = os.environ.get('TOUR_CACHE', os.path.join(REPO, 'work', 'tour'))   # all temporary files (never deleted by this script)
CHM13 = 'https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa'
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
MAX_READS, DONOR_BP, R_PAD = 6, 260_000, 40_000
HERO_SLUG, HERO_READ = 'HG002_2_chr3_MATERNAL_R', '3788f60a-2f28-4b28-917b-f1ecae3645e3'
TECH_SLUG = 'HG002_2_chr19_MATERNAL_R'


def rc(s): return s.translate(COMP)[::-1]


def run(cmd, **kw): return subprocess.run(cmd, check=True, **kw)


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


def contig(end): return f'CHM13_{end}_end'


def fetch_end(end, lens, n):
    chrom, arm = end[:-1], end[-1]; L = lens[chrom]; n = int(min(n, L))
    fa = os.path.join(CACHE, 'chm13', f'{end}_{n}.fa'); os.makedirs(os.path.dirname(fa), exist_ok=True)
    if not os.path.exists(fa):
        reg = f'{chrom}:1-{n}' if arm == 'p' else f'{chrom}:{L - n + 1}-{L}'
        with open(fa + '.tmp', 'w') as fh: run(['samtools', 'faidx', CHM13, reg], stdout=fh, cwd=os.path.dirname(fa))
        os.replace(fa + '.tmp', fa)
    s = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
    return s if arm == 'p' else rc(s)


def write_fa(path, recs):
    with open(path, 'w') as fh:
        for nm, s in recs: fh.write(f'>{nm}\n' + ''.join(s[i:i + 80] + '\n' for i in range(0, len(s), 80)))


def minimap(target_fa, query_fa, preset, extra=()):
    p = run(['minimap2', '-a', '-x', preset, '--secondary=no', '-Y', '--MD', '-t', '6', *extra, target_fa, query_fa],
            capture_output=True, text=True)
    return p.stdout


def sam_to_bam(sam_texts, header, out, rg_of=None):
    """merge SAM texts (records only used) into one coordinate-sorted, indexed BAM with the given header"""
    tmp = os.path.join(CACHE, 'tmp', os.path.basename(out) + '.unsorted.bam'); os.makedirs(os.path.dirname(tmp), exist_ok=True); stats = []
    with pysam.AlignmentFile(tmp, 'wb', header=header) as fo:
        for txt in sam_texts:
            if not txt.strip(): continue
            path = os.path.join(CACHE, 'tmp', 'part.sam'); open(path, 'w').write(txt)
            with pysam.AlignmentFile(path) as f:
                for r in f:
                    if r.is_unmapped or r.is_secondary: continue
                    a = pysam.AlignedSegment(fo.header); a.query_name = r.query_name; a.flag = r.flag
                    a.reference_id = fo.header.get_tid(r.reference_name); a.reference_start = r.reference_start
                    a.mapping_quality = r.mapping_quality; a.cigar = r.cigar; a.query_sequence = r.query_sequence
                    a.query_qualities = r.query_qualities
                    tags = [t for t in r.get_tags() if t[0] in ('NM', 'MD', 'AS', 'SA', 'tp', 'de')]
                    if rg_of: tags.append(('RG', rg_of(r.query_name)))
                    a.set_tags(tags); fo.write(a)
                    if not r.is_supplementary:
                        stats.append(dict(name=r.query_name, ref=r.reference_name, start=r.reference_start + 1,
                                          end=r.reference_end, mapq=r.mapping_quality, aligned=r.query_alignment_length,
                                          qlen=r.infer_read_length(), nm=r.get_tag('NM')))
    run(['samtools', 'sort', '-o', out, tmp], capture_output=True); run(['samtools', 'index', out])
    return stats


def header(contigs, rgs=(), comment=''):
    h = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': l} for n, l in contigs]}
    if rgs: h['RG'] = [{'ID': r, 'SM': r} for r in rgs]
    if comment: h['CO'] = [comment]
    return h


def diffs_bed(D, R, Dseq, Rseq, path):
    """display track: single-base differences between the two CHM13 ends inside their alignments, on both contigs"""
    td, tr = os.path.join(CACHE, 'd.fa'), os.path.join(CACHE, 'r.fa')
    write_fa(td, [(contig(D), Dseq)]); write_fa(tr, [(contig(R), Rseq[:300_000])])
    p = run(['minimap2', '-c', '--cs', '-x', 'asm20', '-t', '6', td, tr], capture_output=True, text=True)
    rows, n_snv, sites = [], 0, []
    for l in p.stdout.splitlines():
        f = l.split('\t'); qn, ql, qs, qe, strand, tn, tl, ts, te, nmatch, alen = f[:11]
        qs, qe, ts, te, nmatch, alen = map(int, (qs, qe, ts, te, nmatch, alen))
        if alen < 2000 or nmatch / alen < 0.85: continue
        ident = nmatch / alen
        rows.append((contig(D), ts, te, f'shared with {R[:-1]} ({ident:.1%})', '170,170,170'))
        rows.append((contig(R), qs, qe, f'shared with {D[:-1]} ({ident:.1%})', '170,170,170'))
        cs = next(x[5:] for x in f[12:] if x.startswith('cs:Z:'))
        t, q = ts, (qs if strand == '+' else qe)
        for op in re.findall(r'(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+)', cs):
            if op[0] == ':':
                k = int(op[1:]); t += k; q += k if strand == '+' else -k
            elif op[0] == '*':
                tb, qb = op[1].upper(), op[2].upper()          # target (D) base; query (R) base on the target strand
                qpos = q if strand == '+' else q - 1
                r_own = qb if strand == '+' else qb.translate(COMP)    # the R contig's own base at qpos
                d_on_r = tb if strand == '+' else tb.translate(COMP)   # the D base expressed on the R contig's strand
                rows.append((contig(D), t, t + 1, f'{D[:-1]}:{tb} | {R[:-1]}:{qb}', '213,94,0'))
                rows.append((contig(R), qpos, qpos + 1, f'{R[:-1]}:{r_own} | {D[:-1]}:{d_on_r}', '0,114,178'))
                sites.append(dict(d_pos=t, d_base=tb, d_alt=qb, r_pos=qpos, r_base=r_own, r_alt=d_on_r))
                n_snv += 1; t += 1; q += 1 if strand == '+' else -1
            elif op[0] == '+':
                q += (len(op) - 1) if strand == '+' else -(len(op) - 1)
            else:
                t += len(op) - 1
    with open(path, 'w') as fh:
        fh.write('track name="CHM13 differences" itemRgb="On"\n')
        for c_, a, b, lab, col in sorted(rows, key=lambda x: (x[0], x[1])):
            fh.write(f'{c_}\t{a}\t{b}\t{lab}\t0\t.\t{a}\t{b}\t{col}\n')
    json.dump(sites, open(path.replace('.bed', '.sites.json'), 'w'))
    return n_snv


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--release', default=os.path.join(REPO, 'work', 'tour', 'release', 'read_evidence_v1'))
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True); os.makedirs(OUT, exist_ok=True)   # outputs are overwritten in place
    C = pd.read_csv(os.path.join(REPO, 'evidence', 'claims.tsv'), sep='\t', low_memory=False)
    C = C[C.definitive_verdict.astype(str).str.startswith('definitive')].copy()
    E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False, usecols=['end_id', 'telo_len'])
    MB = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'mosaic_blocks.tsv'), sep='\t', low_memory=False)
    MB = MB[MB.level == 'raw']
    lens = None
    claims, reads_all = [], {}
    # 1. molecules per claim
    for c in C.itertuples():
        d = os.path.join(a.release, c.cohort, c.slug)
        R = pd.read_csv(os.path.join(d, f'{c.slug}.reads.tsv'), sep='\t', low_memory=False)
        R = R[(R.tier == 'definitive') & (R.counted == True)].sort_values('read_len', ascending=False)  # noqa: E712
        pick = R.groupby('dataset', sort=False).head(1)
        pick = pd.concat([pick, R[~R.read_id.isin(pick.read_id)]]).head(MAX_READS)
        if c.slug == HERO_SLUG and HERO_READ not in set(pick.read_id):
            pick = pd.concat([R[R.read_id == HERO_READ], pick]).head(MAX_READS)
        if c.slug == TECH_SLUG:      # technology chapter: one molecule from every dataset
            pick = R.groupby('dataset', sort=False).head(1)
        with pysam.AlignmentFile(os.path.join(d, f'{c.slug}.chm13.bam')) as b:
            if lens is None: lens = dict(zip(b.references, b.lengths))
            seqs = {}
            for r in b.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_name in seqs: continue
                if r.query_name in set(pick.read_id): seqs[r.query_name] = rc(r.query_sequence) if r.is_reverse else r.query_sequence
        mols = []
        for x in pick.itertuples():
            t = tip_frame(seqs[x.read_id], x.platform, x.telo_len); lo, hi = int(x.switch_ci_read_lo), int(x.switch_ci_read_hi)
            mols.append(dict(read_id=x.read_id, dataset=x.dataset, platform=x.platform, read_len=len(t), telo_len=int(x.telo_len),
                             switch_lo=lo, switch_hi=hi, seq=t))
        claims.append(dict(slug=c.slug, claim_id=c.claim_id, cohort=c.cohort, individual=c.individual, arm=c.arm,
                           allele=c.allele_id, D=c.donor_D, P=c.proximal_P, R=c.recipient_R, haplotype=c.haplotype,
                           n_definitive=int(c.n_definitive), n_datasets=int(c.n_definitive_datasets),
                           verdict=c.definitive_verdict, release_dir=d, mols=mols))
    # 2. one reference per (D, R) pair
    pairs = {}
    for c in claims:
        k = (c['D'], c['R']); pairs.setdefault(k, []).append(c)
    refinfo = {}
    for (D, R), cl in pairs.items():
        key = f'{D}__{R}'; od = os.path.join(CACHE, 'refs_pairs', key); os.makedirs(od, exist_ok=True)   # cache only
        nR = max(m['read_len'] - m['switch_hi'] for c in cl for m in c['mols']) + R_PAD
        Dseq, Rseq = fetch_end(D, lens, DONOR_BP), fetch_end(R, lens, max(nR, 120_000))
        fa = os.path.join(od, 'ref.fa'); write_fa(fa, [(contig(D), Dseq), (contig(R), Rseq)])
        run(['bgzip', '-f', '-i', fa]); run(['samtools', 'faidx', fa + '.gz'])
        n_snv = None   # display differences replaced by the formal claim model (website/tour_models.py)
        refinfo[key] = dict(key=key, D=D, R=R, donor_contig=contig(D), recipient_contig=contig(R), donor_len=len(Dseq),
                            recipient_len=len(Rseq), n_diff_sites=n_snv, chm13_len_D=lens[D[:-1]], chm13_len_R=lens[R[:-1]])
        for name, s in ((contig(D), Dseq), (contig(R), Rseq)):
            write_fa(os.path.join(CACHE, f'{key}.{name}.fa'), [(name, s)])
        print('ref', key, len(Dseq), len(Rseq), 'diff sites', n_snv, flush=True)

    def align_parts(c, mols, out, key):
        ri = refinfo[key]; tD = os.path.join(CACHE, f'{key}.{ri["donor_contig"]}.fa'); tR = os.path.join(CACHE, f'{key}.{ri["recipient_contig"]}.fa')
        sams = []
        for preset, sel in (('map-ont', lambda m: not str(m['platform']).startswith('HiFi')), ('map-hifi', lambda m: str(m['platform']).startswith('HiFi'))):
            sub = [m for m in mols if sel(m)]
            if not sub: continue
            q1, q2 = os.path.join(CACHE, 'q1.fa'), os.path.join(CACHE, 'q2.fa')
            write_fa(q1, [(f'{m["read_id"]}/part1', m['seq'][:m['switch_lo']]) for m in sub])
            write_fa(q2, [(f'{m["read_id"]}/part2', m['seq'][m['switch_hi']:]) for m in sub])
            sams += [minimap(tD, q1, preset), minimap(tR, q2, preset)]
        ds = {m['read_id']: m['dataset'] for m in mols}
        hdr = header([(ri['donor_contig'], ri['donor_len']), (ri['recipient_contig'], ri['recipient_len'])], sorted(set(ds.values())),
                     f'{c["claim_id"]}: definitive molecules cut at the switch; part1 -> donor end, part2 -> recipient end; tip frame')
        return sam_to_bam(sams, hdr, out, rg_of=lambda n: ds[n.rsplit('/', 1)[0]])

    # 3. per-claim files
    for c in claims:
        key = f'{c["D"]}__{c["R"]}'; od = os.path.join(OUT, 'claims', c['slug']); os.makedirs(od, exist_ok=True)
        c['ref'] = key
        c['parts_alignments'] = align_parts(c, c['mols'], os.path.join(od, 'parts.bam'), key)
        # own assembly: filter the release own.bam to the shown molecules
        keep = {m['read_id'] for m in c['mols']}
        src = os.path.join(c['release_dir'], f'{c["slug"]}.own.bam')
        unsorted = os.path.join(CACHE, 'tmp', f'{c["slug"]}.own.unsorted.bam')
        with pysam.AlignmentFile(src) as fi, pysam.AlignmentFile(unsorted, 'wb', template=fi) as fo:
            for r in fi.fetch(until_eof=True):
                if r.query_name in keep: fo.write(r)
        run(['samtools', 'sort', '-o', os.path.join(od, 'own.bam'), unsorted], capture_output=True)
        run(['samtools', 'index', os.path.join(od, 'own.bam')])
        ownfa = os.path.join(CACHE, 'own', f'{c["slug"]}.fa'); os.makedirs(os.path.dirname(ownfa), exist_ok=True)
        shutil.copyfile(os.path.join(c['release_dir'], f'{c["slug"]}.own_windows.fa'), ownfa)
        meta = {}
        for l in open(ownfa):
            if l.startswith('>'):
                f = l[1:].split(); meta[f[0]] = dict(kv.split('=', 1) for kv in f[1:] if '=' in kv)
        with open(os.path.join(od, 'own_blocks.bed'), 'w') as fh:
            fh.write('track name="catalogue blocks" itemRgb="On"\n')
            for w, mt in meta.items():
                eid = mt.get('end_id'); e = E[E.end_id == eid]
                tl = int(e.telo_len.iloc[0]) if len(e) and e.telo_len.iloc[0] > 0 else 0
                if tl: fh.write(f'{w}\t0\t{tl}\ttelomere\t0\t.\t0\t{tl}\t40,40,40\n')
                for bl in MB[MB.end_id == eid].itertuples():
                    if pd.isna(bl.tip_start): continue
                    s0, e0 = int(bl.tip_start), int(min(bl.tip_end, 100_000))
                    lab = f'home ({bl.source_arm})' if bl.is_home else f'{bl.best_end}-like ({bl.label})'
                    col = '0,114,178' if bl.is_home else ('213,94,0' if bl.best_end == c['D'] else '150,150,150')
                    fh.write(f'{w}\t{s0}\t{e0}\t{lab}\t0\t.\t{s0}\t{e0}\t{col}\n')
        with pysam.AlignmentFile(os.path.join(od, 'own.bam')) as b:
            own = [dict(name=r.query_name, ref=r.reference_name, start=r.reference_start + 1, end=r.reference_end,
                        mapq=r.mapping_quality, nm=r.get_tag('NM') if r.has_tag('NM') else None, aligned=r.query_alignment_length)
                   for r in b.fetch(until_eof=True) if not r.is_secondary and not r.is_supplementary]
        c['own_alignments'] = own; c['windows'] = meta
        c['claimed_window'] = next((w for w, m in meta.items() if m.get('claimed') == 'yes'), list(meta)[0])
        print('claim', c['slug'], len(c['mols']), 'molecules', flush=True)
    # 4. chapter extras
    ex = os.path.join(OUT, 'extra'); os.makedirs(ex, exist_ok=True)
    hero = next(c for c in claims if c['slug'] == HERO_SLUG); m = next(m for m in hero['mols'] if m['read_id'] == HERO_READ)
    ri = refinfo[hero['ref']]; q = os.path.join(CACHE, 'hero.fa'); write_fa(q, [(HERO_READ, m['seq'])])
    hero_whole = {}
    for which in ('donor', 'recipient'):
        nm = ri[f'{which}_contig']; t = os.path.join(CACHE, f'{hero["ref"]}.{nm}.fa')
        hero_whole[which] = sam_to_bam([minimap(t, q, 'map-ont')], header([(ri['donor_contig'], ri['donor_len']),
                                       (ri['recipient_contig'], ri['recipient_len'])]), os.path.join(ex, f'hero_whole_{which}.bam'))
    heropart = os.path.join(ex, 'hero_parts.bam')
    hero_parts = align_parts(hero, [m], heropart, hero['ref'])
    tech = next(c for c in claims if c['slug'] == TECH_SLUG); tech_files = {}
    for mm in tech['mols']:
        f = os.path.join(ex, f'tech_{mm["dataset"]}.bam'); align_parts(tech, [mm], f, tech['ref']); tech_files[mm['dataset']] = os.path.basename(f)
    for c in claims:
        c.pop('release_dir')
        for mm in c['mols']: mm.pop('seq')
    json.dump(dict(refs=refinfo, claims=claims, hero=dict(slug=HERO_SLUG, read=HERO_READ, whole=hero_whole, parts=hero_parts),
                   tech=dict(slug=TECH_SLUG, files=tech_files)),
              open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)
    tot = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(OUT) for f in fs)
    big = max(((os.path.getsize(os.path.join(dp, f)), os.path.join(dp, f)) for dp, _, fs in os.walk(OUT) for f in fs))
    print(f'tour data {tot / 1e6:.1f} MB; largest file {big[0] / 1e6:.1f} MB {os.path.relpath(big[1], OUT)}')


def build_combo():
    """Combined view per claim: one reference holding the person's own end windows (hap1/hap2, tip frame) and the CHM13
    donor and recipient ends, so molecules -> own assembly, own assembly -> CHM13 and molecules -> CHM13 show together.
      claims/<slug>/combo.fa.gz (+.fai, .gzi)  own windows + CHM13_<D>_end + CHM13_<R>_end
      claims/<slug>/asm.bam (+.bai)            each own window aligned to the two CHM13 ends (minimap2 asm20; primary and
                                               supplementary pieces >= 2 kb), records named after the window
    Also keeps only the primary alignment in extra/hero_whole_*.bam (aligning a whole molecule to one end alone makes
    spurious short supplementary hits)."""
    J = json.load(open(os.path.join(OUT, 'tour_data.json')))
    for c in J['claims']:
        od = os.path.join(OUT, 'claims', c['slug']); ri = J['refs'][c['ref']]
        own = [(r.name, r.sequence) for r in pysam.FastxFile(os.path.join(CACHE, 'own', f'{c["slug"]}.fa'))]
        ends = [(ri[k], next(pysam.FastxFile(os.path.join(CACHE, f'{c["ref"]}.{ri[k]}.fa'))).sequence)
                for k in ('donor_contig', 'recipient_contig')]
        fa = os.path.join(od, 'combo.fa'); write_fa(fa, own + ends)
        run(['bgzip', '-f', '-i', fa]); run(['samtools', 'faidx', fa + '.gz'])
        # each catalogue block of each window, placed on the CHM13 end its label names (home -> recipient end, a block
        # whose best CHM13 end is the donor -> donor end; blocks labelled with any other end are not placed here)
        ownd = dict(own); targets = {ri['donor_contig']: ends[0][1], ri['recipient_contig']: ends[1][1]}
        bl = [l.rstrip('\n').split('\t') for l in open(os.path.join(od, 'own_blocks.bed')) if not l.startswith('track')]
        tf = {}
        for nm, s in ends:
            tf[nm] = os.path.join(CACHE, f'{c["ref"]}.{nm}.fa')
        pieces, sams = [], []
        for w, s0, e0, lab, *_ in bl:
            s0, e0 = int(s0), int(e0)
            if lab == 'telomere' or e0 - s0 < 2000: continue
            if lab.startswith('home'): tgt, other = ri['recipient_contig'], ri['donor_contig']
            elif lab.startswith(c['D'] + '-like'): tgt, other = ri['donor_contig'], ri['recipient_contig']
            else:
                pieces.append(dict(window=w, win_start=s0 + 1, win_end=e0, label=lab, placed=None)); continue
            q = os.path.join(CACHE, 'tmp', 'block.fa'); nm = f'{w}:{s0 + 1}-{e0}'
            write_fa(q, [(nm, ownd[w][s0:e0])])
            sam = minimap(tf[tgt], q, 'asm10'); sams.append(sam)
            ident = {}
            for which in (tgt, other):
                pf = run(['minimap2', '-c', '-x', 'asm10', '--secondary=no', tf[which], q], capture_output=True, text=True).stdout
                best = max((l.split('\t') for l in pf.splitlines()), key=lambda f: int(f[10]), default=None)
                ident[which] = (round(int(best[9]) / int(best[10]), 4), int(best[3]) - int(best[2])) if best else (None, 0)
            pieces.append(dict(window=w, win_start=s0 + 1, win_end=e0, label=lab, placed=tgt,
                               ident_placed=ident[tgt][0], aligned_placed=ident[tgt][1],
                               other=other, ident_other=ident[other][0], aligned_other=ident[other][1]))
        hdr = header([(n, len(s)) for n, s in own + ends], [n for n, _ in own],
                     f'{c["claim_id"]}: catalogue blocks of the person\'s own end windows placed on the CHM13 end their label '
                     f'names (minimap2 asm10)')
        sam_to_bam(sams, hdr, os.path.join(CACHE, 'tmp', 'asm_keep.bam'), rg_of=lambda n: n.rsplit(':', 1)[0])
        run(['samtools', 'sort', '-o', os.path.join(od, 'asm.bam'), os.path.join(CACHE, 'tmp', 'asm_keep.bam')], capture_output=True)
        run(['samtools', 'index', os.path.join(od, 'asm.bam')])
        c['combo'] = dict(contigs=[[n, len(s)] for n, s in own + ends], asm_pieces=pieces)
        print('combo', c['slug'], len(pieces), 'pieces', flush=True)
    for which in ('donor', 'recipient'):
        f = os.path.join(OUT, 'extra', f'hero_whole_{which}.bam'); tmp = os.path.join(CACHE, 'tmp', f'hero_{which}.bam')
        with pysam.AlignmentFile(f) as fi, pysam.AlignmentFile(tmp, 'wb', template=fi) as fo:
            for r in fi.fetch(until_eof=True):
                if not r.is_supplementary and not r.is_secondary: fo.write(r)
        shutil.copyfile(tmp, f); run(['samtools', 'index', f])
    json.dump(J, open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)


def model_agreement():
    """Per molecule part, against the formal claim model (claims/<slug>/model.json, website/tour_models.py):
    part1 on the donor end: letters at the diagnostic SNV sites (donor letter / comparator letter / other) and bp of donor-only
    sequence covered; part2 on the recipient end (only when the proximal comparator P is the recipient R): letters at the
    sites (comparator letter / donor letter / other) and bp of proximal-only sequence covered. Primary + supplementary pieces."""
    J = json.load(open(os.path.join(OUT, 'tour_data.json')))
    for c in J['claims']:
        ri = J['refs'][c['ref']]; m = json.load(open(os.path.join(OUT, 'claims', c['slug'], 'model.json')))
        pr = m['P'] == m['R']; c['model'] = dict(D=m['D'], P=m['P'], R=m['R'], n_sites=len(m['sites']), p_is_r=pr,
                                                 donor_only_bp=sum(b - a for a, b in m['private']['donor']),
                                                 proximal_only_bp=sum(b - a for a, b in m['private']['proximal']))
        dsite = {s['d_pos']: (s['d_letter'], s['p_letter']) for s in m['sites']}
        psite = {s['p_pos']: (s['p_letter'], s['d_letter']) for s in m['sites']}
        calls = {}
        with pysam.AlignmentFile(os.path.join(OUT, 'claims', c['slug'], 'parts.bam')) as b:
            for r in b.fetch(until_eof=True):
                if r.is_secondary: continue
                rid, part = r.query_name.rsplit('/', 1)
                k = calls.setdefault(rid, {}).setdefault(part, dict(own=0, other=0, n=0, private_bp=0))
                on_d = r.reference_name == ri['donor_contig']
                if part == 'part1' and not on_d: continue
                if part == 'part2' and (on_d or not pr): continue
                sm = dsite if on_d else psite; priv = m['private']['donor' if on_d else 'proximal']
                seq = r.query_sequence
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    if rp in sm:
                        own, alt = sm[rp]; x = seq[qp].upper(); k['n'] += 1; k['own'] += x == own; k['other'] += x == alt
                k['private_bp'] += sum(max(0, min(r.reference_end, b_) - max(r.reference_start, a_)) for a_, b_ in priv)
        c['site_calls'] = calls
        tot = {p_: dict(own=sum(v.get(p_, {}).get('own', 0) for v in calls.values()),
                        other=sum(v.get(p_, {}).get('other', 0) for v in calls.values()),
                        n=sum(v.get(p_, {}).get('n', 0) for v in calls.values()),
                        private_bp=sum(v.get(p_, {}).get('private_bp', 0) for v in calls.values())) for p_ in ('part1', 'part2')}
        c['site_totals'] = tot
    json.dump(J, open(os.path.join(OUT, 'tour_data.json'), 'w'), indent=1, default=str)


def clean_headers():
    """Drop @PG header lines (they record local command paths) from every BAM under OUT, then re-index."""
    n = 0
    for root, _, fs in os.walk(OUT):
        for x in fs:
            if not x.endswith('.bam'): continue
            f = os.path.join(root, x)
            h = run(['samtools', 'view', '--no-PG', '-H', f], capture_output=True, text=True).stdout
            if '@PG' not in h: continue
            hf = os.path.join(CACHE, 'tmp', 'hdr.sam'); tmp = os.path.join(CACHE, 'tmp', 'reheader.bam')
            os.makedirs(os.path.dirname(hf), exist_ok=True)
            with open(hf, 'w') as fh: fh.write(''.join(l + '\n' for l in h.splitlines() if not l.startswith('@PG')))
            with open(tmp, 'wb') as fh: run(['samtools', 'reheader', '--no-PG', hf, f], stdout=fh)
            shutil.copyfile(tmp, f); run(['samtools', 'index', f]); n += 1
    print('clean_headers: rewrote', n, 'BAM headers', flush=True)


if __name__ == '__main__':
    import sys
    if '--sites-only' in sys.argv:
        model_agreement(); clean_headers(); sys.exit()
    if '--combo-only' not in sys.argv: main()
    build_combo(); model_agreement(); clean_headers()   # then: python website/tour_models.py; python website/tour_build.py --sites-only
