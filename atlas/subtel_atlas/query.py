"""Query engine: telomeric ends of contigs / windows / reads -> arm, CHM13 mosaic, S5/S10/S20 states, structural allele, haplogroup.

Every step reuses the catalogue functions with the catalogue parameters: telomere.find_block (boundary), the catalogue window rule
(window_length), tracks.* (one minimap2 pass of window[s0:] per end, per-base identity tracks, paralog-family labels, S5/S10/S20),
alleles.core/match_allele (structures and allele matching), arms.* (name / probe / bin-based mosaic rule), haplo (sketch ANI to
representative medoids + alignment check), and the boundary-shift test of pipeline/catalogue_08e_boundary_stable.py. An end that is
in the catalogue is therefore re-described exactly as in data/catalogue/ends.tsv when the same boundary and arm are used and the
minimap2 backend is chosen (tests/validation/VALIDATION.md).
"""
import os, re, io, json, time, gzip, collections
from types import SimpleNamespace
import numpy as np, pandas as pd
from . import __version__, CATALOGUE_VERSION
from .common import WINDOW, rc, group_of, CHM13_ENDS
from .telomere import find_block, CAN
from . import bins as B, mosaic as M, alleles as AL, arms as AR, haplo as H, refs as R, sketch as S, tracks as TK

MIN_TELO = 500                 # an end is telomeric if its clean (CCCTAA)n block is >= 500 bp (Release 2 tip definition)
MAXJUNK = 8000
HG_ALN_IDENT_MIN, HG_ALN_FRAC_MIN = 0.99, 0.90   # ~1st percentile of catalogue member-to-medoid alignments (identity 0.9906, fraction 0.923)
HG_CLOSE = 0.001
SUB_MIN = 80_000               # catalogue rule (pipeline/catalogue_03): boundary > WINDOW - SUB_MIN -> window extended to keep 80 kb
SHIFT = 500                    # boundary-stability test (as pipeline/catalogue_08e_boundary_stable.py)
_SEED = [m * 3 for m in CAN]

COLUMNS = ['query_id', 'input_name', 'mode', 'side', 'contig_len', 'win_start', 'win_end', 'strand', 'window_len', 'telomere_found',
           'window_extended', 'junk_len', 'telo_len', 'telomere_source', 's0', 'sub_len', 'full_length',
           'home_arm', 'catalogue_arm', 'arm_source', 'arm_confidence', 'probe_arm', 'probe_confidence', 'mosaic_arm', 'mosaic_arm_score_frac',
           'profile', 'structure', 'switches_sub', 'switch_res', 'core_structure', 'core_switches_sub', 'n_core_blocks', 'nativity_class',
           'S5_state', 'S5_family', 'S5_best_end', 'S10_state', 'S10_family', 'S10_best_end', 'S20_state', 'S20_family', 'S20_best_end',
           'S20_best_ident', 'form', 'boundary_stable', 'boundary_stable_S20', 'boundary_stable_allele', 'S20_state_minus500',
           'S20_state_plus500', 'allele_minus500', 'allele_plus500',
           'structural_allele_id', 'sa_is_native', 'sa_nativity_class', 'sa_assignment', 'sa_candidates', 'sa_max_switch_delta', 'sa_nearest_same_structure', 'sa_nearest_delta',
           'sa_tier', 'sa_evidence', 'sa_n_samples', 'sa_r2_freq', 'sa_representative_end_id',
           'haplogroup_id', 'hg_status', 'hg_ani', 'hg_second_id', 'hg_second_ani', 'hg_compared_len', 'hg_aln_ident', 'hg_aln_frac',
           'hg_tier', 'hg_n_samples', 'hg_r2_freq', 'novelty_flags', 'qc_flags']


def has_seed(s, maxjunk=MAXJUNK):
    """Necessary condition for find_block(): three clean CCCTAA units starting within maxjunk bases (fast pre-screen)."""
    head = s[:maxjunk + 17].upper()
    return any(t in head for t in _SEED)


def call_telomere(tip):
    tb = find_block(tip) if has_seed(tip) else None
    return tb


# ------------------------------------------------------------------------------------------------ reference
class Reference:
    """Reference directory (see refs.py): CHM13 ends, structural alleles, haplogroups and their representative medoid windows."""
    def __init__(self, refdir=None, verify=True, catalogue_dir=None):
        self.dir = refdir or R.default_dir()
        if verify:
            ok = R.verify(self.dir, quiet=True)
            if not all(ok.values()):
                bad = ', '.join(f for f, v in ok.items() if not v)
                raise SystemExit(f'reference directory {self.dir} is missing or fails checksums ({bad}); run: subtel-atlas fetch-refs --dir {self.dir}')
        self.chm13_fa = os.path.join(self.dir, R.CHM13_FA)
        if not os.path.exists(self.chm13_fa) and os.path.exists(self.chm13_fa + '.gz'): self.chm13_fa += '.gz'   # test fixtures
        self.refends = M.RefEnds(R.readfa(self.chm13_fa))      # (kept for the bin-based arm rule)
        cat = self.dir if os.path.exists(os.path.join(self.dir, 'alleles_structural.tsv')) else (catalogue_dir or R.repo_catalogue_dir())
        rd = lambda f: pd.read_csv(os.path.join(cat, f), sep='\t', low_memory=False)
        self.SA = rd('alleles_structural.tsv').set_index('allele_id'); self.HG = rd('haplogroups.tsv').set_index('haplogroup_id')
        self.families = TK.load_families(os.path.join(cat, 'paralog_families.tsv'))
        fr = rd('allele_frequencies_by_superpop.tsv'); self.freq = {}
        for r in fr.itertuples(): self.freq.setdefault(r.allele_id, {})[r.superpopulation] = r.frequency
        self.alleles = {}
        for aid, r in self.SA.sort_index().iterrows():
            sw = [int(x) for x in str(r.switches_median).split(',')] if isinstance(r.switches_median, str) and r.switches_median else []
            self.alleles.setdefault(r.arm, []).append(dict(allele_id=aid, structure=r.structure, sw=sw, n=int(r.n_defining)))
        self.rep = {}                                   # end_id -> (window, s0)
        with gzip.open(os.path.join(cat, 'representatives.fa.gz'), 'rt') as fh:
            hdr, buf = None, []
            for line in fh:
                if line[0] == '>':
                    if hdr: self._add_rep(hdr, ''.join(buf))
                    hdr, buf = line[1:].strip(), []
                else: buf.append(line.strip())
            if hdr: self._add_rep(hdr, ''.join(buf))
        self.medoids = {}                               # arm -> [(haplogroup_id, medoid subtelomere)]
        for hid, r in self.HG.sort_index().iterrows():
            if bool(r.representative_in_fasta) and r.medoid_end_id in self.rep:
                w, s0 = self.rep[r.medoid_end_id]; self.medoids.setdefault(r.arm, []).append((hid, w[s0:s0 + H.HG_LEN], r.medoid_end_id))
        self._sk = {}

    def _add_rep(self, hdr, seq):
        kv = dict(t.split('=', 1) for t in hdr.split()[1:] if '=' in t)
        self.rep[kv['end_id']] = (seq, int(kv['subtelomere_start']))

    def medoid_sketch(self, key, ms):
        if key not in self._sk: self._sk[key] = S.sketch(ms, H.K, H.SCALE)
        return self._sk[key]

    def nearest_haplogroups(self, arm, sub):
        """[(ANI, haplogroup_id, medoid subtelomere)] over the represented haplogroups of a catalogue arm, best first; compared over
        the query's subtelomere length (<= 80 kb), as for partial ends in the catalogue."""
        Lp = len(sub); ske = S.sketch(sub, H.K, H.SCALE); out = []
        for hid, ms, mend in self.medoids.get(arm, []):
            skm = self.medoid_sketch(mend, ms) if Lp >= len(ms) else S.sketch(ms[:Lp], H.K, H.SCALE)
            out.append((H.ani(ske, skm), hid, ms))
        return sorted(out, key=lambda x: (-x[0], x[1]))


# ------------------------------------------------------------------------------------------------ inputs
def _end(end_id, name, mode, side, contig_len, ws, we, strand, tip, arm_hint=None, arm_hint_source='', probes=None):
    tb = call_telomere(tip)
    return dict(end_id=end_id, input_name=name, mode=mode, side=side, contig_len=contig_len, win_start=ws, win_end=we, strand=strand,
                seq=tip, telomere_found=tb is not None, junk_len=tb['junk_len'] if tb else 0, telo_len=tb['telo_len'] if tb else 0,
                arm_hint=arm_hint, arm_hint_source=arm_hint_source, probes=probes or [], boundary_source='catalogue_caller',
                window_extended=False)


def window_length(s0, available, tel_to_end=False):
    """Catalogue window rule: 100 kb, or ceil((s0 + 80 kb) / 5 kb) x 5 kb when the boundary lies beyond 20 kb (so that a full-length end
    keeps >= 80 kb of subtelomere), never more than the available sequence."""
    if s0 > WINDOW - SUB_MIN and available > WINDOW and not tel_to_end:
        return min(available, -(-(s0 + SUB_MIN) // 5000) * 5000)
    return min(WINDOW, available)


def _extend(e, full_tip):
    """apply window_length() to an end whose telomere was called on its first 100 kb; full_tip = all available tip-frame sequence"""
    s0 = e['junk_len'] + e['telo_len']; w = window_length(s0, len(full_tip), s0 >= len(e['seq']) - 100)
    e['window_extended'] = w > len(e['seq'])
    if e['window_extended']: e['seq'] = full_tip[:w]
    return e


def contig_ends(name, seq, min_telomere=MIN_TELO, window=WINDOW, name_arms=True, arm_map=None, want_probes=False):
    """Telomeric ends of one contig: the outer `window` bp of each side in tip frame (R ends reverse-complemented) whose telomere
    caller block is >= min_telomere. 1-based inclusive contig coordinates as in ends.tsv."""
    out = []; L = len(seq)
    for side in ('L', 'R'):
        w = min(window, L)
        raw = seq[:w] if side == 'L' else seq[L - w:]
        tip = raw.upper() if side == 'L' else rc(raw.upper())
        if not has_seed(tip): continue
        e = _end(f'{name}#{side}', name, 'contig', side, L, 1 if side == 'L' else L - w + 1, w if side == 'L' else L,
                 '+' if side == 'L' else '-', tip)
        if not e['telomere_found'] or e['telo_len'] < min_telomere: continue
        s0 = e['junk_len'] + e['telo_len']; w2 = window_length(s0, L, s0 >= len(tip) - 100)
        if w2 > w:                                             # long telomere: extended window, as in the catalogue
            raw = seq[:w2] if side == 'L' else seq[L - w2:]
            e.update(seq=raw.upper() if side == 'L' else rc(raw.upper()), win_start=1 if side == 'L' else L - w2 + 1,
                     win_end=w2 if side == 'L' else L, window_extended=True)
        hint, src = (arm_map or {}).get(e['end_id']) or (arm_map or {}).get(name + ':' + side), 'arm_map'
        if not hint and name_arms: hint, src = AR.arm_from_name(name, side), 'contig_name'
        e['arm_hint'], e['arm_hint_source'] = (hint, src) if hint else (None, '')
        if want_probes: e['probes'] = [(off, seq[a:b].upper()) for off, (a, b) in AR.probe_regions(L, side)]
        out.append(e)
    return out


def window_end(name, comment, seq, name_arms=True, arm_map=None):
    """A pre-cut tip-frame window (telomere first). Header tokens arm=<arm> (e.g. arm=20q), side=L|R and junk_len=<bp> telo_len=<bp>
    (a telomere boundary from another caller; otherwise the catalogue caller is run on the window) are honoured."""
    kv = dict(t.split('=', 1) for t in (comment or '').split() if '=' in t)
    side = kv.get('side') or (name[-1] if re.search(r'[#:_|][LR]$', name) else '')
    S_ = seq.upper(); e = _end(name, name, 'window', side, '', '', '', '', S_[:WINDOW])
    if 'junk_len' in kv and 'telo_len' in kv:          # caller-supplied telomere boundary (overrides the caller)
        e.update(junk_len=int(kv['junk_len']), telo_len=int(kv['telo_len']), telomere_found=int(kv['telo_len']) > 0, boundary_source='header')
    _extend(e, S_)
    hint, src = kv.get('arm'), 'header'
    if not hint and arm_map: hint, src = arm_map.get(name), 'arm_map'
    if not hint and name_arms and side: hint, src = AR.arm_from_name(name, side), 'contig_name'
    e['arm_hint'], e['arm_hint_source'] = (hint, src) if hint else (None, '')
    return e


def read_end(name, seq, min_telomere=MIN_TELO):
    """A long read: telomere at the read start (C-strand, CCCTAA) or at its end (G-strand; the read is reverse-complemented)."""
    s = seq.upper(); best = None
    for strand, x in (('+', s), ('-', rc(s))):
        if not has_seed(x): continue
        tb = find_block(x[:WINDOW])
        if tb and tb['telo_len'] >= min_telomere and (best is None or tb['telo_len'] > best[1]['telo_len']): best = (strand, tb, x)
    if best is None: return None
    e = _end(name, name, 'read', '', len(s), '', '', best[0], best[2][:WINDOW])
    return _extend(e, best[2])


def iter_records(path):
    import mappy
    for name, seq, qual, comment in mappy.fastx_read(path, read_comment=True):
        yield name, comment, seq


def iter_bam(path):
    try:
        import pysam
    except ImportError:
        raise SystemExit('BAM/CRAM input needs pysam (pip install pysam); or convert to FASTQ (samtools fastq)')
    with pysam.AlignmentFile(path, 'rb', check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or not r.query_sequence: continue
            s = r.query_sequence
            yield r.query_name, None, (rc(s) if r.is_reverse else s)   # original read orientation


def load_arm_map(path):
    if not path: return None
    m = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2 and not line.startswith('#'): m[f[0]] = f[1]
    return m


# ------------------------------------------------------------------------------------------------ analysis
def _nan(x): return x is None or (isinstance(x, float) and x != x)


def _probe_rows(e, best, clen):
    rows = [AR.probe_row(e['end_id'], off, e['side'], hit, clen) for off, hit in sorted(best.items())]
    P = pd.DataFrame(rows, columns=['end_id', 'offset', 'arm', 'mapq', 'aln', 'consistent', 'dist_from_end', 'informative', 'score'])
    parm, conf = AR.probe_consensus(P[P.informative] if len(P) else P)
    return parm, conf, rows


def _probe_call(e, genome, clen):
    """Probe rule with a mappy asm5 aligner of the CHM13 genome (fallback when minimap2 is not available)."""
    if genome is None or not e['probes']: return '', '', []
    best = {}
    for off, ps in e['probes']:
        for h in genome.map(ps):
            if not h.is_primary: continue
            clen[h.ctg] = h.ctg_len
            AR.best_probe_hit(best, off, h.q_en - h.q_st, h.mapq, h.ctg, h.r_st, h.r_en, '+' if h.strand == 1 else '-')
    return _probe_rows(e, best, clen)


def probe_calls_minimap2(ends, genome_path, threads, exe):
    """Probe rule exactly as the catalogue: minimap2 -x asm5 --secondary=no against the CHM13 genome (FASTA or .mmi), all probes in
    one run. Returns {end_id: (probe arm, confidence, rows)}."""
    import subprocess, tempfile
    q = [(f"{e['end_id']}@{off}", ps) for e in ends for off, ps in e['probes']]
    if not q: return {}
    best, clen = collections.defaultdict(dict), {}
    with tempfile.TemporaryDirectory(prefix='subtel_atlas_') as td:
        fa = os.path.join(td, 'p.fa'); open(fa, 'w').write(''.join(f'>{n}\n{s}\n' for n, s in q))
        p = subprocess.run([exe, '-x', 'asm5', '--secondary=no', '-t', str(threads), genome_path, fa], capture_output=True, text=True)
        if p.returncode: raise RuntimeError('minimap2 (probes) failed: ' + p.stderr[-2000:])
        for line in p.stdout.splitlines():
            f = line.split('\t'); e, off = f[0].rsplit('@', 1); clen[f[5]] = int(f[6])
            AR.best_probe_hit(best[e], int(off), int(f[3]) - int(f[2]), int(f[11]), f[5], int(f[7]), int(f[8]), f[4])
    return {e['end_id']: _probe_rows(e, best.get(e['end_id'], {}), clen) for e in ends if e['probes']}


_CTX = {}
def _finish(i):
    """Mosaic, structural allele and haplogroup of job i (runs in forked worker processes, which inherit _CTX)."""
    e, alns, home, src, conf = _CTX['jobs'][i]; ref = _CTX['ref']
    s0 = e['junk_len'] + e['telo_len']; seq = e['seq']; grp = group_of(home)
    L = len(seq); tr = TK.tracks(max(0, L - s0), alns)
    if L - s0 > 0:
        blocks, rec = TK.mosaic_from_tracks(e['end_id'], tr, s0, L, grp, ref.families, CHM13_ENDS)
    else:
        blocks, rec = [], dict(profile='', structure='', switches_sub='', switch_res='')
    sub_len = L - s0; full = L >= WINDOW and sub_len >= AL.FULL_SUB
    if blocks:
        cbl, struct, sw = AL.core(TK.core_input(blocks))
    else: cbl, struct, sw = [], '', []
    out = dict(rec=rec, blocks=blocks, core=[dict(core_block=i + 1, source=b['hsrc'] if b['source'] == 'home' else b['source'],
                                                  is_home=b['source'] == 'home', sub_start=int(round(b['s'])), sub_end=int(round(b['e'])))
                                             for i, b in enumerate(cbl)],
               struct=struct, sw=sw, sub_len=sub_len, full=full, grp=grp, s0=s0)
    out.update(assign(ref, grp, struct, sw, sub_len, full, seq[s0:s0 + H.HG_LEN]))
    out['shift'] = shift_calls(ref, e['end_id'], tr, s0, L, grp, out) if L - s0 > SHIFT else None
    return out


def _calls_at(ref, e, tr, s0, L, grp, d):
    """core structure, allele, S20 state and family with the boundary moved by d (identity tracks shifted, no re-alignment);
    same computation as pipeline/catalogue_08e_boundary_stable.calls"""
    bl, rec = TK.mosaic_from_tracks(e, tr, s0 + d, L, grp, ref.families, CHM13_ENDS)
    _, st, sw = AL.core(TK.core_input(bl)) if bl else (None, '', [])
    sub_len = L - s0 - d; full = L >= WINDOW and sub_len >= AL.FULL_SUB - max(0, d)
    best, _ = AL.match_allele(st, sw, sub_len, full, ref.alleles.get(grp, []))
    return st, best['allele_id'] if best else '', rec['S20_state'], rec['S20_family']


def shift_calls(ref, e, tr, s0, L, grp, out):
    """boundary stability of the S20 call and the structural allele under a -500 / +500 bp boundary shift"""
    r0 = (out['struct'], out['sa']['allele_id'] if out['sa'] else '', out['rec'].get('S20_state', 'none'), out['rec'].get('S20_family', '') or '')
    r = {d: _calls_at(ref, e, TK.shift_tracks(tr, d), s0, L, grp, d) for d in (-SHIFT, SHIFT)}
    sa = all(r[d][1] == r0[1] for d in r); s20 = all(r[d][2] == r0[2] and (r[d][3] or '') == r0[3] for d in r)
    return dict(boundary_stable=sa and s20, boundary_stable_S20=s20, boundary_stable_allele=sa, S20_state_minus500=r[-SHIFT][2],
                S20_state_plus500=r[SHIFT][2], allele_minus500=r[-SHIFT][1], allele_plus500=r[SHIFT][1])


def assign(ref, grp, struct, sw, sub_len, full, sub):
    """Structural allele and haplogroup of an end, given its catalogue arm, core structure/switches, subtelomere length, full-length
    status and compared subtelomere (seq[s0:s0 + 80 kb]). ref may carry modified .alleles / .medoids (leave-one-out validation)."""
    al = ref.alleles.get(grp, [])
    best, cand = AL.match_allele(struct, sw, sub_len, full, al) if struct else (None, [])
    same = sorted((c for c in al if c['structure'] == struct and len(c['sw']) == len(sw)),
                  key=lambda c: (AL.max_switch_delta(c['sw'], sw), -c['n'], c['allele_id']))
    out = dict(sa=best, sa_cand=cand, sa_near=same[0] if same else None, sa_struct_known=any(c['structure'] == struct for c in al), hg=None)
    if len(sub) >= 5000 and ref.medoids.get(grp):
        nn = ref.nearest_haplogroups(grp, sub)
        a0, h0, ms = nn[0]; ident, frac = H.gc_ident(sub, ms)
        out['hg'] = dict(id=h0, ani=a0, second=(nn[1][1], nn[1][0]) if len(nn) > 1 else ('', np.nan), n_compared=len(nn),
                         compared_len=len(sub), aln_ident=ident, aln_frac=frac, top=[(h, round(a, 5)) for a, h, _ in nn[:3]])
    return out


def analyse(ends, ref, backend='auto', threads=1, genome=None, minimap2=None, log=None):
    """ends: list of end dicts (contig_ends / window_end / read_end). Returns a list of result records (input order)."""
    log = log or (lambda *a: None)
    if not ends: return []
    if backend == 'auto': backend = 'minimap2' if B.find_minimap2(minimap2) else 'mappy'
    exe = B.find_minimap2(minimap2) if backend == 'minimap2' else None
    if backend == 'minimap2' and not exe: raise SystemExit('minimap2 not found on PATH (use --aligner mappy or set SUBTEL_ATLAS_MINIMAP2)')
    t0 = time.time()
    # query names <end_id>@0 and the > 1 kb rule as in pipeline/catalogue_04t (minimap2 breaks score ties with a hash of the name)
    tq = [(e['end_id'] + '@0', e['seq'][e['junk_len'] + e['telo_len']:]) for e in ends if len(e['seq']) - e['junk_len'] - e['telo_len'] > 1000]
    TA = TK.align_minimap2(ref.chm13_fa, tq, threads, exe) if backend == 'minimap2' else TK.align_mappy(ref.chm13_fa, tq, threads)
    TA = {k[:-2]: v for k, v in TA.items()}
    log(f'  track alignment ({backend}): {len(tq)} subtelomeres, {time.time() - t0:.1f} s'); t0 = time.time()
    bt = B.bin_table([(e['end_id'], e['seq'], e['junk_len'] + e['telo_len']) for e in ends], ref.chm13_fa, backend, threads, exe)
    log(f'  bins aligned ({backend}): {len(bt)} queries, {time.time() - t0:.1f} s')
    fixm = bt.bin.isin(M.FIX)
    fixq = {t: {r.end_id: r for r in bt[bt.bin == t].itertuples()} for t in M.FIX}
    Bw = bt[~fixm].copy(); Bw['k'] = Bw.bin.astype(int)
    G = {e: g for e, g in Bw.groupby('end_id', sort=False)}
    clen = {}; jobs = []; armrec = []; pcalls = {}
    if genome is not None and isinstance(genome, str):                 # path: minimap2 executable (exact) or a mappy index
        if exe: pcalls = probe_calls_minimap2(ends, genome, threads, exe); genome = None
        else:
            import mappy
            genome = mappy.Aligner(genome, preset='asm5', n_threads=threads)
    for e in ends:
        s0 = e['junk_len'] + e['telo_len']; g = G.get(e['end_id'])
        ma = AR.mosaic_arm(g, s0) if g is not None else None
        parm, pconf, prow = pcalls[e['end_id']] if e['end_id'] in pcalls else _probe_call(e, genome, clen)
        home, src, conf = AR.decide(e['arm_hint'] or '', parm, pconf, ma[0] if ma else '')
        armrec.append(dict(home_arm=home, arm_source=src if src != 'contig_name' or not e['arm_hint_source'] else e['arm_hint_source'],
                           arm_confidence=conf, probe_arm=parm, probe_confidence=pconf, probe_detail=prow,
                           mosaic_arm=ma[0] if ma else '', mosaic_arm_score_frac=ma[1] if ma else np.nan))
        jobs.append((e, TA.get(e['end_id'], []), home, src, conf))
    t1 = time.time(); _CTX['ref'] = ref; _CTX['jobs'] = jobs
    if threads > 1 and len(jobs) > 1 and hasattr(os, 'fork'):
        import multiprocessing as mp
        with mp.get_context('fork').Pool(min(threads, len(jobs))) as pool: res = pool.map(_finish, range(len(jobs)), chunksize=1)
    else: res = [_finish(i) for i in range(len(jobs))]
    _CTX.clear()
    log(f'  mosaics, alleles and haplogroups: {len(jobs)} ends, {time.time() - t1:.1f} s')
    return [record(e, a, r, ref, backend) for e, a, r in zip(ends, armrec, res)]


def record(e, a, r, ref, backend):
    rec = r['rec']; grp = r['grp']; sa = r['sa']; hg = r['hg']; nov, qc = [], []
    o = dict(query_id=e['end_id'], input_name=e['input_name'], mode=e['mode'], side=e['side'], contig_len=e['contig_len'],
             win_start=e['win_start'], win_end=e['win_end'], strand=e['strand'], window_len=len(e['seq']), window_extended=bool(e.get('window_extended', False)),
             telomere_found=e['telomere_found'],
             junk_len=e['junk_len'], telo_len=e['telo_len'], telomere_source=e.get('boundary_source', ''), s0=r['s0'], sub_len=r['sub_len'], full_length=r['full'],
             home_arm=a['home_arm'], catalogue_arm=grp, arm_source=a['arm_source'], arm_confidence=a['arm_confidence'],
             probe_arm=a['probe_arm'], probe_confidence=a['probe_confidence'], mosaic_arm=a['mosaic_arm'],
             mosaic_arm_score_frac=a['mosaic_arm_score_frac'],
             profile=rec.get('profile', ''), structure=rec.get('structure', ''), switches_sub=rec.get('switches_sub', ''),
             switch_res=rec.get('switch_res', ''), core_structure=r['struct'], core_switches_sub=','.join(map(str, r['sw'])),
             n_core_blocks=len(r['core']), nativity_class=TK.nativity_class(r['struct'], grp, ref.families))
    for t in ('S5', 'S10', 'S20'):
        o[f'{t}_state'] = rec.get(f'{t}_state', 'none'); o[f'{t}_family'] = rec.get(f'{t}_family', '') or ''; o[f'{t}_best_end'] = rec.get(f'{t}_best_end', '')
    o['S20_best_ident'] = rec.get('S20_best_ident', np.nan)
    o['form'] = f"{a['home_arm']}>{o['S20_family']}" if o['S20_state'] == 'foreign' else ''
    sh = r.get('shift') or {}
    for k in ('boundary_stable', 'boundary_stable_S20', 'boundary_stable_allele', 'S20_state_minus500', 'S20_state_plus500',
              'allele_minus500', 'allele_plus500'): o[k] = sh.get(k, '')
    if sa is not None:
        S_ = ref.SA.loc[sa['allele_id']]
        o.update(structural_allele_id=sa['allele_id'],
                 sa_assignment='full_match' if r['full'] else ('partial_unique' if len(r['sa_cand']) == 1 else 'partial_ambiguous'),
                 sa_candidates=','.join(c['allele_id'] for c in r['sa_cand']), sa_max_switch_delta=AL.max_switch_delta(sa['sw'], r['sw']),
                 sa_is_native=bool(S_.is_native), sa_nativity_class=S_.get('nativity_class', '') if isinstance(S_.get('nativity_class', ''), str) else '',
                 sa_tier=S_.tier, sa_evidence=S_.evidence, sa_n_samples=int(S_.n_samples), sa_r2_freq=S_.r2_freq,
                 sa_representative_end_id=S_.representative_end_id)
    else:
        o.update(structural_allele_id='', sa_is_native='', sa_nativity_class='', sa_assignment='unassigned' if r['struct'] else 'not_assessed', sa_candidates='',
                 sa_max_switch_delta='', sa_tier='', sa_evidence='', sa_n_samples='', sa_r2_freq='', sa_representative_end_id='')
        if r['struct'] and grp != 'unassigned':
            if not r['full']: nov.append('sa_no_compatible_allele')
            elif r['sa_struct_known']: nov.append('sa_novel_switch_position')
            else: nov.append('sa_novel_structure')
    n_ = r['sa_near']
    o['sa_nearest_same_structure'] = n_['allele_id'] if n_ else ''
    o['sa_nearest_delta'] = AL.max_switch_delta(n_['sw'], r['sw']) if n_ else ''
    if hg is not None:
        Hr = ref.HG.loc[hg['id']]; ok = hg['ani'] >= H.HG_ANI
        o.update(haplogroup_id=hg['id'] if ok else '', hg_status='assigned' if ok else 'below_threshold', hg_ani=round(hg['ani'], 5),
                 hg_second_id=hg['second'][0], hg_second_ani=round(hg['second'][1], 5) if hg['second'][0] else np.nan,
                 hg_compared_len=hg['compared_len'], hg_aln_ident=round(hg['aln_ident'], 5) if not _nan(hg['aln_ident']) else np.nan,
                 hg_aln_frac=round(hg['aln_frac'], 4), hg_tier=Hr.tier if ok else '', hg_n_samples=int(Hr.n_samples) if ok else '',
                 hg_r2_freq=Hr.r2_freq if ok else '', hg_nearest_id=hg['id'], hg_top=hg['top'])
        if not ok: nov.append('hg_below_threshold')
        if ok and hg['second'][0] and hg['second'][1] >= H.HG_ANI and hg['ani'] - hg['second'][1] < HG_CLOSE: nov.append('hg_close_second')
        if _nan(hg['aln_ident']) or hg['aln_ident'] < HG_ALN_IDENT_MIN or hg['aln_frac'] < HG_ALN_FRAC_MIN: nov.append('hg_alignment_check_failed')
    else:
        o.update(haplogroup_id='', hg_status='not_assessed', hg_ani=np.nan, hg_second_id='', hg_second_ani=np.nan, hg_compared_len=0,
                 hg_aln_ident=np.nan, hg_aln_frac=np.nan, hg_tier='', hg_n_samples='', hg_r2_freq='', hg_nearest_id='', hg_top=[])
    if grp == 'unassigned': nov.append('arm_unassigned')
    if not e['telomere_found']: qc.append('telomere_not_found')
    elif e['telo_len'] < 1000: qc.append('short_telomere')
    if e['junk_len'] > 1000: qc.append('junk_gt_1kb')
    if not r['full']: qc.append('partial_window')
    if a['arm_source'] not in ('contig_name', 'header', 'arm_map'): qc.append('arm_inferred')
    if a['arm_confidence'] == 'low': qc.append('arm_low_confidence')
    if a['arm_confidence'] == 'conflict': qc.append('arm_conflict')
    if len(r['core']) >= 4: qc.append('complex_structure')
    if sa is not None and r['full'] and len(r['sa_cand']) > 1: qc.append('sa_multiple_candidates')
    if r.get('shift') and not r['shift']['boundary_stable']: qc.append('boundary_unstable')
    if e['mode'] == 'read': qc.append('experimental_read_mode')
    if backend == 'mappy': qc.append('mappy_backend_approximation')
    o['novelty_flags'] = ','.join(nov); o['qc_flags'] = ','.join(qc)
    o['blocks'] = [{k: b.get(k) for k in ('block', 'label', 'source', 'is_home', 'sub_start', 'sub_end', 'n_grid', 'ident', 'mean_tied',
                                          'chm13_dist_start', 'chm13_dist_end', 'next_switch_res')} for b in r['blocks']]
    o['core_blocks'] = r['core']; o['probe_detail'] = a['probe_detail']
    o['sa_frequency_by_superpop'] = ref.freq.get(sa['allele_id'], {}) if sa is not None else {}
    o['hg_frequency_by_superpop'] = ref.freq.get(o['haplogroup_id'], {}) if o['haplogroup_id'] else {}
    return o


# ------------------------------------------------------------------------------------------------ output
def _clean(v):
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, (np.floating, float)): return None if v != v else float(v)
    if isinstance(v, (np.bool_,)): return bool(v)
    if isinstance(v, dict): return {str(k): _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)): return [_clean(x) for x in v]
    return v


def write_tsv(records, fh, meta):
    tail = [meta['aligner'], meta['tool_version'], meta['catalogue_version']]
    fh.write('\t'.join(COLUMNS + ['aligner', 'tool_version', 'catalogue_version']) + '\n')
    for r in records:
        fh.write('\t'.join(['' if _nan(r.get(c)) else str(_clean(r.get(c))) for c in COLUMNS] + tail) + '\n')


def write_json(records, fh, meta):
    json.dump(_clean(dict(meta, ends=records)), fh, indent=1, sort_keys=True); fh.write('\n')


def meta(backend, args, n_input, n_ends):
    return dict(tool='subtel-atlas', tool_version=__version__, catalogue_version=CATALOGUE_VERSION, aligner=backend,
                parameters=dict(window_bp=WINDOW, bin_bp=B.BIN, min_telomere=args.get('min_telomere', MIN_TELO), near_tie=M.NEAR,
                                smoothing=M.SMOOTH, min_core_block_bp=AL.MIN_BLOCK, switch_tolerance_bp=AL.SW_TOL,
                                haplogroup_ani=H.HG_ANI, haplogroup_compared_bp=H.HG_LEN, minimap2_args=' '.join(B.MM2_ARGS),
                                track_minimap2_args=' '.join(TK.MM2_TRACK_ARGS), track_grid_bp=TK.GRID, track_window_bp=TK.WIN,
                                track_tie_margin=TK.MARGIN, track_min_identity=TK.MIN_ID, track_min_run_bp=TK.MIN_RAW),
                input=args.get('input'), n_input_records=n_input, n_ends=n_ends)
