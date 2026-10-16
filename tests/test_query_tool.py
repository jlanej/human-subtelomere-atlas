"""Offline tests of the subtel-atlas query tool (atlas/subtel_atlas). Run: python -m pytest tests/test_query_tool.py

Fixtures (tests/fixtures, < 0.4 MB): four catalogue windows (SA-20q-01 representative, a 6p end, CHM13 1p, a partial 20q window) and
six CHM13 v2.0 ends reduced to their telomeric 150 kb (the remaining 350 kb masked with N, so coordinates match the full ends).
Golden outputs (expected_windows.<backend>.tsv) were produced by this version; a change in them means a change in results.
The catalogue tables and representative windows are read from data/catalogue. With a verified full reference directory
($SUBTEL_ATLAS_REFS, from `subtel-atlas fetch-refs`) test_full_reference_matches_catalogue also compares with data/catalogue/ends.tsv.
"""
import os, io, sys, gzip, json, random, shutil, subprocess
import pytest
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(REPO, 'tests', 'fixtures')
CATDIR = os.path.join(REPO, 'data', 'catalogue')


def family_of(arm):
    """Published paralog family token of a CHM13 arm (data/catalogue/paralog_families.tsv); the arm itself if it is in none."""
    F = pd.read_csv(os.path.join(CATDIR, 'paralog_families.tsv'), sep='\t')
    return next((f for f, m in zip(F.family, F.members) if arm in str(m).split(',')), arm)
sys.path.insert(0, os.path.join(REPO, 'atlas'))
from subtel_atlas import query as Q, bins as B, alleles as AL, arms as AR, telomere as T, sketch as S, refs as R, mosaic as M
from subtel_atlas.common import rc, readfa, group_of, group_ends, CHM13_ENDS

HAVE_MM2 = B.find_minimap2() is not None
BACKENDS = ['mappy'] + (['minimap2'] if HAVE_MM2 else [])


@pytest.fixture(scope='session')
def refdir(tmp_path_factory):
    d = tmp_path_factory.mktemp('refs'); shutil.copy(os.path.join(FIX, 'chm13_ends_subset.fa.gz'), d / 'chm13_ends500k.fa.gz')
    for f in R.CAT_FILES: shutil.copy(os.path.join(REPO, 'data', 'catalogue', f), d / f)
    return str(d)


@pytest.fixture(scope='session')
def ref(refdir):
    return Q.Reference(refdir, verify=False)


@pytest.fixture(scope='session')
def windows():
    return [(n, c, s) for n, c, s in Q.iter_records(os.path.join(FIX, 'windows.fa.gz'))]


def cli(*args):
    env = dict(os.environ, PYTHONPATH=os.path.join(REPO, 'atlas') + os.pathsep + os.environ.get('PYTHONPATH', ''))
    return subprocess.run([sys.executable, '-m', 'subtel_atlas'] + list(args), capture_output=True, text=True, env=env)


# ---------------------------------------------------------------- primitives
def test_telomere_caller_synthetic():
    rnd = random.Random(1); sub = ''.join(rnd.choice('ACGT') for _ in range(20000))
    tip = 'GATTACA' * 10 + 'CCCTAA' * 500 + sub
    b = T.find_block(tip)
    # the block start is extended outward while preceding bases continue the repeat phase (here the final A of GATTACA)
    assert (b['junk_len'], b['telo_len']) == (69, 3002)
    assert Q.has_seed(tip) and not Q.has_seed(sub) and T.find_block(sub) is None
    assert Q.call_telomere(rc(tip)) is None                     # G-strand orientation is not a tip-frame telomere


def test_tiling():
    t = B.tile(1000, 100_000)
    assert t[0] == ('0', 1000, 6000) and t[-4:] == [('S5', 1000, 6000), ('S10', 1000, 11000), ('S20', 1000, 21000), ('t20', 1000, 20000)]
    assert [x for x in t if x[0] == '19'][0] == ('19', 96000, 100000)      # last bin kept (>= 2 kb)
    assert not [x for x in B.tile(1000, 22500) if x[0] == '4']           # 1.5 kb remainder dropped


def test_score_cs():
    assert B.score_cs(':100') == (100, 0.0)
    AS, de = B.score_cs(':50*ac:49+ggg')                          # 99 matches, 1 mismatch (-4), 3 bp insertion (-min(6+6, 26+3))
    assert AS == 99 - 4 - 12 and de == float(f'{2 / 101:.4g}')


def test_sketch_deterministic():
    rnd = random.Random(7); s = ''.join(rnd.choice('ACGT') for _ in range(50000))
    a, b = S.sketch(s), S.sketch(rc(s))
    assert (a == b).all() and 3500 < len(a) < 6500                   # canonical k-mers; ~1/10 of 50 kb
    assert abs(float(S.ani_from_jaccard(1.0)) - 1.0) < 1e-12


def test_allele_matching_rules():
    al = [dict(allele_id='SA-x-01', structure='3q|home', sw=[30000], n=10), dict(allele_id='SA-x-02', structure='3q|home', sw=[45000], n=5),
          dict(allele_id='SA-x-03', structure='home', sw=[], n=50)]
    assert AL.match_allele('3q|home', [40000], 92000, True, al)[0]['allele_id'] == 'SA-x-02'     # both within 10 kb: nearest median
    assert AL.match_allele('3q|home', [37500], 92000, True, al)[0]['allele_id'] == 'SA-x-01'     # equidistant: more defining ends
    assert AL.match_allele('3q|home', [40000], 92000, True, al, legacy=True)[0]['allele_id'] == 'SA-x-01'   # previous rule
    assert AL.match_allele('3q|home', [50001], 92000, True, al)[0]['allele_id'] == 'SA-x-02'
    assert AL.match_allele('3q|home', [55001], 92000, True, al) == (None, [])                    # > 10 kb from every median
    assert AL.match_allele('3q|home', [20000], 92000, True, al)[0]['allele_id'] == 'SA-x-01'     # exactly 10 kb: within tolerance
    best, cand = AL.match_allele('3q', [], 25000, False, al)       # partial (25 kb observed): both switches lie beyond 25 kb - 10 kb
    assert best['allele_id'] == 'SA-x-01' and [c['allele_id'] for c in cand] == ['SA-x-01', 'SA-x-02']
    best, cand = AL.match_allele('3q', [], 45000, False, al)            # at 45 kb the 30 kb switch would have been observed
    assert [c['allele_id'] for c in cand] == ['SA-x-02']
    best, cand = AL.match_allele('home', [], 35000, False, al)
    assert best['allele_id'] == 'SA-x-03' and len(cand) == 1


def test_arm_rules():
    assert AR.arm_from_name('chr20', 'R') == '20q' and AR.arm_from_name('chr1_MATERNAL', 'L') == '1p'
    assert AR.arm_from_name('HG002#1#chrX_PATERNAL', 'L') == 'Xp' and AR.arm_from_name('PAN027_pat|chr13|L|1000', 'L') == '13p'
    assert AR.arm_from_name('CM094060.1', 'L') is None and AR.arm_from_name('chr1_chr2_fusion', 'L') is None and AR.arm_from_name('chrUn', 'L') is None
    assert group_of('21p') == 'acroP' and group_of('Yq') == 'PAR2' and group_ends('PAR1') == ['chrXp', 'chrYp']
    assert AR.decide('', '', '', '15p') == ('15p', 'mosaic', 'mosaic_group') and AR.decide('', '', '', '') == ('unassigned', 'none', 'none')


def test_pipeline_imports_shared_code():
    """The catalogue build uses the package implementation (no second copy of the methods)."""
    want = {'catalogue_04_align_bins.py': 'subtel_atlas.bins', 'catalogue_05_arms.py': 'subtel_atlas import arms',
            'catalogue_06_mosaic.py': 'subtel_atlas.mosaic', 'catalogue_07_haplogroups.py': 'subtel_atlas import haplo',
            'catalogue_08_alleles.py': 'subtel_atlas import alleles', 'catalogue_telomere.py': 'subtel_atlas.telomere',
            'catalogue_sketch.py': 'subtel_atlas.sketch', 'catalogue_common.py': 'subtel_atlas.common',
            'catalogue_04t_align_tracks.py': 'subtel_atlas.tracks', 'catalogue_06t_mosaic_tracks.py': 'subtel_atlas import tracks',
            'catalogue_08e_boundary_stable.py': 'subtel_atlas import tracks'}
    for f, s in want.items():
        src = open(os.path.join(REPO, 'pipeline', f)).read()
        assert s in src, f
        assert 'def label_bins' not in src and 'def find_block' not in src and 'def summarise' not in src


def test_window_rule():
    """catalogue window rule (pipeline/catalogue_03): 100 kb, extended to ceil((s0 + 80 kb) / 5 kb) x 5 kb when s0 > 20 kb"""
    assert Q.window_length(19000, 2_000_000) == 100_000 and Q.window_length(20000, 2_000_000) == 100_000
    assert Q.window_length(25000, 2_000_000) == 105_000 and Q.window_length(25001, 2_000_000) == 110_000
    assert Q.window_length(25000, 101_000) == 101_000 and Q.window_length(25000, 2_000_000, tel_to_end=True) == 100_000
    assert Q.window_length(5000, 60_000) == 60_000


def test_checksum_table_matches_release():
    X = R.expected(); sums = {l.split()[1]: l.split()[0] for l in open(os.path.join(REPO, 'data', 'catalogue', 'SHA256SUMS'))}
    for f, h in X['catalogue_files'].items(): assert sums[f'data/catalogue/{f}'] == h, f
    assert list(X['chm13_ends']) == CHM13_ENDS and X['catalogue_version'] == Q.CATALOGUE_VERSION


def test_unverified_reference_is_refused(refdir):
    with pytest.raises(SystemExit): Q.Reference(refdir)                  # fixture CHM13 subset fails the checksums


# ---------------------------------------------------------------- end-to-end on fixtures
@pytest.mark.parametrize('backend', BACKENDS)
def test_windows_golden(ref, windows, backend):
    ends = [Q.window_end(n, c, s) for n, c, s in windows]
    recs = Q.analyse(ends, ref, backend, threads=2)
    buf = io.StringIO(); Q.write_tsv(recs, buf, Q.meta(backend, {}, len(windows), len(recs)))
    assert buf.getvalue() == open(os.path.join(FIX, f'expected_windows.{backend}.tsv')).read()


def test_known_biology_on_fixtures(ref, windows):
    """SA-20q-01 representative: 20q end whose outer ~31 kb best match CHM13 3q, labelled by its paralog family (3q with 19p);
    switch within 1 kb of the catalogue median."""
    be = BACKENDS[-1]; f3q = family_of('3q')
    A = pd.read_csv(os.path.join(CATDIR, 'alleles_structural.tsv'), sep='\t').set_index('allele_id')
    r = Q.analyse([Q.window_end(*windows[0])], ref, be)[0]
    assert r['catalogue_arm'] == '20q' and r['S20_state'] == 'foreign' and r['form'] == f'20q>{f3q}' and r['S20_family'] == f3q
    assert r['S20_best_end'] == 'chr3q'
    assert r['core_structure'] == f'{f3q}|home' and abs(int(r['core_switches_sub']) - int(A.loc['SA-20q-01', 'switches_median'])) <= 1000
    assert r['structural_allele_id'] == 'SA-20q-01' and r['sa_tier'] == 'A' and r['s0'] == 7177
    c = Q.analyse([Q.window_end(*windows[2])], ref, be)[0]                # CHM13 1p: native, own
    assert c['core_structure'] == 'home' and c['S20_state'] == 'own' and c['haplogroup_id'] == 'HG-1p-02'
    p = Q.analyse([Q.window_end(*windows[3])], ref, be)[0]                # short window -> prefix matching
    assert p['full_length'] is False and p['sa_assignment'].startswith('partial') and 'partial_window' in p['qc_flags']


def test_boundary_stability_and_nativity(ref, windows):
    """Per-end -500/+500 bp boundary-shift calls (same computation as pipeline/catalogue_08e) and the allele's nativity class."""
    r = Q.analyse([Q.window_end(*windows[0])], ref, BACKENDS[-1])[0]
    assert isinstance(r['boundary_stable_allele'], bool) and r['allele_minus500'].startswith('SA-20q-')
    assert r['boundary_stable'] == (r['boundary_stable_allele'] and r['boundary_stable_S20'])
    assert ('boundary_unstable' in r['qc_flags']) == (not r['boundary_stable'])
    A = pd.read_csv(os.path.join(CATDIR, 'alleles_structural.tsv'), sep='\t').set_index('allele_id')
    assert r['sa_is_native'] is False
    if 'nativity_class' in A: assert r['sa_nativity_class'] == A.loc['SA-20q-01', 'nativity_class'] != 'native'


def test_contig_mode_finds_both_ends(ref, windows):
    """A synthetic contig = window A (tip frame) + spacer + reverse complement of window B: both telomeric ends are found,
    oriented, named and analysed exactly like the windows themselves."""
    (na, _, a), (nb, _, b) = windows[0], windows[2]
    rnd = random.Random(3); spacer = ''.join(rnd.choice('ACGT') for _ in range(5000))
    contig = a + spacer + rc(b)
    ends = Q.contig_ends('chr20', contig)
    assert [e['end_id'] for e in ends] == ['chr20#L', 'chr20#R'] and ends[0]['seq'] == a[:100000] and ends[1]['seq'] == b[:100000]
    assert ends[0]['arm_hint'] == '20p' and ends[1]['win_end'] == len(contig)
    assert Q.contig_ends('ctg', spacer * 3) == []                               # no telomere -> no end
    ends = Q.contig_ends('ctg', contig, name_arms=False, arm_map={'ctg#L': '20q'})
    assert ends[0]['arm_hint'] == '20q' and ends[1]['arm_hint'] is None


def test_read_mode_simulated(ref, windows):
    """A simulated HiFi-like read (0.5% substitutions) from the 20q window, in G-strand orientation, is oriented and assigned."""
    rnd = random.Random(11); s = list(windows[0][2][:60000])
    for i in range(len(s)):
        if rnd.random() < 0.005: s[i] = rnd.choice('ACGT'.replace(s[i], ''))
    read = rc(''.join(s))
    e = Q.read_end('read1', read)
    assert e is not None and e['strand'] == '-' and abs(e['telo_len'] - 7177) < 200
    e['arm_hint'] = '20q'
    r = Q.analyse([e], ref, BACKENDS[-1])[0]
    assert r['core_structure'] == f"{family_of('3q')}|home" and 'experimental_read_mode' in r['qc_flags']
    assert Q.read_end('noise', ''.join(rnd.choice('ACGT') for _ in range(20000))) is None


def test_json_output_deterministic(ref, windows):
    ends = [Q.window_end(n, c, s) for n, c, s in windows[:2]]
    outs = []
    for _ in range(2):
        buf = io.StringIO(); Q.write_json(Q.analyse(ends, ref, 'mappy'), buf, Q.meta('mappy', {}, 2, 2)); outs.append(buf.getvalue())
    assert outs[0] == outs[1]
    d = json.loads(outs[0]); assert d['catalogue_version'] == 'v0.1' and len(d['ends'][0]['blocks']) >= 2


# ---------------------------------------------------------------- CLI
def test_cli_help_and_version():
    for sub in ([], ['fetch-refs'], ['query'], ['reads'], ['verify-refs']):
        p = cli(*sub, '--help'); assert p.returncode == 0 and 'usage' in p.stdout
    p = cli('--version'); assert 'catalogue v0.1' in p.stdout


def test_cli_query(refdir, tmp_path):
    out = tmp_path / 'o.tsv'
    p = cli('query', '--windows', os.path.join(FIX, 'windows.fa.gz'), '--refs', refdir, '--no-verify-refs', '--aligner', 'mappy', '-o', str(out))
    assert p.returncode == 0, p.stderr
    assert out.read_text() == open(os.path.join(FIX, 'expected_windows.mappy.tsv')).read()
    p = cli('query', os.path.join(FIX, 'windows.fa.gz'), '--refs', str(tmp_path / 'missing'))
    assert p.returncode != 0 and 'fetch-refs' in p.stderr


# ---------------------------------------------------------------- optional: full reference
@pytest.mark.skipif(not os.environ.get('SUBTEL_ATLAS_REFS'), reason='set SUBTEL_ATLAS_REFS to a fetch-refs directory')
def test_full_reference_matches_catalogue(windows):
    ref = Q.Reference(os.environ['SUBTEL_ATLAS_REFS'])
    E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False).set_index('end_id')
    recs = Q.analyse([Q.window_end(n, c, s) for n, c, s in windows], ref, BACKENDS[-1])
    for r in recs:
        e = E.loc[r['query_id']]
        assert r['s0'] == e.junk_len + e.telo_len and r['S20_state'] == e.S20_state and r['core_structure'] == e.core_structure
        assert r['structural_allele_id'] == e.structural_allele_id
