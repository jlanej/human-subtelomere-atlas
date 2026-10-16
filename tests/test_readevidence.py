"""Unit tests for subtel_atlas.readevidence on a tiny synthetic fixture (generated deterministically; < 1 MB).

Synthetic CHM13-like ends in tip frame (both p ends):
  chr1p (recipient R = proximal comparator P) = telomere + 10 kb P-private + 30 kb shared block S + 40 kb proximal
  chr2p (donor D)                             = telomere + 10 kb D-private + S with 1.5% substitutions + 20 kb tail
Carrier read = telomere + D-private + D copy of S[:15 kb] + P copy of S[15 kb:] + 30 kb proximal (switch at sub 25,000);
native read  = telomere + P-private + S + 30 kb proximal.
"""
import os, sys, json, shutil
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'atlas'))
from subtel_atlas import readevidence as RE   # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEL = 'CCCTAA' * 500


def _rand(rng, n): return ''.join(rng.choice(list('ACGT'), n))


def _mutate(rng, s, rate):
    s = list(s); idx = rng.choice(len(s), int(rate * len(s)), replace=False)
    for i in idx: s[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)


@pytest.fixture(scope='module')
def fx(tmp_path_factory):
    rng = np.random.default_rng(7)
    S = _rand(rng, 30000); Dblk = _mutate(rng, S, 0.015)
    Ppriv, Dpriv, Pprox, Dtail = _rand(rng, 10000), _rand(rng, 10000), _rand(rng, 40000), _rand(rng, 20000)
    P = TEL + Ppriv + S + Pprox; D = TEL + Dpriv + Dblk + Dtail
    d = tmp_path_factory.mktemp('re')
    fa = d / 'ends.fa'; fa.write_text(f'>chr1p\n{P}\n>chr2p\n{D}\n')
    gen = d / 'genome.fa'; gen.write_text(f'>chr1\n{P + _rand(rng, 20000)}\n>chr2\n{D}\n')
    carrier = TEL[:1500] + Dpriv + Dblk[:15000] + S[15000:] + Pprox[:30000]
    native = TEL[:1500] + Ppriv + S + Pprox[:30000]
    err = np.random.default_rng(11)
    reads = {'carrier': _mutate(err, carrier, 0.005), 'native': _mutate(err, native, 0.005)}
    ends = RE.Chm13Ends(str(fa))
    claim = dict(claim_id='synthetic', end_id='synthetic', R='chr1p', D='chr2p', P='chr1p', recipient_arm='1p', pattern='DP',
                 switch_sub=25000, switch_sub_distal=None, p_switch=3000 + 10000 + 15000, d_lo=3000, d_hi=28000, window_sub_len=80000,
                 core_structure='2p|home', contig='chr1')
    crit = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json'))
    return dict(ends=ends, claim=claim, crit=crit, reads=reads, genome=str(gen), d=d)


def test_parse_cs():
    R, Q, mm, diff = RE.parse_cs(':3*ag:2+tt:1-c:2', 10, 0)
    assert list(R) == [10, 11, 12, 13, 14, 15, 16, 18, 19]
    assert list(Q) == [0, 1, 2, 3, 4, 5, 8, 9, 10]
    assert mm == [(13, 3, 'A', 'G')]


def test_changepoint_dp_and_null():
    cp = RE.changepoint([1] * 30 + [0] * 40, 'DP', 0.05, 1.92)
    assert cp['cuts'] == [30] and cp['ci'][0][0] <= 30 <= cp['ci'][0][1]
    assert (cp['ll'] - cp['ll_allP']) / RE.LN10 > 20 and (cp['ll'] - cp['ll_allD']) / RE.LN10 > 20
    cp0 = RE.changepoint([0] * 70, 'DP', 0.05, 1.92)
    assert (cp0['ll'] - cp0['ll_allP']) <= 0
    cp2 = RE.changepoint([0] * 20 + [1] * 25 + [0] * 30, 'PDP', 0.05, 1.92)
    assert cp2['cuts'] == [20, 45]


def test_binom_tail():
    assert RE._log10_binom_tail(0, 10, 0.05) == 0.0
    assert abs(RE._log10_binom_tail(10, 10, 0.5) - 10 * np.log10(0.5)) < 1e-9


def test_pairmap_sites_and_private(fx):
    ctx = RE.ClaimContext(fx['claim'], fx['ends'], fx['crit'], mask_cache=None)
    s = ctx.pm.summary()
    assert s['n_blocks'] == 1 and 24000 <= s['homologous_bp'] <= 30500   # D window ends at d_hi + 10 kb
    assert 150 <= s['n_snv_sites'] <= 450                # ~450 substitutions before filters
    assert s['d_private_bp'] >= 9000                     # D-private sequence of the donor block is recognised


def test_carrier_and_native_reads(fx):
    ctx = RE.ClaimContext(fx['claim'], fx['ends'], fx['crit'], mask_cache=None)
    r, u, t = RE.evaluate_read(ctx, 'carrier', fx['reads']['carrier'], side='L')
    assert r['terminal'] and r['switch_pass'], r
    assert r['switch_ci_sub_lo'] - 2000 <= 25000 <= r['switch_ci_sub_hi'] + 2000 and r['switch_concordant']
    assert r['seg_D_purity'] >= 0.95 and r['seg_P_purity'] >= 0.95
    n, _, _ = RE.evaluate_read(ctx, 'native', fx['reads']['native'], side='L')
    assert n['terminal'] and not n['switch_pass']
    assert n['donor_interval_units'] >= 10 and n['donor_interval_P_purity'] >= 0.95
    r2, _, _ = RE.evaluate_read(ctx, 'carrier', fx['reads']['carrier'], side='L')
    assert {k: v for k, v in r.items()} == {k: v for k, v in r2.items()}       # deterministic


def test_orientation_from_telomere(fx):
    t, o, tb = RE.orient(RE.rc(fx['reads']['carrier']), None, fx['crit']['thresholds'])
    assert o == '-' and tb['telo_len'] >= 1000 and t.startswith('CCCTAA')


def test_finalize_verdicts(fx):
    ctx = RE.ClaimContext(fx['claim'], fx['ends'], fx['crit'], mask_cache=None)
    r, _, _ = RE.evaluate_read(ctx, 'carrier', fx['reads']['carrier'], side='L')
    n, _, _ = RE.evaluate_read(ctx, 'native', fx['reads']['native'], side='L')
    L = 100000
    anc = [dict(chrom='chr1', chrom_len=L, strand='+', ref_start=30000, ref_end=55000, q_start=0, q_end=25000, mapq=60, supp=False, aln_len=25000)]
    whole = [dict(chrom='chr1', chrom_len=L, strand='+', ref_start=0, ref_end=55000, q_start=0, q_end=r['read_len'], mapq=60, supp=False,
                  aln_len=r['read_len'])]
    assert RE.finalize(ctx, dict(r), whole, anc)['verdict'] == 'conclusive'
    assert RE.finalize(ctx, dict(r), whole, [dict(anc[0], mapq=3)])['verdict'] == 'consistent'
    fold = whole + [dict(whole[0], strand='-', q_start=1000, q_end=5000, aln_len=4000)]
    assert RE.finalize(ctx, dict(r), fold, anc)['verdict'] == 'consistent'
    wn = [dict(whole[0], q_end=n['read_len'], aln_len=n['read_len'])]
    assert RE.finalize(ctx, dict(n), wn, anc)['verdict'] == 'contradicts'
    v = RE.claim_verdict([dict(verdict='conclusive', switch_concordant=True, dataset='a'),
                          dict(verdict='conclusive', switch_concordant=True, dataset='b')], fx['crit'], same_hap_key='contig')
    assert v['verdict'] == 'confirmed' and v['multi_dataset']


@pytest.mark.skipif(not RE.find_minimap2(), reason='minimap2 executable not on PATH')
def test_run_claim_genome_anchor(fx):
    reads = [(k, v, None) for k, v in fx['reads'].items()]
    rows, doc = RE.run_claim(fx['claim'], reads, fx['ends'], fx['crit'], fx['genome'], side='L', dataset='synthetic')
    v = {r['read']: r['verdict'] for r in rows}
    assert v == {'carrier': 'conclusive', 'native': 'contradicts'}, [(r['read'], r['anchor_mapq'], r['anchor_len']) for r in rows]
    assert doc['criteria_version'] == 'v1' and len(doc['sites']) == doc['context']['units_DP']['n_snv_sites']


def test_packaged_criteria_match_repository():
    a = json.load(open(os.path.join(REPO, 'evidence', 'criteria_v1.json')))
    b = RE.load_criteria()
    assert RE.criteria_hash(a) == RE.criteria_hash(b) and a == b
