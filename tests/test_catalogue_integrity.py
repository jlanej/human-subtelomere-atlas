"""Integrity tests for the released allele catalogue (data/catalogue). Run: python -m pytest tests/test_catalogue_integrity.py
(or python tests/test_catalogue_integrity.py). Needs only pandas and the released files."""
import os, gzip, json, hashlib
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAT = os.path.join(REPO, 'data', 'catalogue')
rd = lambda f: pd.read_csv(os.path.join(CAT, f), sep='\t', low_memory=False)

def test_sha256sums():
    for line in open(os.path.join(CAT, 'SHA256SUMS')):
        h, p = line.rstrip('\n').split('  ', 1)
        assert hashlib.sha256(open(os.path.join(REPO, p), 'rb').read()).hexdigest() == h, p

def test_ids_and_keys():
    E, SA, HG = rd('ends.tsv'), rd('alleles_structural.tsv'), rd('haplogroups.tsv')
    assert E.end_id.is_unique and SA.allele_id.is_unique and HG.haplogroup_id.is_unique
    assert set(E.structural_allele_id.dropna()) <= set(SA.allele_id)
    assert set(E.haplogroup_id.dropna()) <= set(HG.haplogroup_id)
    assert set(SA.representative_end_id) <= set(E.end_id) and set(HG.medoid_end_id) <= set(E.end_id)
    core = E[E.set != 'pedigree_aux']
    n = core.groupby('structural_allele_id').size()
    assert (SA.set_index('allele_id').n_ends == n.reindex(SA.allele_id).fillna(0).values).all()
    # every structural allele is used by its representative with the same core structure
    rep = E.set_index('end_id').loc[SA.representative_end_id]
    assert (rep.core_structure.values == SA.structure.values).all()

def test_frequencies():
    F = rd('allele_frequencies_by_superpop.tsv')
    g = F.groupby(['allele_type', 'arm', 'superpopulation']).agg(n=('n_haplotype_ends', 'sum'), d=('n_arm_ends_assessed', 'first'),
                                                                 dn=('n_arm_ends_assessed', 'nunique'), f=('frequency', 'sum'))
    assert (g.n == g.d).all() and (g.dn == 1).all() and ((g.f - 1).abs() < 0.01).all()   # frequencies are rounded to 4 decimals

def test_representatives_and_json():
    SA, HG = rd('alleles_structural.tsv'), rd('haplogroups.tsv')
    heads = [l[1:].split()[0] for l in gzip.open(os.path.join(CAT, 'representatives.fa.gz'), 'rt') if l.startswith('>')]
    ids = [i for l in gzip.open(os.path.join(CAT, 'representatives.fa.gz'), 'rt') if l.startswith('>')
           for t in l.split() if t.startswith('ids=') for i in t[4:].split(',')]
    want = list(SA[SA.representative_in_fasta].allele_id) + list(HG[HG.representative_in_fasta].haplogroup_id)
    assert sorted(ids) == sorted(want) and len(heads) == len(set(heads))
    J = json.load(open(os.path.join(CAT, 'catalogue.json')))
    ids = {a['id'] for arm in J['arms'].values() for a in arm['structural_alleles']}
    assert ids == set(SA.allele_id)

if __name__ == '__main__':
    for f in [test_sha256sums, test_ids_and_keys, test_frequencies, test_representatives_and_json]:
        f(); print('ok', f.__name__)
