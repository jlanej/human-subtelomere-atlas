#!/usr/bin/env python3
"""Fetch the outer 500 kb of every chromosome end of the T2T great-ape primary assemblies (Yoo et al. 2025,
doi:10.1038/s41586-025-08816-3) from NCBI (Datasets v2 sequence reports + E-utilities efetch with seq_start/seq_stop).
Y chromosomes are skipped. Outputs: work/allele_ages/apes/<species>.ends500k.fa (sequence names <species>|<chr>|<L/R>,
L = first 500 kb, R = last 500 kb, both in assembly orientation) and analyses/allele_ages/results/ape_ends_provenance.tsv.
Idempotent (skips species whose fasta exists). Requests are rate-limited to < 3 per second (NCBI policy without API key).
"""
import os, sys, json, time, datetime, urllib.request, hashlib
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
OUT = os.path.join(REPO, 'work', 'allele_ages', 'apes'); os.makedirs(OUT, exist_ok=True)
PROV = os.path.join(REPO, 'analyses', 'allele_ages', 'results', 'ape_ends_provenance.tsv')
W = 500_000
ACC = {'chimpanzee': 'GCA_028858775.2', 'bonobo': 'GCA_029289425.2', 'gorilla': 'GCA_029281585.2',
       'sumatran_orangutan': 'GCA_028885655.2', 'bornean_orangutan': 'GCA_028885625.2'}
DS = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/sequence_reports?page_size=200'
EF = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text&seq_start={}&seq_stop={}'


def get(url, tries=4):
    for t in range(tries):
        try:
            time.sleep(0.4)
            return urllib.request.urlopen(url, timeout=120).read().decode()
        except Exception as e:
            if t == tries - 1:
                raise
            time.sleep(5 * (t + 1))


rows = []
for sp, acc in ACC.items():
    rep = json.loads(get(DS.format(acc)))['reports']
    chs = [r for r in rep if r.get('role') == 'assembled-molecule' and not r['chr_name'].startswith('Y')]
    fa = os.path.join(OUT, f'{sp}.ends500k.fa')
    have = os.path.exists(fa)
    out = None if have else open(fa + '.tmp', 'w')
    for r in chs:
        L = int(r['length']); g = r['genbank_accession']
        for end, (s, e) in (('L', (1, min(W, L))), ('R', (max(1, L - W + 1), L))):
            if not have:
                txt = get(EF.format(g, s, e))
                seq = ''.join(txt.split('\n')[1:]).strip()
                assert len(seq) == e - s + 1, (sp, g, end, len(seq))
                out.write(f'>{sp}|{r["sequence_name"]}|{end}\n')
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + '\n')
            rows.append(dict(species=sp, assembly=acc, sequence_name=r['sequence_name'], chr_name=r['chr_name'],
                             genbank_accession=g, chr_length=L, end=end, start=s, stop=e))
    if not have:
        out.close(); os.replace(fa + '.tmp', fa)
    print(sp, len(chs), 'chromosomes', file=sys.stderr)
import pandas as pd
d = pd.DataFrame(rows); d['fetch_date'] = datetime.date.today().isoformat()
if os.path.exists(PROV):
    old = pd.read_csv(PROV, sep='\t'); d['fetch_date'] = d.species.map(old.drop_duplicates('species').set_index('species').fetch_date).fillna(d.fetch_date)
d['sha256_species_fasta'] = d.species.map({sp: hashlib.sha256(open(os.path.join(OUT, f'{sp}.ends500k.fa'), 'rb').read()).hexdigest() for sp in ACC})
d.to_csv(PROV, sep='\t', index=False)
print('ends', len(d))
