#!/usr/bin/env python3
"""Build data/metadata/samples.tsv: population, superpopulation and sex for HPRC Release 2 samples.

Sources (public):
  1. HPRC Release 2 sample metadata (human-pangenomics/hprc_intermediate_assembly, data_tables/sample)
  2. 1000 Genomes 3,202-sample pedigree/population table (IGSR FTP, 1000G_2504_high_coverage)
  3. IGSR population -> superpopulation table (phase3 20131219.populations.tsv)
Raw downloads are cached in data/metadata/raw/ with sha256 recorded in data/metadata/raw/SHA256SUMS.
Usage: build_samples_tsv.py [--tips TIPS_META_TSV] [--offline]
"""
import argparse, hashlib, os, sys, urllib.request, datetime
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
RAW = os.path.join(REPO, 'data', 'metadata', 'raw')
OUT = os.path.join(REPO, 'data', 'metadata', 'samples.tsv')
URLS = {
    'hprc_release2_sample_metadata.csv': 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/sample/hprc_release2_sample_metadata.csv',
    '20130606_g1k_3202_samples_ped_population.txt': 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt',
    '20131219.populations.tsv': 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv',
}
# populations absent from the 1kGP phase-3 population table (assigned from the HPRC descriptor)
EXTRA_POP = {
    'MKK': ('AFR', 'Maasai in Kinyawa, Kenya (HapMap3)'),
    'ASL': ('AFR', 'African Americans living in St. Louis, Missouri (WashU pedigree)'),
}
# GIAB HG005 (NA26107) has no population code in the HPRC table; GIAB describes the Chinese trio as Han Chinese.
SPECIAL = {'HG002': ('ASJ_GIAB', 'EUR', 'GIAB Ashkenazi Jewish trio son (no 1kGP population code)'),
           'HG005': ('CHN_GIAB', 'EAS', 'GIAB Han Chinese trio son (no 1kGP population code)')}


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fetch(offline):
    os.makedirs(RAW, exist_ok=True)
    log = []
    for fn, url in URLS.items():
        p = os.path.join(RAW, fn)
        if not (offline and os.path.exists(p)):
            tmp = p + '.tmp'
            urllib.request.urlretrieve(url, tmp)
            os.replace(tmp, p)
            log.append((fn, url, datetime.date.today().isoformat()))
    with open(os.path.join(RAW, 'SHA256SUMS.tmp'), 'w') as f:
        for fn in sorted(URLS):
            f.write(f'{sha256(os.path.join(RAW, fn))}  {fn}\n')
    os.replace(os.path.join(RAW, 'SHA256SUMS.tmp'), os.path.join(RAW, 'SHA256SUMS'))
    if log:
        with open(os.path.join(RAW, 'DOWNLOADS.tsv.tmp'), 'w') as f:
            f.write('file\turl\tdownload_date\n')
            for r in log:
                f.write('\t'.join(r) + '\n')
        os.replace(os.path.join(RAW, 'DOWNLOADS.tsv.tmp'), os.path.join(RAW, 'DOWNLOADS.tsv'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tips', default=os.path.join(os.environ.get('ATLAS_INPUT_DIR', '.'), 'tipdup', 'tips_meta.tsv'))
    ap.add_argument('--offline', action='store_true', help='reuse cached raw files')
    a = ap.parse_args()
    fetch(a.offline)
    h = pd.read_csv(os.path.join(RAW, 'hprc_release2_sample_metadata.csv'))
    ped = pd.read_csv(os.path.join(RAW, '20130606_g1k_3202_samples_ped_population.txt'), sep=' ')
    pops = pd.read_csv(os.path.join(RAW, '20131219.populations.tsv'), sep='\t')
    pops = pops.dropna(subset=['Population Code'])
    p2s = dict(zip(pops['Population Code'], pops['Super Population']))
    for k, (sp, _) in EXTRA_POP.items():
        p2s.setdefault(k, sp)
    pedd = ped.set_index('SampleID')
    tips_samples = set(pd.read_csv(a.tips, sep='\t', usecols=['sample'])['sample']) if os.path.exists(a.tips) else set()
    rows = []
    for _, r in h.iterrows():
        s = r['sample_id']
        src = []
        pop = r['population_abbreviation'] if isinstance(r['population_abbreviation'], str) else None
        if pop:
            src.append('HPRC_R2_metadata')
        sp = p2s.get(pop) if pop else None
        if s in pedd.index:
            pp, ps = pedd.loc[s, 'Population'], pedd.loc[s, 'Superpopulation']
            if pop and pp != pop:
                print(f'WARN population mismatch {s}: HPRC {pop} vs 1kGP {pp}', file=sys.stderr)
            pop = pop or pp
            sp = ps if pp == pop else sp
            src.append('1kGP_3202_ped')
        elif sp:
            src.append('IGSR_population_table' if pop not in EXTRA_POP else 'manual_superpop_from_descriptor')
        if s in SPECIAL:
            pop, sp, note = SPECIAL[s]
            src.append('manual:' + note)
        sex = r['sex'] if isinstance(r['sex'], str) else None
        if s in pedd.index:
            ps = {1: 'male', 2: 'female'}.get(int(pedd.loc[s, 'Sex']))
            if sex and ps and ps != sex:
                print(f'WARN sex mismatch {s}: {sex} vs {ps}', file=sys.stderr)
            sex = sex or ps
        rows.append(dict(sample=s, population=pop or 'NA', superpopulation=sp or 'NA', sex=sex or 'NA',
                         source=';'.join(src) or 'NA',
                         population_descriptor=r['population_descriptor'] if isinstance(r['population_descriptor'], str) else 'NA',
                         biosample_id=r['biosample_id'], alternative_id=r['alternative_id'] if isinstance(r['alternative_id'], str) else 'NA',
                         paternal_id=r['paternal_id'] if isinstance(r['paternal_id'], str) else 'NA',
                         maternal_id=r['maternal_id'] if isinstance(r['maternal_id'], str) else 'NA',
                         in_1kgp_3202=s in pedd.index, has_telomeric_tips=s in tips_samples))
    df = pd.DataFrame(rows).sort_values('sample')
    miss = tips_samples - set(df['sample'])
    if miss:
        print('WARN tip samples absent from HPRC table:', sorted(miss), file=sys.stderr)
    # relatedness within the set: parent/child pairs both present
    ss = set(df['sample'])
    rel = df[df.paternal_id.isin(ss) | df.maternal_id.isin(ss)]
    df['related_in_set'] = df['sample'].isin(set(rel['sample']) | set(rel.paternal_id) | set(rel.maternal_id))
    tmp = OUT + '.tmp'
    df.to_csv(tmp, sep='\t', index=False)
    os.replace(tmp, OUT)
    print(f'wrote {OUT}: {len(df)} samples; with tips {df.has_telomeric_tips.sum()}; '
          f'superpop counts (tip samples) {df[df.has_telomeric_tips].superpopulation.value_counts().to_dict()}; '
          f'related_in_set {df.related_in_set.sum()}')


if __name__ == '__main__':
    main()
