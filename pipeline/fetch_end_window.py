"""Fetch the tip-frame window of catalogue ends from their public sources and verify each against the released SHA-256.
Usage: python pipeline/fetch_end_window.py END_ID [END_ID ...] > windows.fa
       python pipeline/fetch_end_window.py --ends-file ids.txt | --sample N [--seed S] | --all
                                           [--threads 4] [--out windows.fa] [--report report.tsv]
Reads data/catalogue/ends.tsv (assembly_url, contig, win_start, win_end, strand, window_sha256). Sources:
  * Release 2, HG002 Q100: bgzipped FASTA with remote .fai/.gzi; CHM13 v2.0: plain FASTA with remote .fai (range-readable copy);
  * pedigree (pedigree_aux): plain FASTA on the UCSC server; the contig column holds the chromosome, the sequence name is taken from
    the committed index copies data/inputs/subtel/pedigree/<sample>_<hap>.v1.1.fai (passed to samtools with --fai-idx).
Regions are range-read with `samtools faidx -r` (htslib with libcurl), one call per source file, --threads files in parallel,
3 attempts each. The window is reverse-complemented for R ends ('-' strand). --sample draws a seeded sample stratified by set
(and, for Release 2, by named/unnamed contig and extended window). Exit status 1 if any window is missing or mismatches."""
import os, sys, argparse, subprocess, hashlib, time, tempfile, concurrent.futures as cf, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S3 = os.environ.get('ATLAS_S3_HTTPS', 'https://human-pangenomics.s3.amazonaws.com')
PED_FAI = os.path.join(REPO, 'data', 'inputs', 'subtel', 'pedigree')
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def stratum(r):
    if r['set'] != 'HPRC_R2': return r['set']
    s = 'R2_named' if str(r['contig']).split('#')[-1].startswith('CM') else 'R2_unnamed'
    return s + ('_extended' if bool(r['window_extended']) else '')

def sample(E, n, seed):
    E = E.assign(stratum=[stratum(r) for r in E.to_dict('records')]); g = E.groupby('stratum'); k = max(1, n // g.ngroups)
    S = pd.concat([x.sample(min(k, len(x)), random_state=seed) for _, x in g])
    if len(S) < n: S = pd.concat([S, E.drop(S.index).sample(n - len(S), random_state=seed)])
    return S

def pedigree_names(sample_, hap):
    fai = os.path.join(PED_FAI, f'{sample_}_{hap}.v1.1.fai')
    return fai, {l.split('\t')[0].split('.')[1]: l.split('\t')[0] for l in open(fai)}

def fetch_source(url, rows):
    """rows: list of (end_id, seqname, start, end, strand); returns {end_id: seq}"""
    extra = []; regs = {}
    if rows[0][5] == 'pedigree_aux':
        fai, names = pedigree_names(rows[0][6], rows[0][7]); extra = ['--fai-idx', fai]
    for e, contig, a, b, st, set_, *_ in rows:
        name = names[contig] if set_ == 'pedigree_aux' else contig
        regs[f'{name}:{a}-{b}'] = (e, st)
    with tempfile.NamedTemporaryFile('w', suffix='.regions', delete=False) as f: f.write('\n'.join(regs) + '\n'); rf = f.name
    out = {}
    for t in range(3):
        r = subprocess.run(['samtools', 'faidx', *extra, '-r', rf, url], capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.count('>') == len(regs): break
        time.sleep(5 * (t + 1))
    os.unlink(rf)
    name = None; buf = []
    for line in r.stdout.splitlines() + ['>']:
        if line.startswith('>'):
            if name in regs:
                e, st = regs[name]; s = ''.join(buf).upper(); out[e] = s.translate(comp)[::-1] if st == '-' else s
            name = line[1:].strip(); buf = []
        else: buf.append(line)
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('ends', nargs='*'); ap.add_argument('--ends-file'); ap.add_argument('--sample', type=int)
    ap.add_argument('--seed', type=int, default=20261016); ap.add_argument('--all', action='store_true'); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--out'); ap.add_argument('--report'); a = ap.parse_args()
    E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False,
                    usecols=['end_id', 'set', 'sample', 'haplotype', 'assembly_url', 'contig', 'win_start', 'win_end', 'strand',
                             'window_sha256', 'window_len', 'window_extended'])
    if a.all: S = E
    elif a.sample: S = sample(E, a.sample, a.seed)
    else:
        ids = list(a.ends) + ([l.strip() for l in open(a.ends_file) if l.strip()] if a.ends_file else [])
        S = E.set_index('end_id').loc[ids].reset_index()
    S = S.copy(); S['stratum'] = [stratum(r) for r in S.to_dict('records')]
    S['url'] = S.assembly_url.str.replace('s3://human-pangenomics', S3, regex=False)
    jobs = {u: [(r.end_id, r.contig, int(r.win_start), int(r.win_end), r.strand, r.set, r.sample, r.haplotype) for r in g.itertuples()]
            for u, g in S.groupby('url')}
    got = {}
    with cf.ThreadPoolExecutor(max(1, a.threads)) as ex:
        for res in ex.map(lambda u: fetch_source(u, jobs[u]), sorted(jobs)): got.update(res)
    rep = []; fo = open(a.out, 'w') if a.out else sys.stdout
    for r in S.itertuples():
        s = got.get(r.end_id); ok = s is not None and hashlib.sha256(s.encode()).hexdigest() == r.window_sha256
        rep.append(dict(end_id=r.end_id, set=r.set, stratum=r.stratum, source=r.url, window_len=len(s) if s else 0,
                        status='OK' if ok else ('MISSING' if s is None else 'MISMATCH')))
        if s: fo.write(f'>{r.end_id}\n' + '\n'.join(s[i:i + 100] for i in range(0, len(s), 100)) + '\n')
    if a.out: fo.close()
    R = pd.DataFrame(rep)
    if a.report: R.to_csv(a.report, sep='\t', index=False)
    sys.stderr.write(R.groupby(['stratum', 'status']).size().to_string() + '\n')
    sys.exit(0 if (R.status == 'OK').all() else 1)

if __name__ == '__main__':
    main()
