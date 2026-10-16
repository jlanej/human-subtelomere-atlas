#!/usr/bin/env python3
"""Fetch telomere-reaching reads by region from public, indexed long-read alignments (HTTP range reads; samtools >= 1.16).

  hg002  : HG002 Q100 v1.1 read alignments to the HG002 Q100 v1.1 diploid assembly (human-pangenomics S3), one BAM slice per
           (dataset, contig end): the outer --flank bp of the contig end, primary alignments.
  ont1kg : 1000 Genomes ONT Vienna CRAMs aligned to CHM13 v2.0 (Schloissnig, Pani et al. 2025, Nature, doi:10.1038/s41586-025-09290-7), one
           slice per (sample, CHM13 chromosome end); decoding needs the CHM13 v2.0 FASTA (--chm13).

Output: work/readevidence/<source>/<dataset or sample>/<contig>_<side>.bam (+ .bai) and a manifest TSV with the URL, region,
read count and SHA-256 of the slice. Inputs are given on the command line or come from the catalogue; no personal paths.
"""
import argparse, os, sys, subprocess, hashlib, csv, shutil, urllib.request
from concurrent.futures import ThreadPoolExecutor

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
Q100 = 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/'
Q100_SETS = {'ont_ul_dorado': 'ont_r10_ul_dorado/hg002v1.1_ont_r10_ul_dorado.pri.bam',
             'hifi_revio_3cells': 'hifi_revio_3cell/hg002v1.1_hifi_revio_3cells.bam',
             'ont_duplex': 'ont_r10_duplex/hg002v1.1_ont_r10_duplex.pri.bam'}
ONT1KG = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/t2t/{s}.t2t.cram'


def sha(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def index_for(url, cache, local_dirs=()):
    """local copy of the .bai/.crai (downloaded once into cache)"""
    base = os.path.basename(url) + ('.crai' if url.endswith('.cram') else '.bai')
    for d in local_dirs:
        for c in (os.path.join(d, base),):
            if os.path.exists(c): return c
    os.makedirs(cache, exist_ok=True); fn = os.path.join(cache, base)
    if not os.path.exists(fn):
        tmp = fn + '.part'; urllib.request.urlretrieve(url + ('.crai' if url.endswith('.cram') else '.bai'), tmp); os.replace(tmp, fn)
    return fn


def fetch(url, idx, region, out, ref=None, tries=3):
    if os.path.exists(out) and os.path.exists(out + '.bai'): return out
    os.makedirs(os.path.dirname(out), exist_ok=True); tmp = out + '.part.bam'
    cmd = ['samtools', 'view', '-b', '-F', '0x900', '-o', tmp]
    if ref: cmd += ['-T', ref]
    cmd += ['-X', url, idx, region]
    for i in range(tries):
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(out))
        if r.returncode == 0: break
        if i == tries - 1: raise RuntimeError(f'{region}: {r.stderr[-300:]}')
    subprocess.run(['samtools', 'index', tmp], check=True)
    os.replace(tmp, out); os.replace(tmp + '.bai', out + '.bai')
    return out


def count(fn):
    return int(subprocess.run(['samtools', 'view', '-c', fn], capture_output=True, text=True, check=True).stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('source', choices=['hg002', 'ont1kg'])
    ap.add_argument('--ends', required=True, help='TSV with columns key, contig, side, contig_len (hg002) or sample, chm13_end (ont1kg)')
    ap.add_argument('--datasets', default='ont_ul_dorado,hifi_revio_3cells', help='hg002 datasets')
    ap.add_argument('--flank', type=int, default=5000, help='bases from the chromosome end to query (reads reaching the end overlap it)')
    ap.add_argument('--chm13', help='CHM13 v2.0 FASTA (bgzip + .fai/.gzi) for CRAM decoding')
    ap.add_argument('--chm13-fai', help='.fai of CHM13 v2.0 (chromosome lengths); default <chm13>.fai')
    ap.add_argument('--index-dir', action='append', default=[], help='directory with local copies of the public .bai/.crai files')
    ap.add_argument('--out', default=os.path.join(REPO, 'work', 'readevidence'))
    ap.add_argument('-j', '--jobs', type=int, default=6)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.ends), delimiter='\t'))
    jobs = []
    if a.source == 'hg002':
        for ds in a.datasets.split(','):
            url = Q100 + Q100_SETS[ds]; idx = index_for(url, os.path.join(a.out, 'index'), a.index_dir)
            for r in rows:
                L = int(r['contig_len'])
                reg = f"{r['contig']}:1-{a.flank}" if r['side'] == 'L' else f"{r['contig']}:{L - a.flank + 1}-{L}"
                jobs.append((ds, url, idx, reg, os.path.join(a.out, 'hg002', ds, f"{r['contig']}_{r['side']}.bam"), None))
    else:
        fai = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(a.chm13_fai or a.chm13 + '.fai')}
        for r in rows:
            url = ONT1KG.format(s=r['sample']); idx = index_for(url, os.path.join(a.out, 'index'), a.index_dir)
            c = r['chm13_end'][:-1]; L = fai[c]
            reg = f'{c}:1-{a.flank}' if r['chm13_end'].endswith('p') else f'{c}:{L - a.flank + 1}-{L}'
            jobs.append((r['sample'], url, idx, reg, os.path.join(a.out, 'ont1kg', r['sample'], f"{r['chm13_end']}.bam"), a.chm13))

    def go(j):
        ds, url, idx, reg, out, ref = j
        fetch(url, idx, reg, out, ref)
        return dict(set=ds, url=url, region=reg, file=os.path.relpath(out, REPO), n_reads=count(out), sha256=sha(out))
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(go, jobs))
    man = os.path.join(a.out, f'{a.source}_fetch_manifest.tsv')
    with open(man, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]), delimiter='\t'); w.writeheader(); w.writerows(res)
    print(f'{len(res)} slices, {sum(r["n_reads"] for r in res)} records -> {man}', file=sys.stderr)


if __name__ == '__main__':
    main()
