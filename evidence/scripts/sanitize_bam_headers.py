#!/usr/bin/env python3
"""Write copies of BAM files whose @PG/@CO header lines carry no local paths (for publishing evidence BAMs).

Every absolute path in a @PG CL or @CO value is replaced by its file name; the input files are not modified.
Usage: sanitize_bam_headers.py --out <dir> in1.bam [in2.bam ...]   (writes <dir>/<name>.bam + .bai)
"""
import argparse, os, re, sys
import pysam

PATH = re.compile(r'(?:/[^\s/\'"]+)+/([^\s/\'"]+)')


def clean(v):
    return PATH.sub(lambda m: m.group(1), v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True)
    ap.add_argument('bams', nargs='+')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for fn in a.bams:
        out = os.path.join(a.out, os.path.basename(fn))
        if os.path.abspath(out) == os.path.abspath(fn): sys.exit('output would overwrite the input')
        with pysam.AlignmentFile(fn) as f:
            h = f.header.to_dict()
            for pg in h.get('PG', []):
                if 'CL' in pg: pg['CL'] = clean(pg['CL'])
            h['CO'] = [clean(c) for c in h.get('CO', [])]
            with pysam.AlignmentFile(out, 'wb', header=h) as o:
                for r in f.fetch(until_eof=True): o.write(r)
        pysam.index(out)
        print(out, file=sys.stderr)


if __name__ == '__main__':
    main()
