"""Static check of the guided tour (docs/evidence/tour/tour.json): every igv.js stop must point at existing files, every
BAM/BED contig must exist in the stop's reference, every view locus must lie inside its contig, BAMs must pass samtools
quickcheck, and no data file may exceed GitHub's 100 MB limit (warning above 50 MB).

Usage (from the repository root): python website/tour_check.py   -> exit status 1 on any error
"""
import os, sys, json, subprocess, gzip
import pysam

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EV = os.path.join(REPO, 'docs', 'evidence')


def stops(T):
    for ch in T['chapters']:
        if 'gallery' in ch:
            for g in ch['gallery']: yield f"gallery:{g['ref']['id']}", g
        else: yield ch['id'], ch


def main():
    T = json.load(open(os.path.join(EV, 'tour', 'tour.json'))); err, warn, n = [], [], 0
    for name, s in stops(T):
        n += 1; ref = s['ref']
        for k in ('fastaURL', 'indexURL', 'compressedIndexURL'):
            if not os.path.exists(os.path.join(EV, ref[k])): err.append(f'{name}: missing {ref[k]}')
        fai = os.path.join(EV, ref['indexURL'])
        L = dict((l.split('\t')[0], int(l.split('\t')[1])) for l in open(fai)) if os.path.exists(fai) else {}
        for v in s['views']:
            for loc in v['locus']:
                c, rg = loc.rsplit(':', 1); a, b = (int(x.replace(',', '')) for x in rg.split('-'))
                if c not in L or not (1 <= a < b <= L[c]): err.append(f'{name}: locus {loc} outside reference')
        for t in [u for x in s['tracks'] for u in (x['tracks'] if x.get('type') == 'merged' else [x])]:
            f = os.path.join(EV, t['url'])
            if not os.path.exists(f): err.append(f'{name}: missing {t["url"]}'); continue
            if t.get('format') == 'bam':
                if not os.path.exists(os.path.join(EV, t.get('indexURL', t['url'] + '.bai'))): err.append(f'{name}: no index {t["url"]}')
                if subprocess.run(['samtools', 'quickcheck', f]).returncode: err.append(f'{name}: quickcheck {t["url"]}')
                with pysam.AlignmentFile(f) as b:
                    bad = [c for c in b.references if c not in L]
                if bad: err.append(f'{name}: {t["url"]} contigs not in reference: {bad[:3]}')
            elif t.get('format') in ('bed', 'bedgraph'):
                op = gzip.open if f.endswith('.gz') else open
                bad = {l.split('\t')[0] for l in op(f, 'rt') if l.strip() and not l.startswith(('#', 'track', 'browser'))} - set(L)
                if bad: err.append(f'{name}: {t["url"]} contigs not in reference: {sorted(bad)[:3]}')
    # every referenced file must be committed, or GitHub Pages will not serve it (e.g. *.bam / *.bai in .gitignore)
    refs = set()
    for _, s in stops(T):
        refs |= {s['ref'][k] for k in ('fastaURL', 'indexURL', 'compressedIndexURL')}
        for t in [u for x in s['tracks'] for u in (x['tracks'] if x.get('type') == 'merged' else [x])]:
            refs.add(t['url'])
            if t.get('format') == 'bam': refs.add(t.get('indexURL', t['url'] + '.bai'))
    g = subprocess.run(['git', 'ls-files', '--others', '--ignored', '--exclude-standard', '--', os.path.join(EV, 'tour')],
                       cwd=REPO, capture_output=True, text=True)
    if g.returncode == 0:
        ign = {os.path.relpath(os.path.join(REPO, x), EV) for x in g.stdout.split()}
        bad = sorted(refs & ign)
        if bad: err.append(f'{len(bad)} referenced files are git-ignored (not deployed), e.g. {bad[0]}')
    else: warn.append('git not available: deploy (git-ignore) check skipped')
    for root, _, fs in os.walk(os.path.join(EV, 'tour')):
        for x in fs:
            sz = os.path.getsize(os.path.join(root, x))
            if sz > 100e6: err.append(f'{x}: {sz / 1e6:.1f} MB > 100 MB')
            elif sz > 50e6: warn.append(f'{x}: {sz / 1e6:.1f} MB > 50 MB (GitHub warns)')
    print(f'tour stops checked: {n}; errors: {len(err)}; warnings: {len(warn)}')
    for e in err + warn: print(' ', e)
    sys.exit(1 if err else 0)


if __name__ == '__main__':
    main()
