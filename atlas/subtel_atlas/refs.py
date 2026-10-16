"""Reference data for queries: the 48 CHM13 v2.0 500 kb ends and the catalogue v0.1 tables used for assignment.

fetch_refs() builds a reference directory:
  chm13_ends500k.fa   outer 500 kb of each CHM13 v2.0 chromosome end (chr1p .. chrYq; genome orientation: p ends start at the
                      telomere, q ends finish at it), upper case. Extracted by HTTP range reads from the public uncompressed CHM13 v2.0
                      FASTA and its .fai (the source recorded in data/catalogue/PROVENANCE.json), or from a local copy (--chm13).
  alleles_structural.tsv, haplogroups.tsv, allele_frequencies_by_superpop.tsv, representatives.fa.gz, PROVENANCE.json
                      catalogue release files, copied from a repository checkout or downloaded from the repository's raw URL.
Every file is SHA-256-checked against data/reference_checksums.json shipped with the package (one sequence checksum per CHM13 end).
"""
import os, json, gzip, time, shutil, urllib.request
from importlib import resources
from .common import CHM13_ENDS, sha256_str, sha256_file, readfa
from . import CATALOGUE_VERSION, __version__

CHM13_URL = 'https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/chm13v2.0.fa'
RAW_URL = 'https://raw.githubusercontent.com/jlanej/human-subtelomere-atlas/main/data/catalogue/'
END_LEN = 500_000
CAT_FILES = ['alleles_structural.tsv', 'haplogroups.tsv', 'allele_frequencies_by_superpop.tsv', 'representatives.fa.gz', 'paralog_families.tsv',
             'PROVENANCE.json']
CHM13_FA = 'chm13_ends500k.fa'
MANIFEST = 'refs_manifest.json'


def expected():
    return json.loads(resources.files('subtel_atlas').joinpath('data/reference_checksums.json').read_text())


def default_dir():
    return os.environ.get('SUBTEL_ATLAS_REFS') or os.path.join(os.path.expanduser('~'), '.cache', 'subtel-atlas', CATALOGUE_VERSION)


def repo_catalogue_dir():
    """data/catalogue of the repository checkout this package was installed from (editable/source installs), if present."""
    d = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'catalogue'))
    return d if os.path.exists(os.path.join(d, 'SHA256SUMS')) else None


def _get(url, rng=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={'Range': f'bytes={rng[0]}-{rng[1]}'} if rng else {})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if rng and len(data) != rng[1] - rng[0] + 1: raise IOError(f'short range read {len(data)}')
            return data
        except Exception:
            if i == tries - 1: raise
            time.sleep(2 * (i + 1))


def _fai_offset(rec, pos):
    _, _, off, lb, lw = rec
    return off + (pos // lb) * lw + pos % lb


def chm13_ends_remote(url=CHM13_URL, log=print):
    fai = {}
    for line in _get(url + '.fai').decode().splitlines():
        f = line.split('\t'); fai[f[0]] = (f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    out = {}
    for end in CHM13_ENDS:
        c = end[:-1]; rec = fai[c]; L = rec[1]
        a, b = (0, END_LEN) if end.endswith('p') else (L - END_LEN, L)
        raw = _get(url, (_fai_offset(rec, a), _fai_offset(rec, b - 1)))
        out[end] = raw.decode().replace('\n', '').replace('\r', '').upper()
        log(f'  {end}: {len(out[end]):,} bp')
    return out


def chm13_ends_local(path, log=print):
    want = {f'chr{c}' for c in list(range(1, 23)) + ['X', 'Y']}; out = {}
    for n, s in readfa(path):
        if n in want:
            s = s.upper(); out[f'{n}p'] = s[:END_LEN]; out[f'{n}q'] = s[len(s) - END_LEN:]
    return out


def write_fa(fn, seqs, order):
    tmp = fn + '.tmp'
    with open(tmp, 'w') as f:
        for n in order:
            s = seqs[n]; f.write(f'>{n}\n')
            for i in range(0, len(s), 80): f.write(s[i:i + 80] + '\n')
    os.replace(tmp, fn)


def fetch_refs(outdir=None, chm13=None, catalogue=None, chm13_url=CHM13_URL, raw_url=RAW_URL, force=False, log=print):
    """Create/verify a reference directory. chm13: local CHM13 v2.0 FASTA (else range reads from chm13_url); catalogue: local
    data/catalogue directory (else the repository checkout, else raw_url). Returns the manifest dict."""
    outdir = outdir or default_dir(); os.makedirs(outdir, exist_ok=True); X = expected()
    if X['catalogue_version'] != CATALOGUE_VERSION: raise SystemExit('package checksum table does not match its catalogue version')
    fa = os.path.join(outdir, CHM13_FA)
    have = verify(outdir, quiet=True) if not force else {}
    if not have.get(CHM13_FA):
        log(f'CHM13 v2.0 ends: {"local " + chm13 if chm13 else "range reads from " + chm13_url}')
        seqs = chm13_ends_local(chm13, log) if chm13 else chm13_ends_remote(chm13_url, log)
        bad = [e for e in CHM13_ENDS if sha256_str(seqs.get(e, '')) != X['chm13_ends'][e]]
        if bad: raise SystemExit(f'CHM13 end checksum mismatch: {",".join(bad)} (is the source CHM13 v2.0?)')
        write_fa(fa, seqs, CHM13_ENDS)
    src = catalogue or repo_catalogue_dir()
    for f in CAT_FILES:
        dst = os.path.join(outdir, f)
        if have.get(f): continue
        if src and os.path.exists(os.path.join(src, f)):
            shutil.copyfile(os.path.join(src, f), dst + '.tmp'); log(f'  {f}: copied from {src}')
        else:
            open(dst + '.tmp', 'wb').write(_get(raw_url + f)); log(f'  {f}: downloaded from {raw_url}')
        if sha256_file(dst + '.tmp') != X['catalogue_files'][f]:
            raise SystemExit(f'{f}: checksum mismatch against catalogue {CATALOGUE_VERSION}')
        os.replace(dst + '.tmp', dst)
    man = dict(catalogue_version=CATALOGUE_VERSION, tool_version=__version__, chm13_source=chm13 and os.path.basename(chm13) or chm13_url,
               files={f: sha256_file(os.path.join(outdir, f)) for f in [CHM13_FA] + CAT_FILES})
    json.dump(man, open(os.path.join(outdir, MANIFEST), 'w'), indent=1, sort_keys=True)
    ok = verify(outdir)
    if not all(ok.values()): raise SystemExit('reference verification failed')
    return man


def verify(outdir, quiet=False):
    """{file: ok} for every reference file (CHM13 ends checked per sequence)."""
    X = expected(); res = {}
    fa = os.path.join(outdir, CHM13_FA)
    if os.path.exists(fa):
        seqs = dict(readfa(fa))
        res[CHM13_FA] = list(seqs) == CHM13_ENDS and all(sha256_str(seqs[e]) == X['chm13_ends'][e] for e in CHM13_ENDS)
    else: res[CHM13_FA] = False
    for f in CAT_FILES:
        p = os.path.join(outdir, f); res[f] = os.path.exists(p) and sha256_file(p) == X['catalogue_files'][f]
    if not quiet:
        for f, ok in res.items(): print(f'{"OK " if ok else "BAD"} {f}')
    return res
