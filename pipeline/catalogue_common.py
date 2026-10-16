"""Shared paths, constants and helpers for the allele-catalogue pipeline (pipeline/catalogue_*.py).
Input locations (see pipeline/config.sh and pipeline/INPUTS.md):
  * ATLAS_INPUT_DIR set   -> the original analysis home holding the Release 2 tip tables (tipdup/), prior subtelomere results
                             (results/, subtel/), CHM13 (ref/) and HG002 files (hg002/);
  * ATLAS_INPUT_DIR unset -> public mode: the committed snapshot of the small derived inputs (data/inputs/<same relative path>[.gz])
                             and files downloaded from public URLs by pipeline/fetch_public_inputs.sh (work/public/<same path>).
    Optional caches of the analysis home (earlier 100 kb windows, the 20 kb tip FASTA) are simply absent and are re-fetched or
    skipped. ATLAS_SUBSET_ENDS (file with one end_id per line) restricts steps 01 and 03 to those ends (verification runs).
All outputs are repo-relative: large intermediates under work/, release tables under data/catalogue/.
The sequence and arm helpers, the telomere caller, the 5 kb bin method, the mosaic, structure and haplogroup functions live in the
query-tool package (atlas/subtel_atlas) and are imported from there, so the catalogue and the query tool share one implementation.
"""
import os, sys, gzip, hashlib  # noqa: F401  (gzip/hashlib kept for scripts that used them via this module)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'atlas'))          # package source (works without installing it)
from subtel_atlas.common import (WINDOW, rc, readfa, sha256_str, sha256_file, arm_of, ARM_GROUPS, group_of, group_ends,  # noqa: E402,F401
                                 DEFINING_CONFIDENCE)

INPUT = os.environ.get('ATLAS_INPUT_DIR') or None
PUBLIC_DIRS = (os.path.join(REPO, 'data', 'inputs'), os.path.join(REPO, 'work', 'public'))
SUBSET = os.environ.get('ATLAS_SUBSET_ENDS') or None
THREADS = int(os.environ.get('ATLAS_THREADS', '3'))
S3 = os.environ.get('ATLAS_S3_HTTPS', 'https://human-pangenomics.s3.amazonaws.com')
WORK = os.path.join(REPO, 'work')
WIN = os.path.join(WORK, 'windows')          # oriented 100 kb end windows (git-ignored)
CW = os.path.join(WORK, 'catalogue')         # catalogue intermediates (git-ignored)
CAT = os.path.join(REPO, 'data', 'catalogue')
FIGQC = os.path.join(REPO, 'figures', 'qc')
HG002_URL = f'{S3}/T2T/HG002/assemblies/hg002v1.1.fasta.gz'
CHM13_URL = f'{S3}/T2T/CHM13/assemblies/chm13v2.0.fa'   # range-readable public copy (windows are read from a local copy)

def need_input():
    if not INPUT:
        raise SystemExit('ATLAS_INPUT_DIR is not set (source pipeline/config.sh)')
    return INPUT

def inp(*p):
    """Path of an input file given by its analysis-home relative path. Public mode: data/inputs/<p>, data/inputs/<p>.gz,
    work/public/<p>, work/public/<p>.gz (first that exists; otherwise the work/public path, so that existence checks fail cleanly)."""
    if INPUT: return os.path.join(INPUT, *p)
    for d in PUBLIC_DIRS:
        for c in (os.path.join(d, *p), os.path.join(d, *p) + '.gz'):
            if os.path.exists(c): return c
    return os.path.join(PUBLIC_DIRS[1], *p)

def open_text(fn):
    """open a text input that may be gzip-compressed (public snapshot)"""
    return gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn)

def subset_ids():
    """end_ids of ATLAS_SUBSET_ENDS, or None"""
    return set(l.strip() for l in open(SUBSET) if l.strip()) if SUBSET else None

def write_tsv_atomic(df, fn, **kw):
    tmp = fn + '.tmp'; df.to_csv(tmp, sep='\t', index=False, **kw); os.replace(tmp, fn)
