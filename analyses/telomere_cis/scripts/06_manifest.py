"""Checksums and provenance for the telomere_cis track.
Outputs: MANIFEST.tsv (path relative to the track, bytes, md5, sha256 of every file in data/, results/, figures/, scripts/),
         PROVENANCE.json (input files relative to ATLAS_INPUT_DIR with sha256, public source URLs, literature, run date).
"""
import os, sys, json, hashlib, datetime, platform, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tl_common import TRACK, REPO, inp, sha256

def md5(fn):
    h = hashlib.md5()
    with open(fn, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
rows = []
for sub in ('data', 'results', 'figures', 'scripts'):
    for root, _, files in os.walk(os.path.join(TRACK, sub)):
        for f in sorted(files):
            if f.endswith('.pyc'): continue
            p = os.path.join(root, f); rows.append((os.path.relpath(p, TRACK), os.path.getsize(p), md5(p), sha256(p)))
for f in ('README.md', 'run.sh'):
    p = os.path.join(TRACK, f)
    if os.path.exists(p): rows.append((f, os.path.getsize(p), md5(p), sha256(p)))
with open(os.path.join(TRACK, 'MANIFEST.tsv'), 'w') as fo:
    fo.write('path\tbytes\tmd5\tsha256\n')
    for r in sorted(rows): fo.write('\t'.join(map(str, r)) + '\n')
inputs = ['results/subtel_foreign_tips.tsv', 'results/subtel_hg002_ends.tsv', 'results/single_read_zone_by_tip.tsv', 'meta/explorer_seqdata.csv',
          'subtel/foreign/ref/chm13_ends500k.fa', 'hg002/ref_ends.fa', 'hg002/hg002v1.1.fasta.gz.fai',
          'hg002/bams/ont_ul_dorado.bam', 'hg002/bams/revio_3cells.bam', 'hg002/bams/hifi_dcv1.1.bam']
def tool_version(cmd):
    try: return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip().splitlines()[0]
    except Exception: return None
import numpy, pandas, scipy, pysam, matplotlib
prov = dict(
    track='analyses/telomere_cis', run_date_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d'),
    inputs_relative_to_ATLAS_INPUT_DIR={f: (sha256(inp(f)) if os.path.exists(inp(f)) else None) for f in inputs},
    work_intermediates=['work/windows/fetch/<asm>.fa (100 kb tip windows, pipeline/catalogue_02_fetch_windows.sh)',
                        'work/telomere_cis/hg002_ends40k.fa (outer 40 kb of 92 HG002 v1.1 ends, samtools faidx range read)',
                        'work/telomere_cis/fixed_windows.paf, hg002_fixed.paf, hg002_reads.tsv.gz'],
    public_sources={
        'HPRC Release 2 assemblies': 's3://human-pangenomics/working/HPRC/<sample>/assemblies/release2/ (index: HPRC release 2 assembly index)',
        'HG002 Q100 v1.1 assembly': 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/hg002v1.1.fasta.gz',
        'HG002 read alignments to v1.1': 'https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/ '
                                         '(ont_r10_ul_dorado/hg002v1.1_ont_r10_ul_dorado.pri.bam; hifi_revio_3cell/hg002v1.1_hifi_revio_3cells.bam; '
                                         'hifi_dcv1.1/hg002v1.1_hifi_dcv1.1.pri.bam), primary alignments overlapping the outer 20 kb of each end',
        'CHM13 v2.0': 'https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz (outer 500 kb of 48 ends)',
        'HPRC sequencing metadata': 'HPRC data explorer sequencing table (explorer_seqdata.csv)'},
    literature={'Karimian et al. 2024 Science 384:533-539': 'doi:10.1126/science.ado0431',
                'Karimian et al. bioRxiv 2023 (figure source, v1 PDF)': 'doi:10.1101/2023.12.21.572870; '
                'https://www.biorxiv.org/content/biorxiv/early/2024/01/13/2023.12.21.572870.full.pdf',
                'figure source sha256': sha256(os.path.join(REPO, 'work', 'telomere_cis', 'karimian2023_biorxiv_v1.pdf'))
                if os.path.exists(os.path.join(REPO, 'work', 'telomere_cis', 'karimian2023_biorxiv_v1.pdf')) else None},
    software=dict(python=platform.python_version(), numpy=numpy.__version__, pandas=pandas.__version__, scipy=scipy.__version__,
                  pysam=pysam.__version__, matplotlib=matplotlib.__version__, minimap2=tool_version('minimap2 --version'),
                  samtools=tool_version('samtools --version')),
    random_seed=20261015)
json.dump(prov, open(os.path.join(TRACK, 'PROVENANCE.json'), 'w'), indent=1)
print('manifest entries', len(rows))
