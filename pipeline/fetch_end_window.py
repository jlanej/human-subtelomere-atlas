"""Fetch the tip-frame window of any catalogue end from its public source and verify it against the released SHA-256.
Usage: python pipeline/fetch_end_window.py END_ID [END_ID ...] > windows.fa
Reads data/catalogue/ends.tsv (assembly_url, contig, win_start, win_end, strand, window_sha256). Release 2, CHM13 and HG002 Q100
sources are bgzipped FASTA with remote .fai/.gzi and are range-read with `samtools faidx` (htslib with libcurl). The window is
reverse-complemented for R ends ('-' strand). Pedigree (pedigree_aux) sources are plain FASTA and are not supported here."""
import os, sys, subprocess, hashlib, pandas as pd
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S3 = os.environ.get('ATLAS_S3_HTTPS', 'https://human-pangenomics.s3.amazonaws.com')
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
E = pd.read_csv(os.path.join(REPO, 'data', 'catalogue', 'ends.tsv'), sep='\t', low_memory=False).set_index('end_id')
bad = 0
for e in sys.argv[1:]:
    r = E.loc[e]
    if r['set'] == 'pedigree_aux': sys.stderr.write(f'{e}: pedigree sources not supported\n'); bad += 1; continue
    url = r.assembly_url.replace('s3://human-pangenomics', S3)
    out = subprocess.run(['samtools', 'faidx', url, f'{r.contig}:{int(r.win_start)}-{int(r.win_end)}'], capture_output=True, text=True, check=True).stdout
    seq = ''.join(out.splitlines()[1:]).upper()
    if r.strand == '-': seq = seq.translate(comp)[::-1]
    ok = hashlib.sha256(seq.encode()).hexdigest() == r.window_sha256; bad += not ok
    sys.stderr.write(f'{e}\t{len(seq)} bp\tsha256 {"OK" if ok else "MISMATCH"}\n')
    print(f'>{e}\n' + '\n'.join(seq[i:i + 100] for i in range(0, len(seq), 100)))
sys.exit(1 if bad else 0)
