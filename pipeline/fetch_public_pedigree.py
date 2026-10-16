"""Range-read the outer 250 kb of both ends of every chromosome of the WashU pedigree v1.1 assemblies (PAN010, PAN011, PAN027,
PAN028; public plain FASTA on the UCSC server) into work/public/subtel/pedigree/ends250k.fa (7 haplotypes) and
ends250k_PAN028_h2.fa, byte-identical to the files used by the catalogue (re-implementation of the analysis-home script
subtel_ped_fetch_ends.py; sizes and SHA-256 in data/inputs/public_files.tsv). Byte offsets come from the committed index copies
data/inputs/subtel/pedigree/<hap>.v1.1.fai. Header '>HAP|chrN|L|chrom_len' / '>HAP|chrN|R|chrom_len', sequence upper case on one
line, assembly orientation (R windows are not reverse-complemented here; step 03 orients them).
Usage: python pipeline/fetch_public_pedigree.py   (ATLAS_FETCH_JOBS haplotypes in parallel, default 4)"""
import os, sys, time, urllib.request, concurrent.futures as cf
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B = 'https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/'
FILES = {'PAN027_pat': 'PAN027_pat_HiFi_element_final_pat.polished.fasta', 'PAN027_mat': 'PAN027_mat_HiFi_element_final_mat.polished.fasta',
         'PAN010_h1': 'PAN010_hap1_HiFi_element_final_hap1.polished.fasta', 'PAN010_h2': 'PAN010_hap2_HiFi_element_final_hap2.polished.fasta',
         'PAN011_h1': 'PAN011_hap1_HiFi_element_final_XY_hap1.polished.fasta', 'PAN011_h2': 'PAN011_hap2_HiFi_element_final_XY_hap2.polished.fasta',
         'PAN028_h1': 'PAN028_hap1_HiFi_element_final_hap1.polished.fasta', 'PAN028_h2': 'PAN028_hap2_HiFi_element_final_hap2.polished.fasta'}
FL = 250_000
FAI = os.path.join(REPO, 'data', 'inputs', 'subtel', 'pedigree'); OUT = os.path.join(REPO, 'work', 'public', 'subtel', 'pedigree')

def rng(url, a, b):
    for t in range(4):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={'Range': f'bytes={a}-{b}'}), timeout=300).read()
        except Exception:
            if t == 3: raise
            time.sleep(10 * (t + 1))

def fetch(url, off, lb, lw, s, e):  # 0-based half-open sequence coordinates -> sequence (fai line layout)
    a = off + (s // lb) * lw + s % lb; b = off + ((e - 1) // lb) * lw + (e - 1) % lb
    seq = rng(url, a, b).replace(b'\n', b'').replace(b'\r', b'').decode().upper()
    assert len(seq) == e - s, (url, len(seq), e - s)
    return seq

def hap_records(hap):
    recs = []
    for line in open(os.path.join(FAI, f'{hap}.v1.1.fai')):
        name, L, off, lb, lw = line.rstrip('\n').split('\t')[:5]; L, off, lb, lw = int(L), int(off), int(lb), int(lw)
        chrom = name.split('.')[1]
        for end, (s, e) in (('L', (0, FL)), ('R', (L - FL, L))):
            recs.append(f'>{hap}|{chrom}|{end}|{L}\n{fetch(B + FILES[hap], off, lb, lw, s, e)}\n')
    print(hap, 'done', len(recs), 'ends', flush=True); return ''.join(recs)

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    groups = {'ends250k.fa': [h for h in FILES if h != 'PAN028_h2'], 'ends250k_PAN028_h2.fa': ['PAN028_h2']}
    with cf.ThreadPoolExecutor(int(os.environ.get('ATLAS_FETCH_JOBS', '4'))) as ex:
        txt = dict(zip(FILES, ex.map(hap_records, list(FILES))))
    for fn, haps in groups.items():
        p = os.path.join(OUT, fn); open(p + '.tmp', 'w').write(''.join(txt[h] for h in haps)); os.replace(p + '.tmp', p)
