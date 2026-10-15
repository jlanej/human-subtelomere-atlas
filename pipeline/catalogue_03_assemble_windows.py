"""Step 1c: assemble oriented telomere-anchored windows for all ends and write the per-end window table.
Frame: every window is oriented telomere -> centromere (L ends as on the contig, R ends reverse-complemented), upper case.
Sources:
  Release 2 tips : tipdup/tips_meta.tsv + windows from subtel/foreign/r2_100k (reused) or work/windows/fetch (S3 range reads);
                   junk_len / telo_len taken from the tip table (same telomere caller as catalogue_telomere.py).
  CHM13 v2.0     : ref/chm13v2.0.fa.gz (local copy of the T2T analysis set), 100 kb at both ends of chr1-22, X, Y (48 ends).
  HG002 Q100 v1.1: work/windows/hg002.fa (S3 range reads of the outer 100 kb of the 46 chromosomes, 92 ends).
  Pedigree (auxiliary, validation only; excluded from allele definition and frequencies): WashU PAN010/PAN011/PAN027/PAN028
                   v1.1 assemblies, outer 250 kb of every chromosome end previously range-read into subtel/pedigree/ends250k*.fa
                   (assembly orientation); the outer 100 kb in tip frame is used here.
Sanity check: for Release 2 tips, the first 20 kb of each window must equal the earlier tipdup/tips_masked.fa tip window at all
unmasked positions (that file masks the telomere with N).
Outputs: work/windows/ends.fa(.fai) and work/windows/ends_windows.tsv
Usage  : python catalogue_03_assemble_windows.py [hg002-regions|build]
"""
import os, sys, re, subprocess, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import catalogue_common as C, catalogue_telomere as TL

def chrom_regions(fai, keep):
    regs = []
    for line in open(fai):
        c, L = line.split('\t')[:2]; L = int(L)
        if not keep(c): continue
        w = min(C.WINDOW, L)
        regs += [(c, 'L', 1, w, L), (c, 'R', L - w + 1, L, L)]
    return regs
def hg002_regs(): return chrom_regions(C.inp('hg002/hg002v1.1.fasta.gz.fai'), lambda c: re.match(r'chr[0-9XY]+_(MATERNAL|PATERNAL)$', c))
def chm13_regs(): return chrom_regions(C.inp('ref/chm13v2.0.fa.gz.fai'), lambda c: re.match(r'chr[0-9XY]+$', c))

if len(sys.argv) > 1 and sys.argv[1] == 'hg002-regions':
    open(f'{C.WIN}/hg002.regions', 'w').write(''.join(f'{c}:{a}-{b}\n' for c, e, a, b, L in hg002_regs())); sys.exit()

def load(fns):
    d = {}
    for fn in fns:
        for n, s in C.readfa(fn): d[n] = s.upper()
    return d

rows = []; out = open(f'{C.WIN}/ends.fa.tmp', 'w')
def emit(rec, seq):
    rec['window_len'] = len(seq); rec['window_sha256'] = C.sha256_str(seq); rows.append(rec)
    out.write(f">{rec['end_id']}\n"); [out.write(seq[i:i + 100] + '\n') for i in range(0, len(seq), 100)]

# ---- Release 2
M = pd.read_csv(C.inp('tipdup/tips_meta.tsv'), sep='\t').set_index('q')
P = pd.read_csv(f'{C.WIN}/r2_plan.tsv', sep='\t')
E = pd.read_csv(C.inp('meta/explorer_assemblies.csv'))
E['base'] = E.assembly.astype(str).map(os.path.basename)
acc = dict(zip(E.base, E.genbank_accession)); aname = dict(zip(E.base, E.assembly_name))
short = {n: s.upper() for n, s in C.readfa(C.inp('tipdup/tips_masked.fa'))}
nbad = nchk = 0
for asm, g in P.groupby('asm', sort=True):
    seqs = load([f for f in (C.inp(f'subtel/foreign/r2_100k/{asm}.fa'), f'{C.WIN}/fetch/{asm}.fa') if os.path.exists(f)])
    for r in g.itertuples():
        if r.region not in seqs: raise SystemExit(f'missing window {r.region} ({asm}); run catalogue_02_fetch_windows.sh')
        s = seqs[r.region]; s = s if r.side == 'L' else C.rc(s)
        m = M.loc[r.q]; base = os.path.basename(r.s3_path)
        if r.q in short:
            t = short[r.q]; n = min(len(t), len(s)); nchk += 1
            if any(a != b for a, b in zip(s[:n], t[:n]) if b != 'N'): nbad += 1
        sample, hapn, _ = r.contig.split('#', 2)
        emit(dict(end_id=r.end_id, set='HPRC_R2', sample=sample, haplotype=asm.split('.')[1], pansn_hap=int(hapn),
                  assembly=asm, assembly_name=aname.get(base, base.replace('.fa.gz', '')), assembly_accession=acc.get(base, ''),
                  assembly_url=r.s3_path, assembler=m.assembler, hifi_input=m.hifi_input, contig=r.contig, contig_len=r.contig_len,
                  side=r.side, win_start=r.win_start, win_end=r.win_end, strand='+' if r.side == 'L' else '-',
                  contig_arm=m.arm if isinstance(m.arm, str) else '', junk_len=int(m.junk_len), telo_len=int(m.telo_len),
                  telo_source='tip_table', window_source=r.source), s)
print('R2 ends', len(rows), '| 20 kb sanity check: checked', nchk, 'mismatching', nbad)

# ---- CHM13 v2.0 and HG002 Q100 v1.1
def add_t2t(regs, seqs, setname, sample, asmname, accn, url, hapf):
    for c, e, a, b, L in regs:
        s = seqs[f'{c}:{a}-{b}']; s = s if e == 'L' else C.rc(s)
        tb = TL.find_block(s) or dict(junk_len=0, telo_len=0)
        chrom = c.split('_')[0]; hap, hn = hapf(c)
        emit(dict(end_id=f'{sample}#{hn}#{c}#{e}', set=setname, sample=sample, haplotype=hap, pansn_hap=hn, assembly=asmname,
                  assembly_name=asmname, assembly_accession=accn, assembly_url=url, assembler='T2T', hifi_input='', contig=c,
                  contig_len=L, side=e, win_start=a, win_end=b, strand='+' if e == 'L' else '-',
                  contig_arm=chrom[3:] + ('p' if e == 'L' else 'q'), junk_len=tb['junk_len'], telo_len=tb['telo_len'],
                  telo_source='catalogue_telomere', window_source='local' if setname == 'CHM13' else 'fetch'), s)
cr = chm13_regs()
reg = f'{C.WIN}/chm13.regions'; open(reg, 'w').write(''.join(f'{c}:{a}-{b}\n' for c, e, a, b, L in cr))
subprocess.run(f"samtools faidx -r {reg} {C.inp('ref/chm13v2.0.fa.gz')} > {C.WIN}/chm13.fa", shell=True, check=True)
add_t2t(cr, load([f'{C.WIN}/chm13.fa']), 'CHM13', 'CHM13', 'chm13v2.0', 'GCA_009914755.4', C.CHM13_URL, lambda c: ('haploid', 0))
add_t2t(hg002_regs(), load([f'{C.WIN}/hg002.fa']), 'HG002_Q100', 'HG002', 'hg002v1.1', '', C.HG002_URL,
        lambda c: ('pat', 1) if c.endswith('PATERNAL') else ('mat', 2))
PED_URL = 'https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/'
PED_FILES = {'PAN027_pat': 'PAN027_pat_HiFi_element_final_pat.polished.fasta', 'PAN027_mat': 'PAN027_mat_HiFi_element_final_mat.polished.fasta',
             'PAN010_h1': 'PAN010_hap1_HiFi_element_final_hap1.polished.fasta', 'PAN010_h2': 'PAN010_hap2_HiFi_element_final_hap2.polished.fasta',
             'PAN011_h1': 'PAN011_hap1_HiFi_element_final_XY_hap1.polished.fasta', 'PAN011_h2': 'PAN011_hap2_HiFi_element_final_XY_hap2.polished.fasta',
             'PAN028_h1': 'PAN028_hap1_HiFi_element_final_hap1.polished.fasta', 'PAN028_h2': 'PAN028_hap2_HiFi_element_final_hap2.polished.fasta'}
for fn in (C.inp('subtel/pedigree/ends250k.fa'), C.inp('subtel/pedigree/ends250k_PAN028_h2.fa')):
    for n, s in C.readfa(fn):
        hap, chrom, e, L = n.split('|'); L = int(L); s = s.upper(); fl = len(s)
        s = s[:C.WINDOW] if e == 'L' else C.rc(s)[:C.WINDOW]
        a, b = (1, len(s)) if e == 'L' else (L - len(s) + 1, L)
        tb = TL.find_block(s) or dict(junk_len=0, telo_len=0); sample, h = hap.split('_')
        emit(dict(end_id=f'{sample}#{h}#{chrom}#{e}', set='pedigree_aux', sample=sample, haplotype=h, pansn_hap=0, assembly=f'{hap}.v1.1',
                  assembly_name=PED_FILES[hap], assembly_accession='', assembly_url=PED_URL + PED_FILES[hap], assembler='WashU_v1.1',
                  hifi_input='', contig=chrom, contig_len=L, side=e, win_start=a, win_end=b, strand='+' if e == 'L' else '-',
                  contig_arm=chrom[3:] + ('p' if e == 'L' else 'q'), junk_len=tb['junk_len'], telo_len=tb['telo_len'],
                  telo_source='catalogue_telomere', window_source='subtel_pedigree_ends250k'), s)
out.close(); os.replace(f'{C.WIN}/ends.fa.tmp', f'{C.WIN}/ends.fa')
subprocess.run(f'samtools faidx {C.WIN}/ends.fa', shell=True, check=True)
D = pd.DataFrame(rows); assert D.end_id.is_unique
C.write_tsv_atomic(D, f'{C.WIN}/ends_windows.tsv')
print(D.groupby('set').agg(n=('end_id', 'size'), full=('window_len', lambda x: (x == C.WINDOW).sum())).to_string())
