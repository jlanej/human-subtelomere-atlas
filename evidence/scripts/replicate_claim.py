#!/usr/bin/env python3
"""Re-derive one read-evidence verdict from public sources only (no file produced by the authors is needed beyond this repository).

  python evidence/scripts/replicate_claim.py CLAIM_ID [--threads 4] [--chm13 chm13v2.0.fa.gz] [--genome chm13v2.map-ont.mmi]
  (or evidence/replicate_claim.sh CLAIM_ID)

CLAIM_ID is a claim_id of evidence/claims.tsv (e.g. 'HG002#2#chr9_MATERNAL#L', 'HG01258|20q|SA-20q-01') or a control id of
evidence/controls.tsv starting with 'control|' (1000 Genomes ONT strictly native controls).

Steps
  1. criteria: evidence/criteria_v1.json and the packaged copy must both hash to a89c5eaf023ed032 (thresholds frozen before the
     population analysis; evidence/CRITERIA.md).
  2. references (cached in --cache, default ~/.cache/subtel-atlas/replicate):
       CHM13 v2.0 (public analysis-set FASTA, bgzip + .gzi; ~0.98 GB download once; .fai built locally) for CRAM decoding and for the
       genome-wide anchor (minimap2 builds the map-ont index in memory on every run, ~1 min and ~8 GB RAM, unless --genome
       gives a prebuilt .mmi); the 48 CHM13 500 kb ends via `subtel-atlas fetch-refs` (SHA-256-checked against the package).
  3. claim model from the repository catalogue (data/catalogue/ends.tsv + mosaic_blocks.tsv; readevidence.claims_from_catalogue).
  4. the person's assembly end window(s), cut by range from the public assembly FASTA and checked against the catalogue
     window SHA-256 (switch concordance, criteria change C2; and own-assembly placement for the definitive tier).
  5. the cited reads (evidence/reads.tsv.gz rows of this claim), fetched by region from their public alignment file and
     filtered by read ID (samtools view -N), then evaluated with subtel_atlas.readevidence (the module behind
     `subtel-atlas reads-evidence`): stage 1 (telomere, diagnostic units, switch), one genome-wide minimap2 run for the anchor
     and integrity, stage 2 (verdict), assembly concordance.
  6. per-read verdicts and the claim verdict (CRITERIA.md section 4; HG002 phased, others unphased) are printed next to the
     published values; then the reporting-only definitive tier (evidence/scripts/evidence_lib.py).
Exit status 0 when the claim verdict, the definitive verdict and every per-read verdict match the published tables.
HG003/HG004 ONT reads are stored on GRCh38: give --grch38 (GRCh38 no-alt analysis set FASTA) to decode those CRAMs.
"""
import argparse, os, sys, json, csv, gzip, time, subprocess, tempfile, collections, urllib.request, shutil
from concurrent.futures import ThreadPoolExecutor
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'atlas')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subtel_atlas import readevidence as RE, refs  # noqa: E402
import evidence_lib as L  # noqa: E402

CRIT_SHA = 'a89c5eaf023ed032'
T0 = time.time()


def log(*a): print(f'[{time.time() - T0:6.0f} s]', *a, file=sys.stderr, flush=True)


def download(url, dst):
    if os.path.exists(dst): return dst
    os.makedirs(os.path.dirname(dst), exist_ok=True); tmp = dst + '.part'
    log(f'downloading {url}')
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, 'wb') as fh: shutil.copyfileobj(r, fh, 1 << 22)
    os.replace(tmp, dst); return dst


def read_tsv(fn):
    op = gzip.open if fn.endswith('.gz') else open
    with op(fn, 'rt') as fh: return list(csv.DictReader(fh, delimiter='\t'))


def preset_of(platform):
    return 'map-hifi' if str(platform).startswith('HiFi') else 'map-ont'


def fetch_reads(url, region, names, ref, workdir, tag):
    """primary records of `names` in `region` of the public file -> local BAM; returns (path, seconds, command)"""
    os.makedirs(workdir, exist_ok=True)
    nf = os.path.join(workdir, f'{tag}.names'); out = os.path.join(workdir, f'{tag}.bam')
    with open(nf, 'w') as fh: fh.write('\n'.join(sorted(names)) + '\n')
    cmd = ['samtools', 'view', '-b', '-F', '0x900', '-N', nf] + (['-T', ref] if ref else []) + \
          ['-X', url, L.index_url(url), region, '-o', out]
    t = time.time()
    for i in range(3):
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
        if p.returncode == 0: break
        if i == 2: raise RuntimeError(f'{url} {region}: {p.stderr[-400:]}')
        time.sleep(5)
    # provenance string: the exact command, with the local names file shown as <read_ids.txt> and no local output path
    logged = [('<read_ids.txt>' if c == nf else c) for c in cmd[:-2]]
    return out, time.time() - t, ' '.join(logged)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('claim_id')
    ap.add_argument('--cache', default=os.environ.get('SUBTEL_REPLICATE_CACHE', os.path.join(os.path.expanduser('~'), '.cache', 'subtel-atlas', 'replicate')))
    ap.add_argument('--chm13', default=os.environ.get('SUBTEL_CHM13'), help='local CHM13 v2.0 analysis-set FASTA (bgzip with .fai/.gzi); default: download')
    ap.add_argument('--genome', default=os.environ.get('SUBTEL_CHM13_MMI'), help='prebuilt minimap2 map-ont index of CHM13 v2.0 (default: the FASTA)')
    ap.add_argument('--grch38', default=os.environ.get('SUBTEL_GRCH38'), help='GRCh38 FASTA (only for HG003/HG004 ONT CRAMs)')
    ap.add_argument('--workdir', help='scratch directory (default: <cache>/runs/<claim>)')
    ap.add_argument('--no-placement', action='store_true', help='skip own-assembly placement (definitive tier)')
    ap.add_argument('-t', '--threads', type=int, default=4)
    a = ap.parse_args()
    timings = {}
    # ---------------------------------------------------------------- 1. criteria
    crit = RE.load_criteria(os.path.join(REPO, 'evidence', 'criteria_v1.json')); th = crit['thresholds']
    h_repo, h_pkg = RE.criteria_hash(crit), RE.criteria_hash(RE.load_criteria())
    log(f'criteria: evidence/criteria_v1.json {h_repo}, packaged {h_pkg}, required {CRIT_SHA}')
    if not (h_repo == h_pkg == CRIT_SHA): sys.exit('criteria hash mismatch: the frozen thresholds were modified')
    # ---------------------------------------------------------------- claim and cited reads
    is_ctl = a.claim_id.startswith('control|')
    claims = {r['claim_id']: r for r in read_tsv(os.path.join(REPO, 'evidence', 'claims.tsv'))}
    ctls = {r['control_id']: r for r in read_tsv(os.path.join(REPO, 'evidence', 'controls.tsv'))}
    if not is_ctl and a.claim_id not in claims: sys.exit(f'{a.claim_id}: not in evidence/claims.tsv')
    if is_ctl and a.claim_id not in ctls: sys.exit(f'{a.claim_id}: not in evidence/controls.tsv')
    cited = [r for r in read_tsv(os.path.join(REPO, 'evidence', 'reads.tsv.gz')) if r['claim_id'] == a.claim_id]
    if not cited: sys.exit(f'{a.claim_id}: no cited reads in evidence/reads.tsv.gz (nothing to replicate)')
    C = claims.get(a.claim_id) if not is_ctl else claims.get(ctls[a.claim_id]['allele_model']) or \
        next(c for c in claims.values() if c['claim_id'] == ctls[a.claim_id]['allele_model'])
    cohort = C['cohort']
    log(f"claim {a.claim_id} ({cohort}): {len(cited)} cited reads; published v1 {C['v1_verdict'] if not is_ctl else 'control'}")
    wd = a.workdir or os.path.join(a.cache, 'runs', a.claim_id.replace('|', '__').replace('#', '_'))
    os.makedirs(wd, exist_ok=True)
    # ---------------------------------------------------------------- 2. references
    t = time.time()
    chm13 = a.chm13 or os.path.join(a.cache, 'chm13v2.0.fa.gz')
    if not a.chm13:
        for ext in ('.gzi', ''): download(L.CHM13_FA_URL + ext, chm13 + ext)      # the server has no .fai for the bgzip file
    if not os.path.exists(chm13 + '.fai'):
        log('indexing CHM13 (samtools faidx)'); subprocess.run(['samtools', 'faidx', chm13], check=True)
    rdir = os.path.join(a.cache, 'refs')
    if not (os.path.exists(os.path.join(rdir, refs.CHM13_FA)) and all(refs.verify(rdir, quiet=True).values())):
        refs.fetch_refs(rdir, chm13=chm13, catalogue=os.path.join(REPO, 'data', 'catalogue'), log=lambda *x: log(*x))
    ends = RE.Chm13Ends(os.path.join(rdir, refs.CHM13_FA)); genome = a.genome or chm13
    timings['references_s'] = round(time.time() - t, 1)
    # ---------------------------------------------------------------- 3. claim model
    cat = os.path.join(REPO, 'data', 'catalogue')
    model_end = C['model_end']
    model = RE.claims_from_catalogue(os.path.join(cat, 'ends.tsv'), os.path.join(cat, 'mosaic_blocks.tsv'), end_ids=[model_end])[0]
    log(f"model {model_end}: D={model['D']} P={model['P']} R={model['R']} pattern={model['pattern']} switch_sub={model.get('switch_sub')}")
    mask_cache = os.path.join(a.cache, 'masks')
    # ---------------------------------------------------------------- 4. assembly windows (public, SHA-checked)
    import pandas as pd
    E = pd.read_csv(os.path.join(cat, 'ends.tsv'), sep='\t', low_memory=False); Ei = E.set_index('end_id')
    t = time.time()
    win = L.window_from_public(dict(Ei.loc[model_end], end_id=model_end), cwd=wd)
    ctx_hifi = RE.ClaimContext(model, ends, crit, preset='map-hifi', mask_cache=mask_cache)
    asm = RE.assembly_projection(ctx_hifi, win.upper())
    log(f"assembly projection of {model_end}: switch_pass={asm.get('switch_pass')} sub {asm.get('switch_ci_sub_lo')}-{asm.get('switch_ci_sub_hi')}")
    timings['assembly_window_s'] = round(time.time() - t, 1)
    # ---------------------------------------------------------------- 5. reads: fetch by region + read ID, evaluate
    t = time.time()
    groups = collections.OrderedDict()
    for r in cited: groups.setdefault((r['source_url'], r['region'], r['source_reference']), []).append(r)
    fetched = {}; cmds = []
    for i, ((url, region, sref), rows) in enumerate(groups.items()):
        ref = chm13 if sref == 'CHM13v2.0' else (a.grch38 if sref == 'GRCh38' else None)
        if url.endswith('.cram') and not ref: sys.exit(f'{url}: CRAM on {sref} needs --grch38')
        bam, s, cmd = fetch_reads(url, region, {r['read_id'] for r in rows}, ref, os.path.join(wd, 'reads'), f'g{i}')
        fetched[(url, region)] = bam; cmds.append(cmd)
        log(f'fetched {len(rows)} cited reads from {url.split("/")[-1]} {region} in {s:.0f} s')
    timings['fetch_s'] = round(time.time() - t, 1)
    t = time.time()
    ctxs = {}
    def ctx_for(preset):
        if preset not in ctxs: ctxs[preset] = RE.ClaimContext(model, ends, crit, preset=preset, mask_cache=mask_cache)
        return ctxs[preset]
    rows, Q, seqs = [], [], {}
    for (url, region, sref), cr in groups.items():
        want = {r['read_id']: r for r in cr}
        got = {}
        for name, seq, rec in RE.iter_reads(fetched[(url, region)]):
            if name in want and name not in got: got[name] = (seq, rec)
        for name, pub in want.items():
            if name not in got:
                rows.append(dict(read=name, verdict='MISSING', _pub=pub)); continue
            seq, rec = got[name]; seqs[name] = seq
            ctx = ctx_for(preset_of(pub['platform']))
            meta = dict(dataset=pub['dataset'], platform=pub['platform'])
            r, u, tip = RE.evaluate_read(ctx, name, seq, side=pub['orient_side'], meta=meta)
            r['_pub'] = pub; r['_ctx'] = ctx; k = f'{len(rows)}'
            r['_key'] = k; Q.append((f'{k}|w', tip))
            if len(tip) - r['anchor_from'] >= th['anchor_min_bp']: Q.append((f'{k}|a', tip[r['anchor_from']:]))
            rows.append(r)
    timings['stage1_s'] = round(time.time() - t, 1); t = time.time()
    log(f'genome-wide anchor mapping of {len(Q)} segments ({os.path.basename(genome)})')
    G = RE.map_genome(Q, genome, a.threads) if Q else {}
    timings['genome_mapping_s'] = round(time.time() - t, 1)
    for r in rows:
        if r.get('verdict') == 'MISSING': continue
        RE.finalize(r['_ctx'], r, G.get(f"{r['_key']}|w"), G.get(f"{r['_key']}|a")); RE.concordance(r, asm, th)
    # ---------------------------------------------------------------- 6. compare
    ok_reads = 0; out = []
    print(f"\n{'read':40s} {'dataset':28s} {'published':16s} {'replicated':16s} {'conc pub/rep':12s} switch_ci_sub (rep)")
    for r in rows:
        p = r['_pub']
        same = r['verdict'] == p['v1_verdict'] and (r['verdict'] != 'conclusive' or bool(r.get('switch_concordant')) == L.b(p['switch_concordant']))
        ok_reads += same
        print(f"{r['read'][:40]:40s} {str(p['dataset'])[:28]:28s} {p['v1_verdict']:16s} {r['verdict']:16s} "
              f"{str(L.b(p['switch_concordant']))[0]}/{str(bool(r.get('switch_concordant')))[0]:10s} "
              f"{r.get('switch_ci_sub_lo')}-{r.get('switch_ci_sub_hi')}{'' if same else '   <-- differs'}")
        out.append(dict(read_id=r['read'], dataset=p['dataset'], published_verdict=p['v1_verdict'], replicated_verdict=r['verdict'],
                        published_switch_concordant=p['switch_concordant'], replicated_switch_concordant=r.get('switch_concordant'),
                        switch_ci_sub=f"{r.get('switch_ci_sub_lo')}-{r.get('switch_ci_sub_hi')}", anchor_chrom=r.get('anchor_chrom'),
                        anchor_mapq=r.get('anchor_mapq'), anchor_unique_bp=r.get('anchor_unique_bp'), match=same))
    counted = [r for r in rows if r['verdict'] != 'MISSING' and L.b(r['_pub']['counted'])]
    for r in counted: r['contig'] = 'claimed'      # HG002 cited carrier reads are all on the claimed contig (phased set)
    v = RE.claim_verdict(counted, crit, same_hap_key='contig' if cohort == 'HG002_Q100' and not is_ctl else None)
    res = dict(claim_id=a.claim_id, cohort=cohort, criteria_sha=h_repo, n_cited=len(cited), n_reads_match=ok_reads,
               n_missing=sum(r['verdict'] == 'MISSING' for r in rows), replicated=v, timings=timings,
               minimap2=subprocess.run(['minimap2', '--version'], capture_output=True, text=True).stdout.strip(),
               samtools=subprocess.run(['samtools', '--version'], capture_output=True, text=True).stdout.split('\n')[0],
               fetch_commands=cmds)
    if is_ctl:
        k = ctls[a.claim_id]; res['published'] = dict(n_conclusive=int(k['n_conclusive']), n_contradicts=int(k['n_contradicts']))
        n_conc = sum(r['verdict'] == 'conclusive' for r in rows)
        res['match'] = ok_reads == len(cited) and n_conc == int(k['n_conclusive'])
        print(f"\ncontrol {a.claim_id}: conclusive reads replicated {n_conc} (published {k['n_conclusive']}); "
              f"contradicting {sum(r['verdict'] == 'contradicts' for r in rows)} (published {k['n_contradicts']} of all control reads)")
    else:
        res['published'] = dict(v1_verdict=C['v1_verdict'], n_conclusive_concordant=int(C['n_conclusive_concordant']),
                                n_contradicts=int(C['n_contradicts']), n_consistent=int(C['n_consistent']), definitive_verdict=C['definitive_verdict'])
        print(f"\nclaim {a.claim_id}: v1 verdict replicated {v['verdict']} (published {C['v1_verdict']}); conclusive concordant "
              f"{v['n_conclusive_concordant']} ({C['n_conclusive_concordant']}), contradicting {v['n_contradicts']} ({C['n_contradicts']}), "
              f"consistent {v['n_consistent']} ({C['n_consistent']})")
        res['match'] = (ok_reads == len(cited) and v['verdict'] == C['v1_verdict'] and
                        v['n_conclusive_concordant'] == int(C['n_conclusive_concordant']))
        # ------------------------------------------------------------ definitive tier
        if not a.no_placement:
            t = time.time()
            person = C['individual']; pset = 'HG002_Q100' if cohort in ('HG002_Q100', 'HG002_trio') else 'HPRC_R2'
            ids = E[(E['set'] == pset) & ((E['sample'] == ('HG002' if pset == 'HG002_Q100' else person)))].end_id.tolist()
            def cut(e):
                return e, L.window_from_public(dict(Ei.loc[e], end_id=e), cwd=wd)
            with ThreadPoolExecutor(6) as ex: wins = dict(ex.map(cut, ids))
            claimed = C['claimed_end_ids'].split(',')
            conc = [r for r in rows if r.get('verdict') == 'conclusive']
            pl = L.placement({r['read']: seqs[r['read']] for r in conc}, wins)
            for r in conc:
                r.update(L.placement_fields(pl.get(r['read']), claimed)); r['v1_verdict'] = 'conclusive'
                r['tier'] = L.read_tier(r)
            dfn = [r for r in conc if r['tier'] == 'definitive' and L.b(r['_pub']['counted'])]
            dv, nm, ds = L.definitive_claim(v['verdict'], dfn, molecule_key=lambda r: r['_pub']['molecule'])
            timings['placement_s'] = round(time.time() - t, 1)
            tiers_ok = all(r['tier'] == r['_pub']['tier'].replace('not_counted_', '') for r in conc)
            print(f"definitive tier: {len(dfn)} definitive reads ({nm} molecules, datasets {','.join(ds) or '-'}) -> {dv} "
                  f"(published {C['definitive_verdict']}); per-read tiers match: {tiers_ok}; placement on {len(wins)} public end windows")
            res.update(definitive_verdict=dv, n_definitive=len(dfn), tiers_match=tiers_ok)
            res['match'] = res['match'] and dv == C['definitive_verdict'] and tiers_ok
    print(f"\nRESULT: {'MATCH' if res['match'] else 'MISMATCH'}  ({ok_reads}/{len(cited)} reads identical; timings {timings})")
    json.dump(dict(res, reads=out), open(os.path.join(wd, 'replicate_summary.json'), 'w'), indent=1, default=str)
    sys.exit(0 if res['match'] else 1)


if __name__ == '__main__':
    main()
