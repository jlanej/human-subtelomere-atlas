"""H2 step 1 -- reciprocal exchange or one-way (duplicative / telomere-capture / BIR-like) transfer?

For every carrier of a foreign form with >= 10 carrier tips (form R>D: recipient arm R carries a terminal block that matches
CHM13 donor end D best), inspect the DONOR arm's telomeric end in the SAME haplotype assembly:
  donor_status  reciprocal                 donor tip itself is foreign with foreign arm == R (form D>R)
                donor_has_recipient_window some 5 kb window of the donor tip places foreign on chrR (internal/partial)
                donor_replaced_other       donor tip foreign towards another arm X != R (form D>X)
                donor_native               donor tip own (or own-tie) at 20 kb and no window places on chrR
                donor_end_absent           no chromosome-named telomeric tip of arm D in this haplotype assembly
                donor_unassessed           donor tip present but t20 'none'
  donor_native_kb        length of the leading run of own/tie 5 kb windows of the donor tip (kb from the subtelomere start)
  donor_native_covers    donor_native_kb*1000 >= donor-side extent of the transferred block (sw_for_dist; 'through' tips:
                         whole carrier window), i.e. the donor end still carries native sequence over the whole copied span
  blk_*                  the carrier's terminal block (subtel[0:switch]) aligned (mappy asm20) to (i) the same-haplotype donor
                         tip and (ii) the CHM13 donor end: covered fraction and identity (matches / aligned columns)
Recipient's lost terminal sequence elsewhere in the same haplotype:
  rlost_named_arms       other named tips of this haplotype with any 5 kb window placing foreign on chrR
  rlost_unnamed          unplaced/unnamed telomeric tips of this haplotype whose 20 kb subtelomere places on R (tipdup placement)
  rl_*                   direct search: the CHM13 recipient segment that was replaced (telomere end .. junction, sw_own_dist;
                         base- or window-resolved tips only) aligned to every end sequence of the haplotype (named tips,
                         unnamed telomeric tips, 20 kb ends of all contigs from r2scan/) except the carrier's own contig end;
                         best hit identity/covered fraction; the same segment against the other 47 CHM13 ends gives the
                         paralogue background (rl_chm13_para_*)
Breakpoint complementarity (independent of haplotype co-occurrence, which decays for unlinked ends): a reciprocal exchange
between R (break at jR from its telomere) and D (break at jD) also creates a D>R end whose own-side junction is jD and whose
foreign-side coordinate is jR.  For every junction cluster of the tested forms (median jR = sw_own_dist, jD = sw_for_dist of
its base-resolved members) all base- or window-resolved D>R tips in the release are searched for |own_dist - jD| <= 2 kb and
|for_dist - jR| <= 2 kb (results/reciprocal_breakpoint_complement.tsv).
Population baseline for reciprocity: among haplotypes NOT carrying R>D, the frequency of donor tips with form D>R
(Fisher exact test, carriers vs non-carriers).
HG002 Q100 v1.1 (complete diploid, 200 kb ends; subtel/foreign/hg002/place.*) is analysed the same way (dataset HG002_Q100).
Inputs : results/subtel_foreign_tips.tsv, subtel/foreign/r2/place.windows.tsv, subtel/foreign/r2/tips_input.tsv,
         tipdup/tip_placement.tsv, tipdup/tips_masked.fa, r2scan/<asm>.ends.fa, subtel/foreign/ref/chm13_ends500k.fa,
         subtel/foreign/hg002/{place.tips.tsv,place.windows.tsv,hg002_tips.tsv}
Outputs: results/reciprocal_per_carrier.tsv, results/reciprocal_per_form.tsv, results/reciprocal_summary.json
Usage  : python scripts/reciprocal.py
"""
import os, sys, json, collections, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jcommon as J
from scipy.stats import fisher_exact

MINCARRIERS = 10
NAMED_FORMS = ['20q>3q', '3q>11p', '6p>11p', '8p>1p', '16p>9q']

def chm13_telo_len(end, W=100):
    """length of the terminal (TTAGGG)n array of a CHM13 end in telomere frame (first 100 bp window < 50% hexamer cover)"""
    s = J.chm13_teloframe(end); ex, dg, run = J.its_tracks(s[:30000])
    for i in range(0, 30000 - W, W):
        if dg[i:i + W].sum() * 6 < 0.5 * W: return i
    return 30000

def aln_stats(a, seq, minlen=1000):
    """union coverage fraction of seq and identity (matches / aligned columns) of PRIMARY mappy hits >= minlen
    (secondary/paralogous placements excluded)"""
    import re
    n = len(seq); cov = np.zeros(n, bool); m = c = 0
    for h in a.map(seq, cs=True):
        if h.q_en - h.q_st < minlen or not h.is_primary: continue
        cov[h.q_st:h.q_en] = True
        for op, v in re.findall(r'([:*+\-~])([0-9]+|[a-z]+)', h.cs):
            if op == ':': m += int(v); c += int(v)
            elif op == '*': c += 1
            elif op in '+-': c += len(v)
    return (round(cov.mean(), 4) if n else np.nan), (round(m / c, 5) if c else np.nan)

def load_datasets():
    T = pd.read_csv(J.inp('results/subtel_foreign_tips.tsv'), sep='\t', low_memory=False)
    T['dataset'] = 'HPRC_R2'
    W = pd.read_csv(J.inp('subtel/foreign/r2/place.windows.tsv'), sep='\t')
    H = pd.read_csv(J.inp('subtel/foreign/hg002/place.tips.tsv'), sep='\t')
    H['asm'] = 'HG002_Q100.' + H.q.str.split('|').str[1].str.split('_').str[1]; H['sample'] = 'HG002'; H['dataset'] = 'HG002_Q100'
    H['hap'] = H.asm.str.split('.').str[1]
    H['foreign_arm'] = np.where(H.t20_state == 'foreign', H.t20_best_end.str[3:], None)
    H['form'] = np.where(H.t20_state == 'foreign', H.own + '>' + H.foreign_arm.astype(str), None)
    H['switch_subtel'] = H.switch_final - H.s0
    H['contig'] = H.q.str.split('|').str[1]; H['end'] = H.q.str.split('|').str[2]
    HW = pd.read_csv(J.inp('subtel/foreign/hg002/place.windows.tsv'), sep='\t')
    cols = ['q', 'asm', 'sample', 'dataset', 'contig', 'end', 'own', 's0', 'tip_len', 't20_state', 'foreign_arm', 'form', 'form_foreign_end',
            'switch_final', 'switch_res', 'switch_subtel', 'sw_own_dist', 'sw_for_dist', 'foreign_run_kb', 'profile', 'junction', 'junction_cluster']
    H['junction'] = np.where(H.switch_res == 'base', H.sw_own_dist, np.nan); H['junction_cluster'] = None
    A = pd.concat([T[cols], H[[c for c in cols]]], ignore_index=True)
    WW = pd.concat([W, HW], ignore_index=True)
    return A, WW

def windows_by_q(WW):
    D = collections.defaultdict(list)
    for r in WW[['q', 'win', 'best_end', 'state']].itertuples(index=False): D[r.q].append((r.win, r.best_end, r.state))
    for k in D: D[k].sort()
    return D

def classify_donor(dt, R, WB):
    """donor-end status of a haplotype given its named tips of arm D (list of rows) and the recipient arm R"""
    if not dt: return 'donor_end_absent', None
    dt = sorted(dt, key=lambda t: (not (t.t20_state == 'foreign' and t.foreign_arm == R), t.q)); t = dt[0]
    wl = WB.get(t.q, [])
    haswin = any(st == 'foreign' and be == 'chr' + R for w, be, st in wl)
    if t.t20_state == 'foreign' and t.foreign_arm == R: st = 'reciprocal'
    elif haswin: st = 'donor_has_recipient_window'
    elif t.t20_state == 'foreign': st = 'donor_replaced_other'
    elif t.t20_state in ('own', 'tie'): st = 'donor_native'
    else: st = 'donor_unassessed'
    return st, t

def native_kb(wl):
    k = 0
    for w, be, st in wl:
        if st in ('own', 'tie'): k = (w + 1) * 5
        else: break
    return k

# ---------- per-assembly worker: block-vs-donor alignment and lost-sequence search ----------
G = {}
def init(ctx): G.update(ctx)

def hap_end_seqs(asm, dataset, named, unnamed, hg002_seq):
    """all end sequences of one haplotype: name -> (seq, kind)"""
    out = {}
    if dataset == 'HG002_Q100':
        for q, s in hg002_seq.items():
            if q.split('|')[1].endswith(asm.split('.')[1]): out[q] = (s, 'named')
        return out
    for q, s in named.get(asm, {}).items(): out[q] = (s, 'named')
    for q, s in unnamed.get(asm, {}).items(): out[q] = (s, 'unnamed_tip')
    fn = J.inp('r2scan', asm + '.ends.fa')
    if os.path.exists(fn):
        for n, s in J.readfa(fn):
            ctg, end, L = n.split('|')
            key = f'{asm}|{ctg}|{end}'
            if key in out: continue
            out[key + '|ctgend'] = (s.upper() if end == 'L' else J.rc(s.upper()), 'contig_end')
    return out

def work_asm(args):
    import mappy
    asm, dataset, carriers, ends = args
    rows = []
    # one index over all haplotype end sequences for the lost-segment search
    tmp = os.path.join(J.WORK, 'tmp_idx'); os.makedirs(tmp, exist_ok=True)
    fa = os.path.join(tmp, f'{asm}.fa')
    with open(fa, 'w') as f:
        for k, (s, kind) in ends.items(): f.write(f'>{k}\n{s}\n')
    allidx = mappy.Aligner(fa, preset='asm20', n_threads=1, best_n=50)
    for c in carriers:
        d = dict(q=c['q'])
        seq = c['seq']; s0 = c['s0']; sub = seq[s0:]
        x = c['switch_subtel']; blk = sub[:int(x)] if (x == x and x is not None and x > 0) else sub
        d['blk_len'] = len(blk)
        dq = c.get('donor_q')
        if dq and dq in ends:
            ds = ends[dq][0]; ds = ds[c['donor_s0']:] if c.get('donor_s0') else ds
            a = mappy.Aligner(seq=ds, preset='asm20')
            d['blk_vs_hapdonor_cov'], d['blk_vs_hapdonor_ident'] = aln_stats(a, blk)
        d['blk_vs_chm13donor_cov'], d['blk_vs_chm13donor_ident'] = aln_stats(J.aligner('chr' + c['D']), blk)
        # lost recipient segment
        Jd = c['sw_own_dist']
        if c['switch_res'] in ('base', 'window') and Jd == Jd and Jd is not None:
            t0 = G['telo'][c['R']]; lost = J.chm13_teloframe('chr' + c['R'])[t0:int(Jd)]
            d['rl_len'] = len(lost)
            if len(lost) >= 2000:
                own_ctg = c['q'].split('|')[1] + '|' + c['q'].split('|')[2]
                best = None
                n = len(lost); per_t = collections.defaultdict(lambda: np.zeros(n, bool)); idt = collections.defaultdict(lambda: [0, 0])
                import re
                for h in allidx.map(lost, cs=True):
                    if h.q_en - h.q_st < 1000: continue
                    tn = h.ctg
                    if tn.split('|')[1] + '|' + tn.split('|')[2] == own_ctg: continue
                    per_t[tn][h.q_st:h.q_en] = True
                    for op, v in re.findall(r'([:*+\-~])([0-9]+|[a-z]+)', h.cs):
                        if op == ':': idt[tn][0] += int(v); idt[tn][1] += int(v)
                        elif op == '*': idt[tn][1] += 1
                        elif op in '+-': idt[tn][1] += len(v)
                for tn, cv in per_t.items():
                    cov = cv.mean(); ide = idt[tn][0] / idt[tn][1] if idt[tn][1] else 0
                    if best is None or cov * ide > best[1] * best[2]: best = (tn, cov, ide)
                if best:
                    tn = best[0]; kind = ends[tn][1]
                    d.update(rl_best_target=tn, rl_best_kind=kind, rl_best_arm=G['arm_of'].get(tn.replace('|ctgend', ''), ''),
                             rl_best_cov=round(best[1], 4), rl_best_ident=round(best[2], 5))
                cp = G['para'].get((c['R'], int(Jd)))
                if cp is None:
                    cp = (None, 0, 0)
                    for e in J.chm13_ends():
                        if e == 'chr' + c['R']: continue
                        cov, ide = aln_stats(J.aligner(e), lost)
                        if cov == cov and ide == ide and cov * ide > cp[1] * cp[2]: cp = (e, cov, ide)
                    G['para'][(c['R'], int(Jd))] = cp
                d.update(rl_chm13_para_end=cp[0], rl_chm13_para_cov=round(cp[1], 4), rl_chm13_para_ident=round(cp[2], 5) if cp[2] else 0)
        rows.append(d)
    os.remove(fa)
    return rows

def main():
    A, WW = load_datasets(); WB = windows_by_q(WW)
    F = A[A.t20_state == 'foreign']
    nform = F[F.dataset == 'HPRC_R2'].form.value_counts()
    forms = sorted(set(nform[nform >= MINCARRIERS].index) | set(NAMED_FORMS))
    TP = pd.read_csv(J.inp('tipdup/tip_placement.tsv'), sep='\t')
    unn = TP[TP.arm.isna()].copy(); unn['place_arm'] = unn.sub_chr.str[3:] + unn.sub_end
    unn_by_asm = {a: g for a, g in unn.groupby('asm')}
    by_asm_arm = collections.defaultdict(list)
    for r in A.itertuples(): by_asm_arm[(r.asm, r.own)].append(r)
    by_asm = collections.defaultdict(list)
    for r in A.itertuples(): by_asm[r.asm].append(r)
    C = F[F.form.isin(forms)].copy()
    recs = []
    for r in C.itertuples():
        R, D = r.own, r.foreign_arm
        d = dict(q=r.q, dataset=r.dataset, asm=r.asm, sample=r.sample, form=r.form, recipient=R, donor=D, form_foreign_end=r.form_foreign_end,
                 switch_res=r.switch_res, junction_own_dist=r.sw_own_dist, block_donor_dist=r.sw_for_dist, carrier_tip_len=r.tip_len)
        dt = by_asm_arm.get((r.asm, D), [])
        d['donor_tips'] = len(dt)
        # if several named tips of arm D exist (duplicate contigs), a reciprocal one is preferred, else the first by name
        st, t = classify_donor(dt, R, WB)
        if t is None: d['donor_status'] = st
        else:
            wl = WB.get(t.q, []); nk = native_kb(wl)
            ext = r.sw_for_dist if r.switch_res in ('base', 'window') and r.sw_for_dist == r.sw_for_dist else r.tip_len
            d.update(donor_status=st, donor_q=t.q, donor_t20_state=t.t20_state, donor_form=t.form, donor_tip_len=t.tip_len, donor_profile=t.profile,
                     donor_native_kb=nk, donor_native_covers=bool(nk * 1000 >= ext) if ext == ext else None,
                     donor_recipient_window_kb=','.join(str(w * 5) for w, be, s in wl if s == 'foreign' and be == 'chr' + R))
        # R-like sequence elsewhere in this haplotype (placement-based)
        arms = sorted({o.own for o in by_asm[r.asm] if o.q != r.q and any(st == 'foreign' and be == 'chr' + R for w, be, st in WB.get(o.q, []))})
        d['rlost_named_arms'] = ','.join(arms); d['rlost_in_donor'] = D in arms
        g = unn_by_asm.get(r.asm)
        d['rlost_unnamed'] = int((g.place_arm == R).sum()) if g is not None else 0
        recs.append(d)
    X = pd.DataFrame(recs)
    # ---- sequence-level part (parallel per assembly) ----
    need = set(X.q) | set(X.donor_q.dropna())
    seqs = J.load_tip_seqs(need)
    named_by_asm = collections.defaultdict(dict)
    with open(J.inp('subtel/foreign/r2/tips_input.tsv')) as f:
        f.readline(); asms = set(X[X.dataset == 'HPRC_R2'].asm)
        for line in f:
            q, own, s0, src, s = line.rstrip('\n').split('\t')
            a = q.split('|')[0]
            if a in asms: named_by_asm[a][q] = s
    unnamed_by_asm = collections.defaultdict(dict); un = set(unn.q)
    for n, s in J.readfa(J.inp('tipdup/tips_masked.fa')):
        if n in un and n.split('|')[0] in asms: unnamed_by_asm[n.split('|')[0]][n] = s.upper()
    HT = pd.read_csv(J.inp('subtel/foreign/hg002/hg002_tips.tsv'), sep='\t')
    hg002_seq = dict(zip(HT.q, HT.seq)); hg002_s0 = dict(zip(HT.q, HT.s0))
    for q in need:
        if q in hg002_seq: seqs[q] = dict(seq=hg002_seq[q], s0=int(hg002_s0[q]), src='200k')
    arm_of = dict(zip(A.q, A.own))
    telo = {e[3:]: chm13_telo_len(e) for e in J.chm13_ends()}
    s0_of = {q: v['s0'] for q, v in seqs.items()}
    jobs = []; swsub = dict(zip(A.q, A.switch_subtel))
    for (asm, ds), g in X.groupby(['asm', 'dataset']):
        cl = []
        for r in g.itertuples():
            cl.append(dict(q=r.q, seq=seqs[r.q]['seq'], s0=seqs[r.q]['s0'], switch_subtel=swsub[r.q], R=r.recipient, D=r.donor,
                           sw_own_dist=r.junction_own_dist, switch_res=r.switch_res,
                           donor_q=r.donor_q if isinstance(r.donor_q, str) else None, donor_s0=s0_of.get(r.donor_q) if isinstance(r.donor_q, str) else None))
        jobs.append((asm, ds, cl, hap_end_seqs(asm, ds, named_by_asm, unnamed_by_asm, hg002_seq)))
    del named_by_asm, unnamed_by_asm
    if os.environ.get('JUNC_DEBUG_MAX_ASM'):  # debugging only; never set by run.sh
        jobs = sorted(jobs, key=lambda j: j[0])[:int(os.environ['JUNC_DEBUG_MAX_ASM'])]
    ctx = dict(arm_of=arm_of, telo=telo, para={})
    with Pool(J.THREADS, initializer=init, initargs=(ctx,)) as pool:
        out = [row for rows in pool.imap_unordered(work_asm, sorted(jobs, key=lambda j: j[0]), chunksize=1) for row in rows]
    S = pd.DataFrame(out)
    X = X.merge(S, on='q', how='left').sort_values(['form', 'dataset', 'asm', 'q'])
    X['rl_found'] = (X.rl_best_cov >= 0.5) & (X.rl_best_ident >= 0.98)
    X['rl_found_above_paralogue'] = X.rl_found & (X.rl_best_cov * X.rl_best_ident > X.rl_chm13_para_cov.fillna(0) * X.rl_chm13_para_ident.fillna(0) + 0.05)
    X.to_csv(os.path.join(J.RES, 'reciprocal_per_carrier.tsv'), sep='\t', index=False)
    # ---- per form ----
    R2 = A[A.dataset == 'HPRC_R2']
    rows = []
    for f, g in X[X.dataset == 'HPRC_R2'].groupby('form'):
        R, D = f.split('>')
        vc = g.donor_status.value_counts()
        assessed = g[~g.donor_status.isin(['donor_end_absent', 'donor_unassessed'])]
        # baseline: haplotypes whose R tip is assessed (own/tie/foreign) but does NOT carry R>D; their D-end status
        carrier_asm = set(g.asm)
        nc = R2[(R2.own == R) & R2.t20_state.isin(['own', 'tie', 'foreign']) & (R2.form != f) & ~R2.asm.isin(carrier_asm)].asm.unique()
        bst = collections.Counter(classify_donor(by_asm_arm.get((a, D), []), R, WB)[0] for a in nc)
        base_ass = sum(v for k, v in bst.items() if k not in ('donor_end_absent', 'donor_unassessed'))
        n_base_recip = bst.get('reciprocal', 0); n_car_recip = int((assessed.donor_status == 'reciprocal').sum())
        tab = [[n_car_recip, len(assessed) - n_car_recip], [n_base_recip, base_ass - n_base_recip]]
        orr, p = fisher_exact(tab) if len(assessed) and base_ass else (np.nan, np.nan)
        n_car_nat = int((assessed.donor_status == 'donor_native').sum()); n_base_nat = bst.get('donor_native', 0)
        orn, pn = fisher_exact([[n_car_nat, len(assessed) - n_car_nat], [n_base_nat, base_ass - n_base_nat]]) if len(assessed) and base_ass else (np.nan, np.nan)
        n_car_abs = int((g.donor_status == 'donor_end_absent').sum()); n_base_abs = bst.get('donor_end_absent', 0)
        # placement baseline for 'recipient-like sequence elsewhere': non-carrier haplotypes with another named tip (not the R tip)
        # carrying a 5 kb window that places foreign on chrR
        nb_else = sum(any(o.own != R and any(st == 'foreign' and be == 'chr' + R for w, be, st in WB.get(o.q, [])) for o in by_asm[a]) for a in nc)
        ora, pa = fisher_exact([[n_car_abs, len(g) - n_car_abs], [n_base_abs, len(nc) - n_base_abs]]) if len(nc) else (np.nan, np.nan)
        gl = g[g.rl_len >= 2000]
        rows.append(dict(form=f, recipient=R, donor=D, carriers=len(g), samples=g['sample'].nunique(),
                         donor_end_absent=int(vc.get('donor_end_absent', 0)), donor_unassessed=int(vc.get('donor_unassessed', 0)), donor_assessed=len(assessed),
                         donor_native=int(vc.get('donor_native', 0)), reciprocal=n_car_recip, donor_has_recipient_window=int(vc.get('donor_has_recipient_window', 0)),
                         donor_replaced_other=int(vc.get('donor_replaced_other', 0)),
                         frac_native_of_assessed=round(vc.get('donor_native', 0) / len(assessed), 4) if len(assessed) else np.nan,
                         donor_native_covers_block=int((assessed.donor_native_covers == True).sum()),
                         noncarrier_haps=len(nc), noncarrier_donor_assessed=base_ass, noncarrier_donor_absent=n_base_abs,
                         noncarrier_donor_native=n_base_nat, noncarrier_reciprocal_form=n_base_recip,
                         noncarrier_donor_has_recipient_window=bst.get('donor_has_recipient_window', 0), noncarrier_donor_replaced_other=bst.get('donor_replaced_other', 0),
                         reciprocal_frac_carriers=round(n_car_recip / len(assessed), 4) if len(assessed) else np.nan,
                         reciprocal_frac_noncarriers=round(n_base_recip / base_ass, 4) if base_ass else np.nan,
                         fisher_reciprocal_or=orr, fisher_reciprocal_p=p,
                         native_frac_noncarriers=round(n_base_nat / base_ass, 4) if base_ass else np.nan, fisher_native_or=orn, fisher_native_p=pn,
                         absent_frac_carriers=round(n_car_abs / len(g), 4), absent_frac_noncarriers=round(n_base_abs / len(nc), 4) if len(nc) else np.nan,
                         fisher_absent_or=ora, fisher_absent_p=pa,
                         blk_hapdonor_ident_median=g.blk_vs_hapdonor_ident.median(), blk_hapdonor_cov_median=g.blk_vs_hapdonor_cov.median(),
                         blk_chm13donor_ident_median=g.blk_vs_chm13donor_ident.median(),
                         rlost_tested=len(gl), rlost_named_elsewhere=int((gl.rlost_named_arms.fillna('') != '').sum()),
                         carriers_R_window_elsewhere=int((g.rlost_named_arms.fillna('') != '').sum()),
                         noncarrier_R_window_elsewhere=int(nb_else), noncarrier_R_window_elsewhere_frac=round(nb_else / len(nc), 4) if len(nc) else np.nan,
                         rlost_in_donor=int(gl.rlost_in_donor.sum()), rlost_unnamed=int((gl.rlost_unnamed > 0).sum()),
                         rl_found=int(gl.rl_found.sum()), rl_found_above_paralogue=int(gl.rl_found_above_paralogue.sum()),
                         rl_len_median=gl.rl_len.median()))
    PF = pd.DataFrame(rows).sort_values('carriers', ascending=False)
    # ---- breakpoint complementarity ----
    TOL = 2000; bc = []
    B = R2[R2.t20_state == 'foreign']
    for f in forms:
        R, D = f.split('>')
        g = B[(B.form == f) & (B.switch_res == 'base') & B.junction_cluster.notna()]
        rev = B[(B.form == f'{D}>{R}') & B.switch_res.isin(['base', 'window'])]
        for cid, h in g.groupby('junction_cluster'):
            jR, jD = h.sw_own_dist.median(), h.sw_for_dist.median()
            m = rev[((rev.sw_own_dist - jD).abs() <= TOL) & ((rev.sw_for_dist - jR).abs() <= TOL)]
            bc.append(dict(form=f, junction_cluster=cid, cluster_tips=len(h), cluster_samples=h['sample'].nunique(), jR_median=jR, jD_median=jD,
                           reverse_form=f'{D}>{R}', reverse_form_tips=len(rev), reverse_tips_complementary=len(m),
                           reverse_samples_complementary=m['sample'].nunique(), reverse_complementary_q=','.join(sorted(m.q))[:2000]))
    BC = pd.DataFrame(bc); BC.to_csv(os.path.join(J.RES, 'reciprocal_breakpoint_complement.tsv'), sep='\t', index=False)
    PF.to_csv(os.path.join(J.RES, 'reciprocal_per_form.tsv'), sep='\t', index=False)
    XR = X[X.dataset == 'HPRC_R2']; ass = XR[~XR.donor_status.isin(['donor_end_absent', 'donor_unassessed'])]
    gl = XR[XR.rl_len >= 2000]
    summ = dict(forms=len(PF), carriers=len(XR), carriers_samples=int(XR['sample'].nunique()),
                donor_status=XR.donor_status.value_counts().to_dict(), donor_assessed=len(ass),
                donor_native_frac=round((ass.donor_status == 'donor_native').mean(), 4),
                reciprocal=int((ass.donor_status == 'reciprocal').sum()), reciprocal_frac=round((ass.donor_status == 'reciprocal').mean(), 4),
                forms_with_any_reciprocal=int((PF.reciprocal > 0).sum()),
                forms_reciprocal_enriched_p05=int(((PF.fisher_reciprocal_p < 0.05) & (PF.reciprocal_frac_carriers > PF.reciprocal_frac_noncarriers)).sum()),
                forms_reciprocal_enriched_bonf=int(((PF.fisher_reciprocal_p < 0.05 / len(PF)) & (PF.reciprocal_frac_carriers > PF.reciprocal_frac_noncarriers)).sum()),
                forms_native_depleted_bonf=int(((PF.fisher_native_p < 0.05 / len(PF)) & (PF.frac_native_of_assessed < PF.native_frac_noncarriers)).sum()),
                forms_absent_enriched_bonf=int(((PF.fisher_absent_p < 0.05 / len(PF)) & (PF.absent_frac_carriers > PF.absent_frac_noncarriers)).sum()),
                pooled=dict(carriers=dict(assessed=int(PF.donor_assessed.sum()), native=int(PF.donor_native.sum()), reciprocal=int(PF.reciprocal.sum()),
                                          has_recipient_window=int(PF.donor_has_recipient_window.sum()), replaced_other=int(PF.donor_replaced_other.sum()),
                                          absent=int(PF.donor_end_absent.sum()), total=int(PF.carriers.sum())),
                            noncarriers=dict(assessed=int(PF.noncarrier_donor_assessed.sum()), native=int(PF.noncarrier_donor_native.sum()),
                                             reciprocal=int(PF.noncarrier_reciprocal_form.sum()), has_recipient_window=int(PF.noncarrier_donor_has_recipient_window.sum()),
                                             replaced_other=int(PF.noncarrier_donor_replaced_other.sum()), absent=int(PF.noncarrier_donor_absent.sum()), total=int(PF.noncarrier_haps.sum()))),
                donor_native_covers_block=int((ass.donor_native_covers == True).sum()),
                blk_hapdonor_ident_median=float(XR.blk_vs_hapdonor_ident.median()), blk_hapdonor_cov_median=float(XR.blk_vs_hapdonor_cov.median()),
                complement_clusters_tested=len(BC), complement_clusters_with_reverse_form=int((BC.reverse_form_tips > 0).sum()),
                complement_clusters_matched=int((BC.reverse_tips_complementary > 0).sum()),
                complement_clusters_matched_list=BC[BC.reverse_tips_complementary > 0][['junction_cluster', 'cluster_samples', 'reverse_form', 'reverse_tips_complementary', 'reverse_samples_complementary']].to_dict('records'),
                complement_tips_in_matched=int(BC.reverse_tips_complementary.sum()),
                R_window_elsewhere=dict(carriers=int(PF.carriers_R_window_elsewhere.sum()), carriers_total=int(PF.carriers.sum()),
                                        noncarriers=int(PF.noncarrier_R_window_elsewhere.sum()), noncarriers_total=int(PF.noncarrier_haps.sum())),
                rlost_tested=len(gl), rlost_named_elsewhere=int((gl.rlost_named_arms.fillna('') != '').sum()),
                rlost_in_donor=int(gl.rlost_in_donor.sum()), rlost_unnamed=int((gl.rlost_unnamed > 0).sum()),
                rl_found=int(gl.rl_found.sum()), rl_found_above_paralogue=int(gl.rl_found_above_paralogue.sum()),
                hg002=X[X.dataset == 'HG002_Q100'][['q', 'form', 'donor_status', 'donor_form', 'rlost_named_arms', 'rl_best_arm', 'rl_best_cov', 'rl_best_ident']].to_dict('records'),
                params=dict(min_carriers=MINCARRIERS, named_forms=NAMED_FORMS, rl_found='best non-self hit covers >= 50% of the lost CHM13 segment at >= 98% identity',
                            rl_found_above_paralogue='rl_found and cov*ident exceeds the best other-CHM13-end cov*ident by > 0.05'))
    from scipy.stats import chi2_contingency
    pc, pn_ = summ['pooled']['carriers'], summ['pooled']['noncarriers']
    cats = ['native', 'reciprocal', 'has_recipient_window', 'replaced_other']
    chi = chi2_contingency([[pc[c] for c in cats], [pn_[c] for c in cats]])
    summ['pooled']['chi2_donor_status'] = dict(categories=cats, chi2=float(chi[0]), p=float(chi[1]), dof=int(chi[2]))
    summ['pooled']['frac'] = {k: {c: round(v[c] / v['assessed'], 4) for c in cats} for k, v in (('carriers', pc), ('noncarriers', pn_))}
    json.dump(summ, open(os.path.join(J.RES, 'reciprocal_summary.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps({k: v for k, v in summ.items() if k not in ('hg002', 'params', 'pooled', 'complement_clusters_matched_list')}, default=str))
    print(BC[BC.reverse_tips_complementary > 0].to_string())

if __name__ == '__main__':
    main()
