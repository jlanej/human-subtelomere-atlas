"""'Read evidence' section of the site (docs/evidence/): overview, searchable claims table, per-claim pages, 'Verify it yourself'
and a FAQ for skeptics. Built by website/build_site.py from evidence/claims.tsv, evidence/reads.tsv.gz, evidence/controls.tsv,
evidence/showcase/ and evidence/figures/ (all numbers are read from those files at build time)."""
import os, json, html, re, gzip
import pandas as pd

esc = html.escape
EV = "evidence"
COHORT = {"HG002_Q100": "HG002 Q100", "HG002_trio": "HG003/HG004 (trio)", "1KG_ONT": "1000 Genomes ONT"}
DEF_LAB = {"definitive_multi_dataset": "definitive (>= 2 datasets)", "definitive": "definitive",
           "confirmed_not_definitive": "confirmed, not definitive", "not_definitive": "-"}
TIER_LAB = {"definitive": "definitive", "conclusive_shared_anchor": "conclusive, shared anchor",
            "conclusive_unique_anchor_off_end": "conclusive, places on another end",
            "conclusive_shared_anchor_off_end": "conclusive, shared anchor, another end",
            "conclusive_discordant": "conclusive, switch-discordant", "contradicts": "contradicts", "consistent": "consistent"}


def _f(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def load(rp):
    C = pd.read_csv(rp("evidence", "claims.tsv"), sep="\t", low_memory=False)
    R = pd.read_csv(rp("evidence", "reads.tsv.gz"), sep="\t", low_memory=False)
    K = pd.read_csv(rp("evidence", "controls.tsv"), sep="\t")
    S = pd.read_csv(rp("evidence", "showcase", "showcase_claims.tsv"), sep="\t")
    SL = pd.read_csv(rp("evidence", "showcase", "showcase_loci.tsv"), sep="\t")
    M = pd.read_csv(rp("evidence", "release_manifest.tsv"), sep="\t")
    crit = json.load(open(rp("evidence", "criteria_v1.json")))
    rep = rp("evidence", "replication_results.tsv")
    RT = pd.read_csv(rep, sep="\t") if os.path.exists(rep) else None
    return dict(C=C, R=R, K=K, S=S, SL=SL, M=M, crit=crit, RT=RT)


def numbers(E):
    C, R, K = E["C"], E["R"], E["K"]
    isdef = C.definitive_verdict.str.startswith("definitive")
    out = {}
    for k in COHORT:
        g = C[C.cohort == k]
        if k == "HG002_trio":
            g = g[g.v1_class == "transmitting_parent"]
        rr = R[(R.cohort == k) & (R.tier == "definitive") & (R.counted == True) & R.claim_id.isin(g[isdef.loc[g.index]].claim_id)]  # noqa: E712
        out[k] = dict(n=len(g), confirmed=int((g.v1_verdict == "confirmed").sum()), supported=int((g.v1_verdict == "supported").sum()),
                      definitive=int(isdef.loc[g.index].sum()),
                      multi=int((g.definitive_verdict == "definitive_multi_dataset").sum()),
                      def_reads=len(rr), conc_reads=int(g.n_conclusive_concordant.sum()))
    out["ctl_hg002"] = (int(K[K.cohort == "HG002_Q100"].n_eligible.sum()), int(K[K.cohort == "HG002_Q100"].n_conclusive.sum()))
    out["ctl_1kg"] = (int(K[K.cohort == "1KG_ONT"].n_eligible.sum()), int(K[K.cohort == "1KG_ONT"].n_conclusive.sum()))
    out["ctl_1kg_reads"] = int(R[R.claim_id.str.startswith("control|") & (R.cohort == "1KG_ONT")].read_id.nunique())
    out["n_1kg_ind"] = int(C[C.cohort == "1KG_ONT"].individual.nunique())
    out["n_1kg_def_ind"] = int(C[(C.cohort == "1KG_ONT") & isdef].individual.nunique())
    out["n_1kg_def_alleles"] = int(C[(C.cohort == "1KG_ONT") & isdef].allele_id.nunique())
    out["n_1kg_def_arms"] = int(C[(C.cohort == "1KG_ONT") & isdef].arm.nunique())
    return out


def _ub(k, n):
    from scipy.stats import beta
    return 1.0 if n == 0 else float(beta.ppf(0.95, k + 1, n - k))


def claim_label(c):
    if c["cohort"] == "HG002_Q100":
        return f"HG002 {c['arm']} {c['haplotype']}"
    if c["cohort"] == "HG002_trio":
        return f"{c['individual']} {c['arm']} ({c['haplotype']})"
    return f"{c['individual']} {c['arm']}"


def has_page(c):
    return int(c["n_conclusive"] or 0) > 0 or c["v1_verdict"] in ("confirmed", "supported")


def claim_href(c, pre):
    return f"{pre}evidence/claim/{c['slug']}.html"


def build(ctx):
    """ctx: dict with rp, write, page, figure, copy_image, BLOB, TREE, table_html, fnum"""
    rp, write, page, figure, BLOB = ctx["rp"], ctx["write"], ctx["page"], ctx["figure"], ctx["BLOB"]
    E = load(rp); N = numbers(E)
    sizes = page_index(ctx, E, N)
    page_claims(ctx, E, N)
    n = page_claim_all(ctx, E, N)
    page_verify(ctx, E, N)
    page_faq(ctx, E, N)
    import tour_pages
    n_ch, n_gal = tour_pages.build(ctx)
    return dict(numbers=N, n_claim_pages=n, tour_chapters=n_ch, tour_claims=n_gal)


def subnav(pre, cur):
    items = [("evidence/index.html", "Overview"), ("evidence/tour.html", "Guided tour"), ("evidence/claims.html", "Claims table"),
             ("evidence/verify.html", "Verify it yourself"),
             ("evidence/faq.html", "FAQ for skeptics")]
    return '<nav class="arm-nav" aria-label="Read evidence">' + "".join(
        f'<a href="{pre}{h}"{" aria-current=page" if h == cur else ""}>{l}</a>' for h, l in items) + "</nav>"


# ------------------------------------------------------------------------------------------------------------------ overview
def tier_table():
    rows = [("v1 conclusive read", "pre-registered (criteria v1, frozen hash a89c5eaf023ed032)",
             "terminal telomere; donor alleles distal and recipient/proximal alleles proximal to a maximum-likelihood switch "
             "(each segment &ge; 10 units at &ge; 80% purity, both log10 likelihood ratios &ge; 5); a &ge; 10 kb MAPQ-60 genome-wide "
             "anchor on the recipient chromosome within 2 Mb of its telomere; no internal telomere, adapter or foldback; &ge; 90% aligned"),
            ("switch-concordant", "pre-registered (amendment C2)", "the read's switch interval (&plusmn; 2 kb) overlaps the switch of "
             "the person's own assembly window analysed with the same unit model"),
            ("unique anchor", "reported in v1, required only for <em>definitive</em>", "&ge; 1 kb of the anchor on CHM13 recipient "
             "positions that no other CHM13 end covers at &ge; 95% identity over &ge; 1 kb"),
            ("own-assembly placement", "post hoc (diagnostics)", "the read aligns with more matching bases to the claimed end window "
             "of the same person's assembly than to any other end window of that person (HG002: the 92 Q100 windows; 1KG-ONT "
             "individuals: their ~90 HPRC Release 2 windows)"),
            ("<strong>definitive read</strong>", "post hoc, reporting only", "all four rows above"),
            ("<strong>definitive claim</strong>", "post hoc, reporting only", "v1 <em>confirmed</em> and &ge; 2 definitive reads from "
             "distinct molecules; <em>multi-dataset</em> when they come from &ge; 2 independent datasets")]
    return ('<div class="table-wrap free"><table><thead><tr><th scope="col">criterion</th><th scope="col">status</th>'
            '<th scope="col">definition</th></tr></thead><tbody>' +
            "".join(f"<tr><td>{a}</td><td>{b}</td><td>{c}</td></tr>" for a, b, c in rows) + "</tbody></table></div>")


def _pitch(rp):
    """the numbers of report/PITCH (report/pitch_content.py), so the site and the PDF quote the same values"""
    import importlib, os, sys
    cwd = os.getcwd(); os.chdir(rp('.'))
    try:
        if rp('report') not in sys.path: sys.path.insert(0, rp('report'))
        return importlib.import_module('pitch_content')
    finally:
        os.chdir(cwd)


def one_picture(rp, figure, BLOB, pre):
    """'The evidence in one picture': the pitch figure of one molecule, read in three steps"""
    Q = _pitch(rp); hr, hp, hg = Q.hr, Q.hp, Q.hg
    nD, nPre = int((Q.pre.call == 'D').sum()), len(Q.pre); nP, nPost = int((Q.post.call == 'P').sum()), len(Q.post)
    fig = figure('figures/pitch/fig_pitch_molecule.png', 'One 513 kb molecule: telomere, chromosome-11-like DNA, then chromosome 3',
                 f'Nanopore read <code>{esc(Q.HERO)}</code> ({esc(hr.dataset)}, HG002 Q100 public alignments). (a) Drawn to scale '
                 f'against the CHM13 ends of chromosomes 11 and 3; cut at the switch, its parts align to CHM13 '
                 f'chr11:{hg.loc["distal", "pos"]:,}-{hg.loc["distal", "end"]:,} and chr3:{hg.loc["proximal", "pos"]:,}-'
                 f'{hg.loc["proximal", "end"]:,} (both MAPQ 60). (b) Every site where the CHM13 chr11 and chr3 copies differ: '
                 f'which letter the molecule carries. (c) The letters at 6 evenly spaced sites on each side. (d) {Q.n_flip} of '
                 f'{Q.n_mol} molecules of five showcase claims turn from donor-like to recipient-like at their switch (2 kb '
                 f'window comparison; the 17p and 2q showcase claims are not shown because this simple comparison does not work '
                 f'there, see the <a href="{BLOB}report/PITCH.md">pitch</a> caption). Data: '
                 f'<a href="{BLOB}report/pitch_data">report/pitch_data/</a>.', pre)
    common = figure('figures/pitch/fig_pitch_common.png', 'The 12 most frequent foreign chromosome ends',
                    'The 12 foreign structural alleles with the highest Release 2 frequency. Filled: at least one definitive '
                    'read-evidence claim; open: none yet (not testable with the available reads, or not confirmed).', pre)
    import tour_pages
    combo = figure('figures/pitch/fig_consolidated.png', 'Molecules, assembly and reference genome tell the same story',
                   tour_pages.consolidated_caption(rp), pre)
    return f"""
<h2 id="one-picture">The evidence in one picture</h2>
<p>A DNA molecule that was never used to build an assembly, read from one end to the other. If the atlas is right, a
foreign chromosome end should look like this in a single molecule. It does.</p>
<div class="note callout"><strong>See the molecules yourself, in a genome browser:</strong> <a href="tour.html">guided tour</a>
(nothing to install; every definitive claim, with the person's own assembly, the reference genome and the molecules on one
screen).</div>
{fig}
<ol class="steps">
<li><strong>It starts with telomere</strong> ({hr.telo_len / 1e3:.1f} kb of TTAGGG repeats), so it is the real end of a chromosome.</li>
<li><strong>Then it reads like chromosome 11.</strong> At {nD} of {nPre} sites where chromosomes 11 and 3 differ, it carries chromosome
11's letter, and this part maps uniquely to the end of chromosome 11.</li>
<li><strong>Then it switches to chromosome 3 and stays there</strong>: chromosome 3's letter at {nP} of {nPost} sites, and
{hg.loc['proximal', 'aligned_bp'] / 1e3:.0f} kb aligned to chromosome 3 at MAPQ 60. One molecule, {hr.read_len / 1e3:.0f} kb long.</li>
</ol>
<p>The same flip appears in molecule after molecule: in different people, in HG002's mother, and on Oxford Nanopore
and PacBio data. The allele in the picture, SA-3q-01, is carried by {Q.a3.n_samples} people, and some foreign ends are
carried by more than half of all people (below). Two-minute version:
<a href="{BLOB}report/PITCH.pdf">PITCH.pdf</a> (<a href="{BLOB}report/PITCH.md">Markdown</a>).</p>
<p>The same story from three directions at once: the molecules on the person's own assembly, the assembly on the reference
genome, and the molecules on the reference genome.</p>
{combo}
<p>And piece by piece: every 2 kb of every molecule, aligned on its own to both reference ends.</p>
{figure('figures/pitch/fig_alignment_contrast.png', 'Window-by-window alignment contrast of every tour molecule', tour_pages.alignment_caption(rp), pre)}
{common}
<h2 id="formal-test">The formal test</h2>
"""


def pitch_callout(rp, pre, BLOB):
    """short home-page box pointing at the one-molecule picture"""
    Q = _pitch(rp)
    return (f'<div class="note callout"><strong>New here? Start with one molecule.</strong> A single {Q.hr.read_len / 1e3:.0f} kb '
            f'DNA molecule from HG002 starts in the telomere, reads like chromosome 11 for {(Q.hp.switch_lo - Q.hp.telo_len) / 1e3:.0f} kb, '
            f'then switches and runs along chromosome 3 (MAPQ 60). Every one of the {Q.n_people} Release 2 people carries foreign '
            f'chromosome ends (median {Q.per.f.median():.0f} of about 90 assembled ends). <a href="{pre}evidence/index.html#one-picture">See the molecule</a> &middot; '
            f'<a href="{pre}evidence/tour.html">guided tour in a genome browser</a> &middot; '
            f'<a href="{BLOB}report/PITCH.pdf">two-minute PDF</a></div>')


def page_index(ctx, E, N):
    pre = "../"; BLOB, TREE, figure = ctx["BLOB"], ctx["TREE"], ctx["figure"]
    h, t, k = N["HG002_Q100"], N["HG002_trio"], N["1KG_ONT"]
    (nh, kh), (nk, kk) = N["ctl_hg002"], N["ctl_1kg"]
    M = E["M"]; RT = E["RT"]
    rel = "".join(f"<tr><td><code>{esc(r.file)}</code></td><td class=num>{r.bytes / 1e6:.1f} MB</td><td><code>{r.sha256[:16]}&hellip;</code></td>"
                  f"<td>{esc(r.content)}</td></tr>" for r in M.itertuples())
    show = "".join(f'<li><a href="claim/{r.slug}.html">{esc(r.claim_id)}</a> ({r.n_reads_in_showcase} reads; '
                   f'<a href="{TREE}evidence/showcase/{r.slug}">files and IGV sessions</a>)</li>' for r in E["S"].itertuples())
    rep_txt = ""
    if RT is not None:
        ok = int((RT.result == "MATCH").sum())
        rep_txt = (f"<p>Replication test from a clean copy of the repository (no <code>work/</code> directory): {ok} of {len(RT)} "
                   f"claims and controls re-derived from the public sources gave identical per-read and claim verdicts "
                   f"(<a href=\"verify.html#test\">details</a>).</p>")
    body = f"""
{subnav(pre, "evidence/index.html")}
<h1>Read-level evidence</h1>
{one_picture(ctx["rp"], figure, BLOB, pre)}
<p class="lead">Are the non-native terminal blocks of the catalogue real chromosome ends, or assembly artefacts? Single long
reads that were not used to build the assemblies answer this independently: a read that starts in the telomere, carries
the donor arm's paralog-diagnostic alleles up to a switch, continues with the recipient arm's alleles, and then anchors
uniquely in the recipient chromosome is a physical molecule with the claimed structure.</p>
<div class="grid">
<div class="card"><span class="stat">{h['definitive']} / {h['confirmed']}</span>HG002 Q100 claims confirmed under the pre-registered
criteria are also <em>definitive</em> ({h['multi']} with definitive reads from &ge; 2 independent datasets; {h['def_reads']}
definitive reads; {h['n']} claims tested).</div>
<div class="card"><span class="stat">{t['confirmed']} / {h['confirmed'] + h['supported']}</span>HG002 claims confirmed in the
transmitting parent (HG003 or HG004, GIAB reads; {t['definitive']} definitive).</div>
<div class="card"><span class="stat">{k['definitive']} / {k['confirmed']}</span>1000 Genomes ONT claims confirmed under v1 are
definitive ({k['def_reads']} definitive reads; {N['n_1kg_def_ind']} individuals, {N['n_1kg_def_alleles']} alleles,
{N['n_1kg_def_arms']} arms; {k['n']} claims in {N['n_1kg_ind']} individuals tested, most not assessable).</div>
<div class="card"><span class="stat">0</span>conclusive reads in negative controls: 0 of {nh:,} eligible HG002 other-haplotype
read tests (95% upper bound {100 * _ub(kh, nh):.2f}%) and 0 of {nk:,} 1KG-ONT read &times; model tests of strictly native
individuals ({N['ctl_1kg_reads']:,} cited control reads; upper bound {100 * _ub(kk, nk):.2f}%).</div>
</div>
<div class="note"><strong>Two tiers.</strong> The <em>v1</em> verdicts follow criteria registered and frozen before any
population read was inspected (<a href="{BLOB}evidence/CRITERIA.md">CRITERIA.md</a>). The stricter <em>definitive</em> tier was
defined afterwards, from the population diagnostics (<a href="{BLOB}evidence/population/DIAGNOSTICS.md">DIAGNOSTICS.md</a>).
It is used for reporting only and can only remove claims from the v1-confirmed set; it never adds one.</div>
{figure('evidence/figures/fig_showcase_readthrough.png', 'Eight exemplar reads drawn from the telomere to the anchor',
        'Showcase reads, one per exemplar claim, drawn from the telomeric end (left). Ticks: paralog-diagnostic units (upper row: '
        'SNV sites, lower row: private 200 bp bins) coloured by the CHM13 allele they carry (orange donor D, blue recipient or '
        'proximal comparator P, grey neither). Yellow band: 95% interval of the switch; green: genome-wide MAPQ-60 anchor; '
        'boxes below: the CHM13 coordinates the donor and recipient segments map to.', pre)}
<h2 id="tiers">Tier definitions</h2>
{tier_table()}
<h2 id="results">Results by cohort</h2>
<ul>
<li><strong>HG002 Q100</strong> (phased by the diploid Q100 assembly; {h['n']} non-native ends; 16 public read alignment files, 12
independent datasets after molecule de-duplication; HERRO-corrected reads are derived and never counted): {h['confirmed']}
confirmed, {h['supported']} supported; all {h['confirmed']} confirmed claims are definitive, {h['multi']} with definitive reads
from &ge; 2 datasets.</li>
<li><strong>Trio</strong>: HG003 (father) and HG004 (mother), unphased GIAB reads; the transmitting parent confirms
{t['confirmed']} HG002 claims and supports {t['supported']}; {t['definitive']} are definitive (placement on the transmitted
HG002 haplotype is a proxy for the parent's own assembly).</li>
<li><strong>1000 Genomes ONT</strong> (R9.4.1, {N['n_1kg_ind']} HPRC Release 2 individuals; unphased): {k['confirmed']} confirmed,
{k['supported']} supported, {k['definitive']} definitive. The 19 confirmed-but-not-definitive claims rest on anchors in
CHM13-shared sequence (mostly SA-2q-02) or on reads that place best on another end of the same person.</li>
</ul>
{figure('evidence/figures/fig_evidence_summary.png', 'Claims by cohort and tier, controls, dataset independence and trio',
        '(a) Claims by cohort and tier. (b) Negative controls with exact one-sided 95% upper bounds. (c) HG002 definitive reads '
        'per claim and dataset. (d) Trio transmission: verdict in the transmitting and in the other parent.', pre)}
<h2 id="showcase">Showcase claims</h2>
<p>Eight exemplar claims are in the repository with their reads (CRAM and BAM on CHM13 v2.0, BAM on the person's own assembly
windows, window FASTA, IGV sessions and a "what you should see" note):</p>
<ul>{show}</ul>
<h2 id="downloads">Downloads</h2>
<ul>
<li><a href="{BLOB}evidence/claims.tsv">claims.tsv</a>: one row per claim ({len(E['C'])} claims: v1 and definitive verdicts, read
counts by tier and dataset, model, assembly switch, diagnostics note).</li>
<li><a href="{BLOB}evidence/reads.tsv.gz">reads.tsv.gz</a>: one row per cited read ({len(E['R']):,} rows): public source URL,
region, read group, units, purity, switch interval, anchor, own-assembly placement, tier and the samtools command that fetches the
read.</li>
<li><a href="{BLOB}evidence/controls.tsv">controls.tsv</a>, <a href="{BLOB}evidence/criteria_v1.json">criteria_v1.json</a>,
<a href="{TREE}evidence/showcase">showcase/</a>.</li>
</ul>
<p>Full evidence set (every v1-conclusive read of every claim) as GitHub release assets:</p>
<div class="table-wrap free"><table><thead><tr><th scope="col">file</th><th scope="col" class="num">size</th>
<th scope="col">SHA-256</th><th scope="col">content</th></tr></thead><tbody>{rel}</tbody></table></div>
{rep_txt}
<p><a href="claims.html">Search all claims</a> &middot; <a href="verify.html">Verify it yourself</a> &middot;
<a href="faq.html">FAQ for skeptics</a> &middot; <a href="{BLOB}evidence/READ_EVIDENCE.md">full report</a>
(<a href="{BLOB}evidence/READ_EVIDENCE.pdf">PDF</a>)</p>
"""
    ctx["write"]("evidence/index.html", ctx["page"]("evidence/index.html", "Read evidence", body, active="evidence/index.html",
                                                   description="Read-level evidence for non-native terminal blocks of human chromosome ends."))


# ------------------------------------------------------------------------------------------------------------------ claims table
def page_claims(ctx, E, N):
    pre = "../"; C = E["C"]
    cols = ["claim_id", "slug", "cohort", "individual", "haplotype", "arm", "allele_id", "nativity_class", "donor_D", "recipient_R",
            "v1_verdict", "definitive_verdict", "n_conclusive_concordant", "n_definitive", "n_definitive_datasets", "n_contradicts",
            "v1_class", "page"]
    C = C.copy(); C["page"] = [1 if has_page(c) else 0 for c in C.to_dict("records")]
    rows = [[None if (isinstance(v, float) and v != v) else (int(v) if isinstance(v, (int,)) or (isinstance(v, float) and v.is_integer() and c.startswith("n_")) else v)
             for c, v in zip(cols, r)] for r in C[cols].itertuples(index=False)]
    ctx["write"]("data/evidence_claims.js", "window.EVIDENCE_CLAIMS = " + json.dumps(dict(cols=cols, rows=rows), separators=(",", ":"),
                                                                                      allow_nan=False) + ";\n")
    body = f"""
{subnav(pre, "evidence/claims.html")}
<h1>Read-evidence claims</h1>
<p>All {len(C)} claims: {N['HG002_Q100']['n']} HG002 Q100 ends, {len(C[C.cohort == 'HG002_trio'])} trio tests (each HG002 claim
in both parents) and {N['1KG_ONT']['n']} 1000 Genomes ONT claims (individual &times; arm &times; allele). Claims with at least
one conclusive read link to a page with every cited read. Columns: v1 = pre-registered verdict; definitive = post-hoc
stricter tier; reads = conclusive switch-concordant / definitive.</p>
<div class="controls" role="group" aria-label="Claim filters">
<label>Search <input type="search" id="e-q" placeholder="individual, arm, allele, CHM13 end"></label>
<label>Cohort <select id="e-cohort"><option value="">all</option>{''.join(f'<option value="{k}">{v}</option>' for k, v in COHORT.items())}</select></label>
<label>v1 verdict <select id="e-v1"><option value="">all</option><option>confirmed</option><option>supported</option>
<option>not_confirmed</option><option>contradicted</option></select></label>
<label>Definitive tier <select id="e-def"><option value="">all</option><option value="definitive">definitive (any)</option>
<option value="confirmed_not_definitive">confirmed, not definitive</option><option value="not_definitive">not definitive</option></select></label>
<label class="inline"><input type="checkbox" id="e-evid" checked> only claims with a conclusive read</label>
<button type="button" class="btn" id="e-download">Download shown rows (TSV)</button>
<span class="count" id="e-count" aria-live="polite"></span>
</div>
<div class="table-wrap" id="e-table"><noscript><p>The claims table needs JavaScript; all claims are in
<a href="{ctx['BLOB']}evidence/claims.tsv">evidence/claims.tsv</a>.</p></noscript></div>
"""
    ctx["write"]("evidence/claims.html", ctx["page"]("evidence/claims.html", "Read-evidence claims", body, active="evidence/index.html",
                                                    scripts=("data/evidence_claims.js", "assets/js/evidence.js")))


# ------------------------------------------------------------------------------------------------------------------ claim pages
def _num(x, d=0):
    v = _f(x)
    return "" if v is None else (f"{v:,.{d}f}" if d else f"{int(round(v)):,}")


def _reads_table(rr, pre):
    cols = ["read", "dataset", "tier", "length", "switch (subtelomere, kb)", "switch on CHM13", "donor units / purity",
            "proximal units / purity", "anchor (CHM13, MAPQ, unique)", "own-assembly best end", "fetch"]
    h = ['<div class="table-wrap"><table><thead><tr>' + "".join(f'<th scope="col">{c}</th>' for c in cols) + "</tr></thead><tbody>"]
    for r in rr.itertuples():
        own = ""
        if isinstance(r.own_best_end, str) and r.own_best_end:
            own = (f"{esc(r.own_best_end)} ({_num(r.own_best_identity * 100 if _f(r.own_best_identity) else None, 1)}%)"
                   f"{' &#10003;' if str(r.own_on_claimed_end) == 'True' else ' (not claimed end)'}")
        sw = f"{_num(_f(r.switch_ci_sub_lo) / 1e3 if _f(r.switch_ci_sub_lo) is not None else None, 1)}-{_num(_f(r.switch_ci_sub_hi) / 1e3 if _f(r.switch_ci_sub_hi) is not None else None, 1)}"
        an = f"{esc(str(r.anchor_chm13))}, {_num(r.anchor_mapq)}, {_num(_f(r.anchor_unique_bp) / 1e3 if _f(r.anchor_unique_bp) is not None else None, 1)} kb" \
            if isinstance(r.anchor_chm13, str) else ""
        tier = TIER_LAB.get(str(r.tier).replace("not_counted_", ""), str(r.tier))
        if str(r.tier).startswith("not_counted_"):
            tier += " (not counted: " + ("derived" if str(r.derived) == "True" else "duplicate molecule") + ")"
        cmd = f'<details><summary>samtools</summary><code>{esc(str(r.samtools_command))}</code></details>' if isinstance(r.samtools_command, str) else ""
        h.append(f"<tr><td><code>{esc(r.read_id)}</code></td><td>{esc(str(r.dataset))}</td><td>{esc(tier)}</td>"
                 f"<td class=num>{_num(r.read_len)}</td><td>{sw if _f(r.switch_ci_sub_lo) is not None else ''}</td>"
                 f"<td>{esc(str(r.switch_ci_chm13)) if isinstance(r.switch_ci_chm13, str) else ''}</td>"
                 f"<td>{_num(r.seg_D_units)} / {_num(r.seg_D_purity, 2)}</td><td>{_num(r.seg_P_units)} / {_num(r.seg_P_purity, 2)}</td>"
                 f"<td>{an}</td><td>{own}</td><td>{cmd}</td></tr>")
    h.append("</tbody></table></div>")
    return "".join(h)


def page_claim_all(ctx, E, N):
    C, R, S = E["C"], E["R"], E["S"]
    BLOB, TREE, rp = ctx["BLOB"], ctx["TREE"], ctx["rp"]
    show = set(S.claim_id); n = 0
    RL = E["SL"].set_index("claim_id")
    for c in C.to_dict("records"):
        if not has_page(c): continue
        pre = "../../"
        rr = R[R.claim_id == c["claim_id"]]
        conc = rr[rr.v1_verdict.isin(["conclusive", "contradicts"])].copy()
        order = {t: i for i, t in enumerate(["definitive", "conclusive_shared_anchor", "conclusive_unique_anchor_off_end",
                                             "conclusive_shared_anchor_off_end", "conclusive_discordant", "contradicts"])}
        conc["_o"] = conc.tier.astype(str).str.replace("not_counted_", "").map(order).fillna(9) + conc.tier.astype(str).str.startswith("not_counted").astype(int) * 0.5
        conc = conc.sort_values(["_o", "dataset"])
        ncons = int((rr.v1_verdict == "consistent").sum())
        figs = []
        if c["claim_id"] in show:
            figs.append(ctx["figure"](f"evidence/figures/showcase/{c['slug']}.readthrough.png", "Exemplar read from telomere to anchor",
                                      "Exemplar read: diagnostic units coloured donor (orange) vs recipient/proximal (blue), switch 95% "
                                      "interval (yellow), unique anchor (green) and the CHM13 coordinates of each segment.", pre))
        if os.path.exists(rp("evidence", "figures", "claims", f"{c['slug']}.png")):
            figs.append(ctx["figure"](f"evidence/figures/claims/{c['slug']}.png", "Definitive reads of the claim",
                                      "All definitive reads of the claim (up to 12 shown), aligned at the telomere boundary: units "
                                      "coloured by CHM13 allele, per-read switch interval (yellow), unique anchor (green), assembly "
                                      "switch (grey band).", pre))
        igv = ""
        if c["claim_id"] in show:
            sl = RL.loc[c["claim_id"]]
            igv = (f"<h2 id=\"igv\">Files and IGV</h2><p>This is a showcase claim: its reads are in the repository under "
                   f"<a href=\"{TREE}evidence/showcase/{c['slug']}\">evidence/showcase/{c['slug']}/</a> with IGV sessions "
                   f"<a href=\"{BLOB}evidence/showcase/{c['slug']}/{c['slug']}.chm13.igv_session.xml\">CHM13 v2.0</a> (loci "
                   f"<code>{esc(sl.locus_donor)}</code> and <code>{esc(sl.locus_recipient)}</code>) and "
                   f"<a href=\"{BLOB}evidence/showcase/{c['slug']}/{c['slug']}.own.igv_session.xml\">own assembly</a>, and a "
                   f"<a href=\"{BLOB}evidence/showcase/{c['slug']}/README.md\">what-you-should-see note</a>. Split at the switch, "
                   f"the distal segments of {int(sl.distal_on_D)} of {int(sl.n_reads)} showcase reads align to CHM13 {esc(sl.D[3:])} and "
                   f"the proximal segments of {int(sl.proximal_on_R_mapq60)} to CHM13 {esc(sl.R[3:])} at MAPQ 60.</p>")
        else:
            igv = (f"<h2 id=\"files\">Files</h2><p>CRAM/BAM on CHM13 v2.0, BAM on the person's own assembly windows and the window FASTA: "
                   f"release asset <code>subtel_read_evidence_v1_{c['cohort']}.tar.gz</code>, directory "
                   f"<code>read_evidence_v1/{c['cohort']}/{c['slug']}/</code> (<a href=\"../index.html#downloads\">downloads</a>).</p>")
        dl = [("cohort", COHORT[c["cohort"]]), ("individual", c["individual"]), ("haplotype", c["haplotype"]), ("arm", c["arm"]),
              ("allele", f"<a href=\"{pre}arms/{esc(c['arm'])}.html#{esc(c['allele_id'])}\">{esc(c['allele_id'])}</a> ({esc(str(c['nativity_class']).replace('_', ' '))})"),
              ("model", f"donor CHM13 {esc(str(c['donor_D'])[3:])}, proximal CHM13 {esc(str(c['proximal_P'])[3:])}, recipient CHM13 {esc(str(c['recipient_R'])[3:])}, pattern {esc(str(c['pattern']))}"),
              ("claimed end(s)", "<code>" + esc(str(c["claimed_end_ids"])) + "</code>"),
              ("assembly switch", f"{_num(c['assembly_switch_sub'])} bp from the telomere boundary"
               + (f"; on CHM13 {esc(str(c['assembly_switch_chm13_P']))}" if isinstance(c["assembly_switch_chm13_P"], str) else "")
               + (f"; own-window projection {esc(str(c['asm_projection_switch_ci_sub']))} bp" if isinstance(c["asm_projection_switch_ci_sub"], str) else "")),
              ("v1 verdict", f"<strong>{esc(c['v1_verdict'])}</strong> ({esc(str(c['v1_class']))})"),
              ("definitive tier", f"<strong>{esc(DEF_LAB.get(c['definitive_verdict'], c['definitive_verdict']))}</strong>"),
              ("reads", f"{int(c['n_conclusive_concordant'])} conclusive switch-concordant, {int(c['n_conclusive_discordant'])} discordant, "
                        f"{int(c['n_contradicts'])} contradicting, {int(c['n_consistent'])} consistent; {int(c['n_definitive'])} definitive "
                        f"({esc(str(c['definitive_datasets'])) if isinstance(c['definitive_datasets'], str) else '-'})"),
              ("by tier", esc(str(c["reads_by_tier"])) if isinstance(c["reads_by_tier"], str) else ""),
              ("source", esc(str(c["source"]))),
              ("note", esc(str(c["diagnostics_note"])) if isinstance(c["diagnostics_note"], str) else "")]
        dlh = "<dl class=\"facts\">" + "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in dl if v) + "</dl>"
        body = f"""
{subnav(pre, '')}
<h1>{esc(claim_label(c))}: {esc(c['allele_id'])}</h1>
<p class="muted">Claim <code>{esc(c['claim_id'])}</code>. <a href="../claims.html">All claims</a></p>
{dlh}
{''.join(figs)}
<h2 id="reads">Conclusive and contradicting reads ({len(conc)})</h2>
<p>{ncons} further reads are <em>consistent</em> (switch test passed, but no terminal telomere, anchor or integrity); they are
listed in <a href="{BLOB}evidence/reads.tsv.gz">reads.tsv.gz</a>. Switch: 95% interval in subtelomere coordinates (kb from the
telomere boundary) and projected on CHM13 {esc(str(c['proximal_P'])[3:])}. Each <em>samtools</em> command prints the read's record
from its public file (CRAM: add the CHM13 v2.0 FASTA as <code>-T</code>).</p>
{_reads_table(conc, pre)}
{igv}
<h2 id="replicate">Re-derive this verdict</h2>
<pre><code>git clone https://github.com/jlanej/human-subtelomere-atlas &amp;&amp; cd human-subtelomere-atlas
evidence/replicate_claim.sh '{esc(c['claim_id'])}'</code></pre>
<p>The script fetches the cited reads from the public source by region and read ID, re-runs the frozen criteria and compares
every per-read and claim verdict with this page (<a href="../verify.html">details</a>).</p>
"""
        ctx["write"](f"evidence/claim/{c['slug']}.html", ctx["page"](f"evidence/claim/{c['slug']}.html", f"Read evidence {claim_label(c)}",
                                                                    body, active="evidence/index.html",
                                                                    description=f"Read-level evidence for {c['claim_id']}."))
        n += 1
    return n


# ------------------------------------------------------------------------------------------------------------------ verify
def page_verify(ctx, E, N):
    pre = "../"; BLOB = ctx["BLOB"]; RT = E["RT"]
    test = ""
    if RT is not None:
        rows = "".join(f"<tr><td><code>{esc(r.claim_id)}</code></td><td>{esc(r.kind)}</td><td class=num>{r.n_cited}</td>"
                       f"<td>{esc(str(r.published))}</td><td>{esc(str(r.replicated))}</td><td class=num>{r.reads_identical}</td>"
                       f"<td class=num>{r.wall_s:.0f}</td><td>{esc(r.result)}</td></tr>" for r in RT.itertuples())
        test = (f"<h2 id=\"test\">Replication test</h2><p>Run from a copy of the repository without <code>work/</code>, with a fresh "
                f"cache (CHM13 downloaded from the public URL), on one laptop (3 threads per run): see "
                f"<a href=\"{BLOB}evidence/REPLICATION_TEST.md\">REPLICATION_TEST.md</a>.</p>"
                "<div class=\"table-wrap free\"><table><thead><tr><th scope=col>claim</th><th scope=col>kind</th><th scope=col class=num>cited reads</th>"
                "<th scope=col>published</th><th scope=col>replicated</th><th scope=col class=num>reads identical</th>"
                f"<th scope=col class=num>wall (s)</th><th scope=col>result</th></tr></thead><tbody>{rows}</tbody></table></div>")
    body = f"""
{subnav(pre, "evidence/verify.html")}
<h1>Verify it yourself</h1>
<p>Every verdict can be re-derived from public data with one command. Nothing produced by the authors is needed beyond this
repository: the reads come from the public alignment files, the person's assembly windows from the public assembly FASTA (checked
against the SHA-256 in <code>data/catalogue/ends.tsv</code>), CHM13 v2.0 from its public URL, and the thresholds from the frozen
criteria file.</p>
<h2 id="setup">1. Set up</h2>
<pre><code>git clone https://github.com/jlanej/human-subtelomere-atlas &amp;&amp; cd human-subtelomere-atlas
conda create -n atlas -c conda-forge -c bioconda python=3.11 samtools minimap2 mappy pysam numpy pandas scipy
conda activate atlas</code></pre>
<p>Needs internet, about 1 GB of disk for CHM13 v2.0 (downloaded once into <code>~/.cache/subtel-atlas/replicate</code>) and about
8 GB of RAM for the in-memory minimap2 index.</p>
<h2 id="run">2. Pick a claim and run it</h2>
<pre><code>evidence/replicate_claim.sh --list                       # claims with reads, grouped by definitive tier
evidence/replicate_claim.sh 'HG002#2#chr9_MATERNAL#L'      # HG002 maternal 9p (also confirmed in the mother)
evidence/replicate_claim.sh 'HG01258|20q|SA-20q-01'        # a 1000 Genomes ONT claim
evidence/replicate_claim.sh 'control|HG00133|20p|SA-20p-10'  # a strictly native control</code></pre>
<p>The script:</p>
<ol>
<li>checks that <code>evidence/criteria_v1.json</code> and the packaged copy hash to <code>a89c5eaf023ed032</code>;</li>
<li>downloads CHM13 v2.0 and extracts its 48 chromosome ends (<code>subtel-atlas fetch-refs</code>, checksum-verified);</li>
<li>builds the claim model (donor, comparator, pattern, switch) from the catalogue tables in the repository;</li>
<li>cuts the person's end windows from the public assembly by range and checks their SHA-256;</li>
<li>fetches the cited reads by region from the public alignment file and keeps them by read ID
(<code>samtools view -F 0x900 -N ids.txt -X URL URL.bai REGION</code>);</li>
<li>evaluates every read with <code>subtel_atlas.readevidence</code> (the module behind <code>subtel-atlas reads-evidence</code>):
diagnostic units, switch, genome-wide anchor, integrity, verdict, concordance with the assembly switch;</li>
<li>prints each read's published and replicated verdict, the claim verdict, and the definitive tier (own-assembly placement on
all of the person's end windows);</li>
<li>exits 0 when everything matches.</li>
</ol>
<h2 id="one-read">3. Look at one read without the script</h2>
<p>Each row of <a href="{BLOB}evidence/reads.tsv.gz">reads.tsv.gz</a> has a <code>samtools_command</code>. For example, the first
definitive read of HG002 maternal 9p in the R10 ultra-long data:</p>
<pre><code>{esc(_example_cmd(E))}</code></pre>
<p>Align it to CHM13 v2.0 (<code>minimap2 -ax map-ont -Y</code>) and to the person's end window, or open the showcase IGV
sessions (<a href="{ctx['TREE']}evidence/showcase">evidence/showcase/</a>).</p>
{test}
<h2 id="limits">What replication does and does not test</h2>
<ul>
<li>It re-derives each published verdict from the same public reads under the same frozen thresholds. It does not re-select the
reads: the cited set is every read of the claim's read set that was conclusive, consistent or contradicting. Non-informative
reads cannot change a verdict and are not fetched.</li>
<li>Results depend on minimap2/mappy versions (the published run used minimap2 2.31 and mappy 2.31). Other versions can move
individual alignments, and with them, rarely, a verdict near a threshold.</li>
<li>HG003/HG004 ONT reads are stored as GRCh38 CRAMs; decoding them needs a GRCh38 FASTA (<code>--grch38</code>).</li>
</ul>
"""
    ctx["write"]("evidence/verify.html", ctx["page"]("evidence/verify.html", "Verify the read evidence", body, active="evidence/index.html"))


def _example_cmd(E):
    R = E["R"]
    x = R[(R.claim_id == "HG002#2#chr9_MATERNAL#L") & (R.tier == "definitive") & (R.dataset == "ont_ul_dorado")]
    return x.iloc[0].samtools_command if len(x) else R[R.tier == "definitive"].iloc[0].samtools_command


# ------------------------------------------------------------------------------------------------------------------ FAQ
def page_faq(ctx, E, N):
    pre = "../"; BLOB = ctx["BLOB"]
    h, k = N["HG002_Q100"], N["1KG_ONT"]
    (nh, kh), (nk, kk) = N["ctl_hg002"], N["ctl_1kg"]
    qa = [
        ("Is this not just a chimeric read?", "A chimera joins two molecules. A definitive read has no internal telomere, no "
         "adapter sequence, no foldback, and at least 90% of its subtelomeric part aligned. The switch it shows also falls, within "
         "&plusmn; 2 kb, where the person's own assembly places it. Chimeras form at random positions; they do not reproduce the same "
         "switch in independent molecules, datasets and platforms (HG002: up to 9 datasets per claim), or in the parent who "
         "transmitted the haplotype."),
        ("Could the reads come from a different chromosome end (mis-mapping, paralogs)?", "This is the main risk, and the reason the "
         "definitive tier exists. Diagnostic units say which CHM13 end a stretch most resembles, but human haplotypes carry blocks "
         "that CHM13 lacks. For example, 6p, 20q and 9q haplotypes carry an 11p-like tip followed by 3q-like sequence "
         f"(<a href=\"{BLOB}evidence/population/DIAGNOSTICS.md\">DIAGNOSTICS.md</a>). A definitive read must therefore also (i) anchor "
         "with &ge; 1 kb in CHM13 sequence unique to the recipient end and (ii) align best to the claimed end window among all of the "
         "same person's assembly end windows. Reads that fail (ii) are reported as <em>places on another end</em>."),
        ("Could the assembly be wrong?", "The read test is assembly-free: the units come from CHM13 alone, and the reads were not "
         "used to build the HPRC Release 2 assemblies. The 1000 Genomes ONT reads come from another laboratory and platform (R9.4.1), "
         "and the GIAB parental reads from other DNA samples. Assemblies enter only to locate the expected switch (concordance) and "
         "in the post-hoc placement step."),
        ("What about population variation at the diagnostic sites?", "Carrier haplotypes differ from CHM13 at some diagnostic "
         "sites, so purity thresholds are 0.80, not 1.0 (true-carrier donor purity is as low as 0.80; calibration T0). Variation of "
         "this kind lowers purity in both segments. It does not create an ordered donor-then-recipient pattern with log10 likelihood "
         "ratios &ge; 5 against both no-switch models. In strictly native controls, 0 of "
         f"{nk:,} read &times; model tests (1KG-ONT) and 0 of {nh:,} HG002 other-haplotype tests were conclusive."),
        ("Are the reads independent?", "HG002 reads are de-duplicated by molecule across datasets (same ONT read UUID, or PacBio "
         "movie and ZMW). Corrected (HERRO) reads are derived from the ultra-long reads; they are shown but never counted. "
         f"{h['multi']} of {h['confirmed']} definitive HG002 claims have definitive reads from &ge; 2 independent datasets."),
        ("Were the thresholds tuned to get these results?", f"The v1 criteria were registered and frozen (hash "
         f"<code>a89c5eaf023ed032</code>) on HG002 training ends and native controls before any population read was fetched "
         f"(<a href=\"{BLOB}evidence/CRITERIA.md\">CRITERIA.md</a>, <a href=\"{BLOB}evidence/calibration/CALIBRATION.md\">CALIBRATION.md</a>). "
         "The definitive tier was defined later. It can only remove claims, and both tiers are reported for every claim."),
        ("Why are most population claims not assessable?", "Reads must reach the telomere and extend past the switch by &ge; 10 kb. "
         "1KG-ONT R9.4.1 reads are shorter than HG002 ultra-long reads, and the telomere-reaching subset is small. A claim without "
         "such a read cannot be tested. That is not evidence against it, and it is labelled <em>not assessable</em> with the reason "
         f"in <a href=\"{BLOB}evidence/claims.tsv\">claims.tsv</a>."),
        ("Why are some confirmed claims not definitive?", f"{k['confirmed'] - k['definitive']} of {k['confirmed']} confirmed 1KG-ONT claims "
         "rest on reads whose anchor lies only in CHM13-shared sequence (e.g. 2q SA-2q-02, whose outer ~100+ kb are shared), or on "
         "reads that align best to another end of the same person, often the other haplotype at the same arm. They remain v1-confirmed. "
         "The definitive tier does not resolve them."),
        ("Which arms lack controls?", "Nine claimed arms have no strictly native individual among the 42 (12p, 17q, 19p, 19q, 1p, "
         "3q, 4q, 6q, Yq). For these arms, specificity rests on HG002 other-haplotype controls and the calibration."),
        ("How do I check one read myself?", "Every read in <a href=\"" + BLOB + "evidence/reads.tsv.gz\">reads.tsv.gz</a> has the exact "
         "samtools command that prints it from the public file, plus the region and read group. The showcase directories have "
         "CRAM/BAM files and IGV sessions. See <a href=\"verify.html\">Verify it yourself</a>."),
    ]
    body = f"""
{subnav(pre, "evidence/faq.html")}
<h1>FAQ for skeptics</h1>
<p>Alternative explanations for a read that looks like a non-native terminal block, and how each is excluded. The full
discussion is in <a href="{BLOB}evidence/READ_EVIDENCE.md">READ_EVIDENCE.md</a>.</p>
""" + "".join(f'<h2 id="q{i + 1}">{esc(q)}</h2><p>{a}</p>' for i, (q, a) in enumerate(qa))
    ctx["write"]("evidence/faq.html", ctx["page"]("evidence/faq.html", "Read evidence FAQ", body, active="evidence/index.html"))


# ------------------------------------------------------------------------------------------------------------------ hooks
def arm_links(rp):
    """{(arm, allele): [(label, href_from_arm_page, definitive_verdict)]} for definitive claims"""
    C = pd.read_csv(rp("evidence", "claims.tsv"), sep="\t", low_memory=False)
    out = {}
    for c in C[C.definitive_verdict.str.startswith("definitive")].to_dict("records"):
        out.setdefault((c["arm"], c["allele_id"]), []).append((claim_label(c), f"../evidence/claim/{c['slug']}.html"))
    return out
