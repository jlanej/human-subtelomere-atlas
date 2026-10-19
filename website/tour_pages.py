"""Guided read tour (docs/evidence/tour.html): igv.js views of the read evidence, from docs/evidence/tour/data/tour_data.json
(written by website/tour_build.py). Writes evidence/tour.html and evidence/tour/tour.json (chapters for assets/js/tour.js)."""
import json, os
from html import escape as esc

D = 'tour/data/'          # relative to docs/evidence/
COL = dict(donor='#D55E00', recip='#0072B2')


def _bam(url, name, height=90, mode='SQUISHED', **kw):
    t = dict(type='alignment', format='bam', url=url, indexURL=url + '.bai', name=name, displayMode=mode, height=height,
             visibilityWindow=1_200_000, showCoverage=False, showSoftClips=False, samplingDepth=50, order=kw.pop('order', None))
    t.update(kw); return {k: v for k, v in t.items() if v is not None}


def _bed(url, name, height=34, order=None):
    return {k: v for k, v in dict(type='annotation', format='bed', url=url, name=name, displayMode='EXPANDED', height=height,
                                  order=order).items() if v is not None}


def _ref(slug, label):
    b = f'{D}claims/{slug}/combo.fa.gz'
    return dict(id=slug, name=label, fastaURL=b, indexURL=b + '.fai', compressedIndexURL=b + '.gzi')


def _loc(contig, a, b): return f'{contig}:{max(1, int(a)):,}-{int(b):,}'.replace(',', '')


def _views(c, ref):
    """Everything / reference / own assembly loci for one claim"""
    wins = [w for w, _ in c['combo']['contigs'] if not w.startswith('CHM13_')]
    cw = c['claimed_window']; others = [w for w in wins if w != cw]
    sw = max(m['switch_hi'] for m in c['mols'])
    d_pieces = [p for p in c['combo']['asm_pieces'] if p.get('placed') == ref['donor_contig']]
    pa = [m for m in c['parts_alignments'] if m['ref'] == ref['donor_contig']]
    d_lo = min([p['start'] for p in pa] + [1]); d_hi = max([p['end'] for p in pa] + [sw + 5000])
    if d_hi - d_lo > 80_000: d_lo = max(1, d_hi - 60_000)
    r_hi = min(ref['recipient_len'], max(sw * 2 + 30_000, 80_000))
    dloc = _loc(ref['donor_contig'], d_lo - 2000, d_hi + 3000); rloc = _loc(ref['recipient_contig'], 1, r_hi)
    wloc = [_loc(cw, 1, 100_000)] + [_loc(w, 1, 100_000) for w in others[:1]]
    return [dict(label='Everything at once', locus=wloc + [dloc, rloc]),
            dict(label='Reference genome (CHM13)', locus=[dloc, rloc]),
            dict(label='Own assembly (both copies)', locus=wloc)]


def _model(c): return _bed(f'{D}claims/{c["slug"]}/model.bed', 'diagnostic sites (claim model)', height=40, order=0)


def _claim_tracks(c, ref, mol_name='molecules'):
    s = c['slug']; b = f'{D}claims/{s}/'
    return [_model(c), _bed(b + 'own_blocks.bed', 'assembly blocks (catalogue)', order=1),
            _bam(b + 'asm.bam', 'assembly blocks \u2192 CHM13', height=60, mode='EXPANDED', order=3),
            _bam(b + 'own.bam', f'{mol_name} \u2192 own assembly', height=110, order=4),
            _bam(b + 'parts.bam', f'{mol_name} \u2192 CHM13 (cut at switch)', height=110, order=5)]


def _label(c):
    hap = {'1': 'paternal', '2': 'maternal'}.get(str(c.get('haplotype', '')), '')
    who = c['individual'] + (f' {hap}' if hap and c['individual'] in ('HG002', 'HG003', 'HG004') else '')
    return f'{who} {c["arm"]}, {c["allele"]}'


def _ident(c, which):
    p = [x for x in c['combo']['asm_pieces'] if x.get('placed') == which]
    return max(p, key=lambda x: x['win_end'] - x['win_start']) if p else None


def _evidence(c, part):
    """plain-language summary of one part's agreement with the claim model"""
    t = c['site_totals'][part]; m = c['model']
    end = m['D'] if part == 'part1' else m['R']; other = m['P'] if part == 'part1' else m['D']
    if t['n']:
        return f'{end}\'s letter at {t["own"]:,} of {t["n"]:,} diagnostic sites ({other}\'s at {t["other"]:,})'
    if t['private_bp']:
        n = sum(1 for v in c['site_calls'].values() if v.get(part, {}).get('private_bp', 0) > 0)
        return f'sequence found only in the {end} end (not in {other}) in {n} of {len(c["mols"])} molecules'
    if part == 'part2' and not m['p_is_r']:
        return f'(its model compares against {m["P"]}, not shown here)'
    return 'no diagnostic sites in view'


def _on_claimed(c):
    return sum(1 for a in c['own_alignments'] if a['ref'] == c['claimed_window']), len(c['mols'])



def consolidated_caption(rp, html=True):
    """Caption of figures/pitch/fig_consolidated.png (report/pitch_consolidated.py), with counts from tour_data.json."""
    J = json.load(open(rp('docs/evidence/tour/data/tour_data.json'))); CL = {c['slug']: c for c in J['claims']}
    a, b = CL['HG002_2_chr3_MATERNAL_R'], CL['HG01258__20q__SA-20q-01']
    code = (lambda s: f'<code>{s}</code>') if html else (lambda s: s)
    return (f'Each panel stacks four tracks, all starting at the telomere: the CHM13 reference end the claim names as donor '
            f'(top), the person\'s own two assembled copies of the chromosome end (middle; claimed copy in bold) and the CHM13 '
            f'end of the chromosome itself (bottom). Shaded bands place each block of the person\'s assembly on the '
            f'reference end its catalogue label names (minimap2 asm10). Black lines are single DNA molecules from public '
            f'sequencing runs that were not used to build the assemblies: whole, on the claimed copy; cut at their switch '
            f'interval (criteria v1), part 1 aligned to the donor end and part 2 to the recipient end. '
            f'(a) HG002 chromosome 3q: {len(a["mols"])} molecules of the claim {code(a["claim_id"])}. '
            f'(b) HG01258 chromosome 20q, allele SA-20q-01: {len(b["mols"])} 1000 Genomes nanopore molecules of the claim '
            f'{code(b["claim_id"])}. Tracks are cut at 100 kb. The interactive version of every definitive claim is the '
            f'guided tour.')

def build(ctx):
    rp = ctx['rp']; J = json.load(open(rp('docs/evidence/tour/data/tour_data.json')))
    R, CL = J['refs'], {c['slug']: c for c in J['claims']}
    ch = []
    hero = CL[J['hero']['slug']]; hr = R[hero['ref']]; hv = _views(hero, hr)
    pd_, pr_ = _ident(hero, hr['donor_contig']), _ident(hero, hr['recipient_contig'])
    hm = next(m for m in hero['mols'] if m['read_id'] == J['hero']['read'])
    k, n = _on_claimed(hero)
    ch.append(dict(id='everything', short='Everything at once', title='The whole argument on one screen',
                   ref=_ref(hero['slug'], _label(hero)), views=hv, tracks=_claim_tracks(hero, hr),
                   lead=('Four panels, left to right: HG002\'s own two chromosome 3 ends (his mother\'s copy, then his '
                         'father\'s), then the reference genome\'s chromosome 11 end and chromosome 3 end. Every panel starts '
                         'at the telomere.'),
                   notes=[f'<b>Molecules on his own assembly</b>: {k} of {n} nanopore molecules (from '
                          f'{len({m["dataset"] for m in hero["mols"]})} datasets) align to the mother\'s copy and run from '
                          f'the telomere straight through the junction.',
                          f'<b>The assembly on the reference</b>: the mother\'s copy starts with a chromosome-11-like block '
                          f'({pd_["win_end"] / 1e3:.0f} kb; {pd_["ident_placed"]:.1%} identical to CHM13 chr11p, '
                          f'{pd_["ident_other"]:.1%} to chr3q), then chromosome 3 ({pr_["ident_placed"]:.1%} identical).',
                          f'<b>The molecules on the reference</b>, cut at their switch: part 1 carries '
                          f'{_evidence(hero, "part1")}; part 2 carries {_evidence(hero, "part2")}.',
                          'The father\'s copy (right of the mother\'s) carries a different allele, SA-3q-04.']))
    sc = hero['site_calls'].get(hm['read_id'], {})
    p1, p2 = sc.get('part1', {}), sc.get('part2', {})
    ch.append(dict(id='one-molecule', short='One molecule', title=f'One molecule, {hm["read_len"] / 1e3:.0f} kb, start to finish',
                   ref=_ref(hero['slug'], _label(hero)), views=[hv[1]],
                   tracks=[_model(hero),
                           _bam(f'{D}extra/hero_parts.bam', 'the molecule, cut at its switch', height=60, mode='EXPANDED', order=2),
                           _bam(f'{D}extra/hero_whole_donor.bam', 'the whole molecule vs chr11 only', height=50, mode='EXPANDED', order=3),
                           _bam(f'{D}extra/hero_whole_recipient.bam', 'the whole molecule vs chr3 only', height=50, mode='EXPANDED', order=4)],
                   lead=(f'A single nanopore read ({esc(hm["dataset"])}, <code>{esc(hm["read_id"])}</code>): '
                         f'{hm["telo_len"] / 1e3:.1f} kb of telomere, then a switch at {hm["switch_lo"] / 1e3:.1f}-'
                         f'{hm["switch_hi"] / 1e3:.1f} kb.'),
                   notes=[f'Top track: cut at the switch. Part 1, aligned to the chromosome 11 end, carries chr11\'s letter at '
                          f'{p1.get("own", 0)} of {p1.get("n", 0)} diagnostic sites; part 2, aligned to the chromosome 3 end, '
                          f'carries chr3\'s letter at {p2.get("own", 0)} of {p2.get("n", 0)}. Zoom in on the orange and blue '
                          f'ticks of the claim-model track to see single sites.',
                          'Bottom two tracks: the same molecule, uncut, aligned to each end alone (primary alignment). '
                          'Each end accounts for only its own part of the molecule.']))
    pat = CL.get('HG002_1_chr3_PATERNAL_R')
    if pat:
        km, nm_ = _on_claimed(hero); kp, np_ = _on_claimed(pat)
        ch.append(dict(id='two-copies', short='Two copies', title='Two copies of chromosome 3, two different ends',
                       ref=_ref(hero['slug'], _label(hero)), views=[hv[2]],
                       tracks=[_bed(f'{D}claims/{hero["slug"]}/own_blocks.bed', 'assembly blocks (catalogue)', order=1),
                               _bam(f'{D}claims/{hero["slug"]}/own.bam', 'molecules of the maternal claim (SA-3q-01)', height=110, order=2),
                               _bam(f'{D}claims/{pat["slug"]}/own.bam', 'molecules of the paternal claim (SA-3q-04)', height=110, order=3)],
                       lead=('HG002 inherited two different chromosome 3 ends. Left: his mother\'s copy (a chromosome-11-like '
                             'block, then 3q). Right: his father\'s copy (a short chromosome-11-like tip block, then 19p-like blocks among 3q '
                             'sequence).'),
                       notes=[f'Each set of molecules lands on its own copy: {km} of {nm_} maternal-claim molecules on the '
                              f'left window, {kp} of {np_} paternal-claim molecules on the right one.',
                              'Both copies start with chromosome-11-like sequence next to the telomere, but it runs about 37 kb '
                              'on the mother\'s copy and about 6 kb on the father\'s, and what follows differs.',
                              'The built-in control: the same chromosome arm, in the same person, read by the same machines, '
                              'gives two different and internally consistent ends.']))
    kid, mom = CL.get('HG002_2_chr9_MATERNAL_L'), CL.get('HG004__HG002_2_chr9_MATERNAL_L')
    if kid and mom:
        v = _views(kid, R[kid['ref']])[1]
        sw_kid = ', '.join(f'{m["switch_lo"] / 1e3:.0f}' for m in kid['mols'])
        sw_mom = ', '.join(f'{m["switch_lo"] / 1e3:.0f}' for m in mom['mols'])
        ch.append(dict(id='inherited', short='Inherited', title='Passed from mother to son', ref=_ref(kid['slug'], _label(kid)),
                       views=[v],
                       tracks=[_model(kid), _bam(f'{D}claims/{kid["slug"]}/parts.bam', 'HG002 (son)', height=110, order=2),
                               _bam(f'{D}claims/{mom["slug"]}/parts.bam', 'HG004 (his mother)', height=110, order=3)],
                       lead=('HG002\'s maternal chromosome 9 end starts with chromosome-15-like DNA. Molecules from HG002 '
                             '(top) and from his mother HG004 (bottom), sequenced separately.'),
                       notes=[f'Son: part 1 carries {_evidence(kid, "part1")}; part 2 carries {_evidence(kid, "part2")}.',
                              f'Mother: part 1 carries {_evidence(mom, "part1")}; part 2 carries {_evidence(mom, "part2")}.',
                              f'Switch lower bounds: son {sw_kid} kb; mother {sw_mom} kb.',
                              'An assembly error is not inherited; a real chromosome end is.']))
    com = [c for c in J['claims'] if c['allele'] == 'SA-20q-01']
    if com:
        v = _views(com[0], R[com[0]['ref']])[1]
        nmol = sum(len(c['mols']) for c in com)
        n1 = sum(1 for c in com for x in c['site_calls'].values() if x.get('part1', {}).get('private_bp', 0) > 0)
        n2 = sum(1 for c in com for x in c['site_calls'].values() if x.get('part2', {}).get('private_bp', 0) > 0)
        ch.append(dict(id='common', short='Common', title=f'The same junction in {len(com)} unrelated people',
                       ref=_ref(com[0]['slug'], 'SA-20q-01'), views=[v],
                       tracks=[_model(com[0])] + [_bam(f'{D}claims/{c["slug"]}/parts.bam', c['individual'], height=45, order=2 + i)
                                                  for i, c in enumerate(com)],
                       lead=('SA-20q-01: the chromosome 20q end starts with a block of 3q/19p-type DNA. About half of all '
                             'assembled 20q ends carry it. One track per person; 1000 Genomes nanopore data (Schloissnig et '
                             'al. 2025), independent of the assemblies.'),
                       notes=[f'{len(com)} people, {nmol} molecules. The two CHM13 ends share no block here, so the evidence is '
                              f'sequence found in only one of them: part 1 covers 3q-only sequence in {n1} of {nmol} molecules, '
                              f'part 2 covers 20q-only sequence in {n2} of {nmol}.',
                              'These are older R9.4.1 reads, so mismatches are denser than in HG002. The switch sits in the same '
                              'region in everyone, because it is one inherited allele.']))
    tc = CL.get(J['tech']['slug'])
    if tc:
        v = _views(tc, R[tc['ref']])[1]; files = sorted(J['tech']['files'].items())
        per = []
        for m in tc['mols']:
            x = tc['site_calls'].get(m['read_id'], {}); a, b = x.get('part1', {}), x.get('part2', {})
            per.append(f'{m["dataset"].replace("_", " ")}: {a.get("own", 0)}/{a.get("n", 0)} and {b.get("own", 0)}/{b.get("n", 0)}')
        ch.append(dict(id='technologies', short='Any technology', title='Every sequencing technology agrees',
                       ref=_ref(tc['slug'], _label(tc)), views=[v],
                       tracks=[_model(tc)] + [_bam(f'{D}extra/{f}', ds.replace('_', ' '), height=40, order=2 + i)
                                              for i, (ds, f) in enumerate(files)],
                       lead=(f'HG002\'s maternal chromosome 19q end starts with chromosome-22-like DNA. One molecule from each '
                             f'of {len(files)} HG002 datasets: ' + ', '.join(ds.replace('_', ' ') for ds, _ in files) + '.'),
                       notes=['Diagnostic sites carrying the expected letter, part 1 (22q) and part 2 (19q): ' + '; '.join(per) + '.',
                              'Different chemistries, error profiles and labs; the same switch.']))
    s17 = [c for c in J['claims'] if c['allele'] == 'SA-17p-02']
    if s17:
        v = _views(s17[0], R[s17[0]['ref']])[1]
        nmol = sum(len(c['mols']) for c in s17)
        n1 = sum(1 for c in s17 for x in c['site_calls'].values() if x.get('part1', {}).get('private_bp', 0) > 0)
        ch.append(dict(id='deep-block', short='Deep block', title='A block from 200 kb inside another chromosome',
                       ref=_ref(s17[0]['slug'], 'SA-17p-02'), views=[v],
                       tracks=[_model(s17[0])] + [_bam(f'{D}claims/{c["slug"]}/parts.bam', _label(c).split(',')[0], height=55, order=1 + i)
                                                  for i, c in enumerate(s17)],
                       lead=('SA-17p-02 (about 1 in 5 assembled 17p ends): the 17p end starts with a copy of sequence that, in '
                             'the reference genome, lies about 200 kb inside the chromosome 11p end.'),
                       notes=[f'{len(s17)} people: ' + ', '.join(_label(c).split(',')[0] for c in s17) + f'. Part 1 covers '
                              f'11p-only sequence in {n1} of {nmol} molecules (left panel, 190-230 kb into 11p).',
                              'This block is rich in tandem repeats, so alignments look ragged.']))
    gal = []
    for c in J['claims']:
        r = R[c['ref']]
        gal.append(dict(slug=c['slug'], label=_label(c), cohort=c['cohort'], verdict=c['verdict'], n_definitive=c['n_definitive'],
                        n_shown=len(c['mols']), ref=_ref(c['slug'], _label(c)), views=_views(c, r), tracks=_claim_tracks(c, r),
                        donor=c['D'], recipient=c['R'], page=f'claim/{c["slug"]}.html',
                        note=f'Part 1 carries {_evidence(c, "part1")}; part 2 carries {_evidence(c, "part2")}.'))
    ch.append(dict(id='gallery', short='All claims', title=f'All {len(gal)} definitive claims', gallery=gal,
                   lead=('Every claim with at least two definitive molecules, in the same four-panel view: own assembly '
                         '(both copies), then the CHM13 donor and recipient ends.'),
                   notes=['Pick a claim, or use the arrows. Up to 6 molecules per claim are shown (one per dataset first); all '
                          'are in the release assets.']))
    os.makedirs(rp('docs/evidence/tour'), exist_ok=True)
    ctx['write']('evidence/tour/tour.json', json.dumps(dict(chapters=ch), separators=(',', ':')))
    body = f"""
<nav class="arm-nav" aria-label="Read evidence"><a href="index.html">Overview</a><a href="tour.html" aria-current=page>Guided tour</a><a href="claims.html">Claims table</a><a href="verify.html">Verify it yourself</a><a href="faq.html">FAQ for skeptics</a></nav>
<h1>Guided tour: see the molecules</h1>
<p class="lead">Real DNA molecules from public sequencing data, in a genome browser, with nothing to install. Use
<b>Next</b> (or the arrow keys). Each stop loads its own data; drag to pan, scroll or use + / &minus; to zoom. Prefer a static picture? <a href="#one-figure">The whole argument in one figure</a>.</p>
<div class="tour-bar"><button id="tour-prev" type="button" aria-label="previous stop">&larr; Previous</button>
<ol id="tour-steps" class="tour-steps"></ol><button id="tour-next" type="button" aria-label="next stop">Next &rarr;</button></div>
<section class="tour-card" aria-live="polite"><h2 id="tour-title"></h2><p id="tour-lead"></p><div id="tour-views" class="tour-views"></div>
<div id="tour-gallery" class="tour-gallery" hidden><label for="tour-claim">Claim</label> <select id="tour-claim"></select>
<button id="gal-prev" type="button">&larr;</button><button id="gal-next" type="button">&rarr;</button> <a id="gal-page" href="#">claim page</a></div>
<ul id="tour-notes" class="tour-notes"></ul></section>
<div id="tour-status" class="muted" role="status"></div>
<div id="igv" class="igv-host"></div>
<h2 id="one-figure">The whole argument in one figure</h2>
{ctx['figure']('figures/pitch/fig_consolidated.png', 'Molecules, assembly and reference genome tell the same story', consolidated_caption(ctx['rp']), '../')}
<details class="note"><summary>What the tracks are</summary><ul>
<li><b>assembly blocks (catalogue)</b>: the atlas's labels for the person's own chromosome end (telomere; blocks named
after the CHM13 end they match best; home = the chromosome's own sequence).</li>
<li><b>diagnostic sites (claim model)</b>: the pre-registered claim model (criteria v1): single-base sites where the donor and
comparator CHM13 ends differ (orange on the donor end, blue on the comparator end; the label gives both letters), sequence found
only in one of them, and their shared blocks in grey (<a href="faq.html">FAQ</a>, <a href="verify.html">verify</a>).</li>
<li><b>assembly blocks &rarr; CHM13</b>: each block of the person's assembly placed on the reference end its label names.</li>
<li><b>molecules &rarr; own assembly</b>: whole molecules aligned to the person's own chromosome ends.</li>
<li><b>molecules &rarr; CHM13 (cut at switch)</b>: each molecule cut at its switch interval (criteria v1); part 1
(telomere to switch) aligned to the donor end, part 2 to the recipient end. Names end in <code>/part1</code>, <code>/part2</code>.</li>
<li>All reference panels are in <b>tip frame</b>: position 1 is the telomeric tip; q-arm ends are shown reverse-complemented
so that every panel reads from the telomere inward. Data and the script that made them: <a href="{ctx['BLOB']}website/tour_build.py">website/tour_build.py</a>.</li>
</ul></details>
<noscript><p>The tour needs JavaScript. The same files are listed on each <a href="claims.html">claim page</a>.</p></noscript>
"""
    ctx['write']('evidence/tour.html', ctx['page']('evidence/tour.html', 'Guided read tour', body, active='evidence/index.html',
                                                    scripts=('assets/js/vendor/igv.min.js', 'assets/js/tour.js'),
                                                    description='Guided igv.js tour of single DNA molecules that span non-native human chromosome ends.'))
    return len(ch), len(gal)
