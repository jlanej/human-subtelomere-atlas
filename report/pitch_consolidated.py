"""'Everything at once': molecules -> own assembly (both copies), own assembly -> CHM13, molecules -> CHM13, in one figure.
From docs/evidence/tour/data/tour_data.json (website/tour_build.py). Writes figures/pitch/fig_consolidated.(png|pdf).
Usage: python report/pitch_consolidated.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DONOR, RECIP, TELO, GREY, OTHER = '#D55E00', '#0072B2', '#222222', '#9A9A9A', '#7a7a7a'
EXAMPLES = [('HG002_2_chr3_MATERNAL_R', 'HG002, chromosome 3 ends', 'his mother\'s copy', 'his father\'s copy'),
            ('HG01258__20q__SA-20q-01', 'HG01258 (1000 Genomes), chromosome 20q ends', 'the copy with the 3q-like block',
             'the other copy')]


def style():
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'savefig.dpi': 300, 'pdf.fonttype': 42, 'axes.titlelocation': 'left', 'axes.spines.top': False,
                         'axes.spines.right': False})


def panel(ax, c, R, title, lab1, lab2, xmax=100_000):
    ri = R[c['ref']]; D, Rr = ri['donor_contig'], ri['recipient_contig']
    wins = [w for w, _ in c['combo']['contigs'] if not w.startswith('CHM13_')]
    w1 = c['claimed_window']; w2 = next((w for w in wins if w != w1), None)
    y = {D: 4.0, w1: 2.6, w2: 1.4, Rr: 0.0}; h = 0.28; s = 1e3
    # tracks
    blocks = {}
    for l in open(os.path.join(REPO, 'docs', 'evidence', 'tour', 'data', 'claims', c['slug'], 'own_blocks.bed')):
        if l.startswith('track'): continue
        f = l.rstrip('\n').split('\t'); blocks.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3], tuple(int(x) / 255 for x in f[8].split(','))))
    for nm, col, lab in ((D, DONOR, f'reference genome: the {ri["D"]} end'), (Rr, RECIP, f'reference genome: the {ri["R"]} end')):
        ax.add_patch(Rectangle((0, y[nm] - h / 2), xmax / s, h, color=col, alpha=0.85, lw=0))
        ax.add_patch(Rectangle((0, y[nm] - h / 2), 0.6, h, color=TELO, lw=0))
        ax.text(-1.5, y[nm], lab.replace(': the ', ':\nthe '), color=col, fontsize=7, va='center', ha='right')
    for w, lab in ((w1, lab1), (w2, lab2)):
        if w is None: continue
        for a, b, name, col in blocks.get(w, []):
            ax.add_patch(Rectangle((a / s, y[w] - h / 2), (min(b, xmax) - a) / s, h, color=col, lw=0))
        ax.text(-1.5, y[w], f'{c["individual"]}\'s own assembly:\n{lab}', color=TELO, fontsize=7, va='center', ha='right',
                fontweight='bold' if w == w1 else 'normal')
    # ribbons: assembly blocks placed on the reference end their label names
    import pysam
    with pysam.AlignmentFile(os.path.join(REPO, 'docs', 'evidence', 'tour', 'data', 'claims', c['slug'], 'asm.bam')) as b:
        for r in b.fetch(until_eof=True):
            w, rng = r.query_name.rsplit(':', 1); a0 = int(rng.split('-')[0]) - 1
            qs, qe = r.query_alignment_start, r.query_alignment_end
            if r.is_reverse: qs, qe = r.infer_read_length() - qe, r.infer_read_length() - qs
            ws, we = a0 + qs, a0 + qe; rs, re_ = r.reference_start, r.reference_end
            if w not in y or r.reference_name not in y or rs > xmax: continue
            col = DONOR if r.reference_name == D else RECIP
            yw, yr = y[w], y[r.reference_name]
            top, bot = (yw - h / 2, yr + h / 2) if yr < yw else (yw + h / 2, yr - h / 2)
            ax.add_patch(Polygon([(ws / s, top), (min(we, xmax) / s, top), (min(re_, xmax) / s, bot), (rs / s, bot)],
                                 color=col, alpha=0.13, lw=0))
    # molecules on the own assembly (above the claimed copy) and their two parts on the reference ends
    own = [a for a in c['own_alignments'] if a['ref'] == w1]
    for i, a in enumerate(sorted(own, key=lambda a: a['start'])):
        yy = y[w1] + h / 2 + 0.12 + i * 0.07
        ax.plot([a['start'] / s, min(a['end'], xmax) / s], [yy, yy], color=TELO, lw=0.9, solid_capstyle='butt')
    for nm in (D, Rr):
        parts = [p for p in c['parts_alignments'] if p['ref'] == nm]
        for i, p in enumerate(sorted(parts, key=lambda p: p['start'])):
            # part 1 pieces above the donor end (below its label), part 2 pieces below the recipient end
            yy = y[nm] - h / 2 - 0.12 - i * 0.07 if nm == Rr else y[nm] + h / 2 + 0.12 + i * 0.07
            ax.plot([p['start'] / s, min(p['end'], xmax) / s], [yy, yy], color=TELO, lw=0.9, solid_capstyle='butt')
    n_on = len({a['name'] for a in own})
    ax.text(xmax / s + 0.8, y[w1] + h / 2 + 0.12, f'{n_on} whole molecules', fontsize=6, va='bottom', ha='left', color=TELO)
    npd = len([p for p in c['parts_alignments'] if p['ref'] == D]); npr = len([p for p in c['parts_alignments'] if p['ref'] == Rr])
    ax.text(xmax / s + 0.8, y[D] + h / 2 + 0.12, f'part 1 of each\nmolecule ({npd})', fontsize=6, va='bottom', ha='left', color=TELO)
    ax.text(xmax / s + 0.8, y[Rr] - h / 2 - 0.12, f'part 2 of each\nmolecule ({npr})', fontsize=6, va='top', ha='left', color=TELO)
    ax.set_xlim(0, xmax / s); ax.set_ylim(-0.75, 4.75); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('distance from the telomere (kb)'); ax.set_title(title)


def main():
    style()
    J = json.load(open(os.path.join(REPO, 'docs', 'evidence', 'tour', 'data', 'tour_data.json')))
    CL = {c['slug']: c for c in J['claims']}
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.9))
    plt.subplots_adjust(left=0.25, right=0.87, top=0.93, bottom=0.135, hspace=0.45)
    for ax, (slug, title, l1, l2), letter in zip(axes, EXAMPLES, 'ab'):
        panel(ax, CL[slug], J['refs'], title, l1, l2)
        ax.text(-0.33, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', ha='left', va='bottom')
    fig.text(0.02, 0.975, 'Molecules, assembly and reference genome tell the same story', fontsize=9, va='top')
    fig.text(0.02, 0.025, 'Black lines: single DNA molecules from public sequencing runs (on the person\'s assembly: whole; on the '
             'reference: cut at their switch). Shaded bands: where each block of the person\'s assembly sits on the reference '
             '(orange: donor end, blue: recipient end). Grey blocks match a third chromosome end, not shown. All tracks start '
             'at the telomere (black).', fontsize=6, va='bottom', wrap=True)
    os.makedirs(os.path.join(REPO, 'figures', 'pitch'), exist_ok=True)
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(REPO, 'figures', 'pitch', f'fig_consolidated.{ext}'))
    return fig


if __name__ == '__main__':
    main()
