"""Digitize published arm-specific relative telomere lengths from Karimian et al. (Science 2024, doi:10.1126/science.ado0431;
preprint doi:10.1101/2023.12.21.572870), Fig. 4e: mean telomere length of each chromosome end relative to the individual's
grand mean, 147 individuals, Telomere Profiling (ONT), ends sorted by mean. Acrocentric p arms, X and Y were excluded by
the authors. Numeric per-arm values are not tabulated in the paper, so they are digitized from the figure raster.

Input : work/telomere_cis/karimian2023_biorxiv_v1.pdf (bioRxiv v1 full text,
        https://www.biorxiv.org/content/biorxiv/early/2024/01/13/2023.12.21.572870.full.pdf; same Fig. 4 as the Science paper).
        If absent, the committed data/karimian2024_fig4e.tsv is kept unchanged.
Method: page 36 raster (1262x1216 px). y axis calibrated from tick marks (+3, +2, -1, -2 kb); the per-end mean is the black
        horizontal bar (width 19-27 px) nearest each of 38 equally spaced columns. The x-axis labels (rank order) were
        transcribed from the figure; the 26th label is partly covered by a vertical line in the source raster and is read as
        18q (consistent with Fig. 4f, where 18q appears and 18p does not).
Output: data/karimian2024_fig4e.tsv (arm, rank_short_to_long, rel_mean_kb_digitized)
"""
import os, sys, io, collections, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tl_common import WORK, DAT, sha256
LABELS = ['17p', '20q', '12p', '17q', '22q', '3q', '9q', '2q', '19p', '7p', '20p', '5q', '16q', '21q', '9p', '16p', '15q', '8q',
          '6p', '11q', '11p', '10p', '1q', '14q', '4p', '18q', '19q', '10q', '6q', '8p', '1p', '5p', '7q', '13q', '2p', '4q',
          '12q', '3p']
pdf = os.path.join(WORK, 'karimian2023_biorxiv_v1.pdf'); out = os.path.join(DAT, 'karimian2024_fig4e.tsv')
if not os.path.exists(pdf):
    print('source PDF absent; keeping committed', out); raise SystemExit(0)
import pypdf
from PIL import Image
page = pypdf.PdfReader(pdf).pages[35]
im = np.array(Image.open(io.BytesIO(page.images[0].data)).convert('RGB')).astype(int)
blk = im.max(axis=2) < 70
y0, y1 = 315, 700                                   # panel e
ax_x = int(np.argmax(blk[320:700, :200].sum(axis=0)))
ticks = [y for y in range(y0, y1) if blk[y, ax_x - 6:ax_x].all()]
assert ticks == [325, 392, 592, 659], ticks
# ticks correspond to +3, +2, -1, -2 kb
tv = np.array([3, 2, -1, -2]); slope, icpt = np.polyfit(tv, np.array(ticks, float), 1)   # px = slope*kb + icpt
runs = []
for y in range(y0, y1):
    row = blk[y, ax_x:1250]; x = 0
    while x < len(row):
        if row[x]:
            s = x
            while x < len(row) and row[x]: x += 1
            if 19 <= x - s <= 27: runs.append((y, (s + x) / 2 + ax_x, x - s))
        else: x += 1
# column centres: from the full-width (25 px) bars, fit equal spacing
full = sorted({round(c) for _, c, w in runs if w >= 24})
c0, c1 = min(full), max(full); step = (c1 - c0) / (len(LABELS) - 1)
by = collections.defaultdict(list)
for y, c, w in runs:
    k = int(round((c - c0) / step))
    if 0 <= k < len(LABELS) and abs(c - (c0 + k * step)) < 8: by[k].append((w, y))
rows = []
for k, a in enumerate(LABELS):
    v = by.get(k, [])
    wmax = max(w for w, _ in v); ys = [y for w, y in v if w == wmax]
    rows.append(dict(arm=a, rank_short_to_long=k + 1, rel_mean_kb_digitized=round((np.mean(ys) - icpt) / slope, 3), bar_px_width=wmax))
D = pd.DataFrame(rows)
assert D.rel_mean_kb_digitized.is_monotonic_increasing or (D.rel_mean_kb_digitized.diff().dropna() > -0.03).all()
D.to_csv(out, sep='\t', index=False)
print('digitized', len(D), 'ends; px/kb', round(-slope, 2), 'zero px', round(icpt, 1), '| source sha256', sha256(pdf)[:16])
