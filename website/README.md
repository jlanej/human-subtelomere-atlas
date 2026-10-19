# Website and atlas figures

The GitHub Pages site in `docs/` and the atlas figures in `figures/atlas/` are generated from the release tables in
`data/catalogue/`, the analyses in `analyses/` and the literature files in `literature/`. Nothing in `docs/` is edited by
hand; serve the `docs/` folder of the main branch (Settings > Pages). The pages use relative links only, so they also
work when opened from a local checkout (`docs/index.html`).

```bash
python pipeline/atlas_palette.py   # figures/atlas/PALETTE.tsv + palette_qc.tsv (only if PALETTE.tsv is missing; --force to redo)
python pipeline/figures_atlas.py   # figures/atlas/*.png|pdf|svg, CAPTIONS.md, figure_stats.json
python website/build_site.py       # docs/ (static HTML/CSS/JS + data/arms/<arm>.js, data/alleles_index.json)
python website/check_site.py       # QC report -> work/site_check/report.json (needs html5lib)
# optional browser-free render test (needs Node.js and jsdom):
npm install --prefix work/site_check jsdom@24
NODE_PATH=work/site_check/node_modules node website/render_check.js work/site_check
```

Requirements: Python 3.9+ with pandas, numpy and matplotlib (figures); Pillow (optional, web-sized image copies);
html5lib (QC). The site has no external runtime resources: plain JavaScript in `website/static/` (no fonts, trackers or
CDNs). The one third-party library is igv.js 3.8.10 (MIT licence, `website/static/vendor/igv.min.js` and `IGV_LICENSE.txt`,
from the npm package `igv`), served from the site itself and used only by the guided tour.

## Guided tour (igv.js)

`docs/evidence/tour.html` walks through the read evidence in an embedded genome browser. Its data are small BAM, BED and
bgzipped FASTA files in `docs/evidence/tour/data/` (one combined reference per claim: the person's own end windows and the
CHM13 v2.0 donor and recipient ends, all in tip frame), built from the v0.1.0 release assets:

```bash
# release assets subtel_read_evidence_v1_{HG002_Q100,HG002_trio,1KG_ONT,tables}.tar.gz unpacked to work/tour/release/
export TOUR_CACHE=work/tour                 # downloads (CHM13 end ranges) and temporary files
python website/tour_build.py               # molecules, parts, own-assembly and assembly-to-CHM13 alignments
python website/tour_models.py              # formal claim models (criteria v1) in tip frame
python website/tour_build.py --sites-only  # per-molecule agreement with the claim model; strips @PG lines from BAMs
python website/build_site.py               # writes evidence/tour.html and evidence/tour/tour.json
python website/tour_check.py               # static check: files, contigs, loci, BAM integrity, file sizes
python report/pitch_consolidated.py        # static four-track figure figures/pitch/fig_consolidated.png (tour data)
```

Molecules are put in tip frame exactly as in the release (`readevidence.orient` with the criteria-v1 telomere motifs,
including the ONT R9 basecalling artefacts); the builder asserts that each molecule reproduces its released telomere
length. A browser render cannot be tested without a display; check the live page after each push.

## Files

| file | role |
|---|---|
| `pipeline/atlas_palette.py` | family-aware CHM13 end palette with CIEDE2000 / deuteranopia QC |
| `pipeline/figures_atlas.py` | F1 genome overview (chords + nativity classes), F2 allele mosaics, F3 superpopulation frequencies, F4 validation, palette key; captions computed from the tables |
| `website/build_site.py` | builds every page; renders `QC.md`, `SCHEMA.md`, `CATALOGUE.md` and the verdict sections of the analysis READMEs with a small built-in Markdown renderer; reads the stability tables (QC.md section 9), `tests/validation/VALIDATION.md` and `data/validation/` (Tan et al. cross-check) at build time and shows an "in progress" notice when they are absent; builds the analysis-form vs revision-2 concordance table from `ends.tsv` |
| `website/tour_build.py`, `tour_models.py`, `tour_pages.py`, `tour_check.py`, `static/tour.js` | guided tour: data, claim models, page and chapter list, static check, carousel |
| `website/static/` | `site.css`, `common.js` (sortable/filterable tables, TSV export, tooltips), `arm.js` (interactive allele mosaic, member mosaic and member table), `alleles.js` (allele and haplogroup search) |
| `website/check_site.py` | strict HTML5 parse, internal links and anchors, repository links, no external runtime resources, JSON/JS payloads, local-path scan, 150 MB size budget (the tour data are most of it) |
| `website/render_check.js` | loads every page with its scripts in jsdom, checks that mosaics and tables render and that filters, keyboard selection and deep links work; writes the generated SVG mosaics for inspection |

## Colour palette

`pipeline/atlas_palette.py` writes `figures/atlas/PALETTE.tsv` (one colour per CHM13 v2.0 end) and
`figures/atlas/palette_qc.tsv` (CIEDE2000 difference of every pair, in normal vision and in simulated deuteranopia after
Machado et al. 2009, with the pair's co-occurrence weight in the mosaics). Members of a paralog family
(`data/catalogue/paralog_families.tsv`) share one hue and differ in lightness; every other end has its own colour. The
assignment was chosen once by a seeded optimisation that maximises the colour difference between CHM13 ends that are
adjacent or co-occur in the catalogue mosaics; it is reused as long as `PALETTE.tsv` exists (`--force` re-optimises).
Home blocks are drawn in the arm's own colour (pooled groups: the member end); other blocks in the colour of their best
single CHM13 end; family-labelled blocks get a dashed outline on the website.
