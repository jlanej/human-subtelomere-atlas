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
html5lib (QC). The site has no runtime dependencies: plain JavaScript in `website/static/` (no third-party libraries,
no fonts, trackers or CDNs).

## Files

| file | role |
|---|---|
| `pipeline/atlas_palette.py` | family-aware CHM13 end palette with CIEDE2000 / deuteranopia QC |
| `pipeline/figures_atlas.py` | F1 genome overview (chords + nativity classes), F2 allele mosaics, F3 superpopulation frequencies, F4 validation, palette key; captions computed from the tables |
| `website/build_site.py` | builds every page; renders `QC.md`, `SCHEMA.md`, `CATALOGUE.md` and the verdict sections of the analysis READMEs with a small built-in Markdown renderer; reads the stability tables (QC.md section 9), `tests/validation/VALIDATION.md` and `data/validation/` (Tan et al. cross-check) at build time and shows an "in progress" notice when they are absent; builds the analysis-form vs revision-2 concordance table from `ends.tsv` |
| `website/static/` | `site.css`, `common.js` (sortable/filterable tables, TSV export, tooltips), `arm.js` (interactive allele mosaic, member mosaic and member table), `alleles.js` (allele and haplogroup search) |
| `website/check_site.py` | strict HTML5 parse, internal links and anchors, repository links, no external runtime resources, JSON/JS payloads, local-path scan, 15 MB size budget |
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
