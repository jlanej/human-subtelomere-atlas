/* Render check for docs/ without a browser: loads pages from file:// in jsdom with their scripts, exercises the
   interactive parts and writes the generated SVG mosaics to an output folder for visual inspection.
   Usage (from the repository root, jsdom installed e.g. with `npm install --prefix work/site_check jsdom@24`):
     NODE_PATH=work/site_check/node_modules node website/render_check.js work/site_check
*/
"use strict";
const fs = require("fs");
const path = require("path");
const { JSDOM, VirtualConsole } = require("jsdom");

const root = path.resolve(__dirname, "..");
const docs = path.join(root, "docs");
const out = path.resolve(process.argv[2] || path.join(root, "work", "site_check"));
fs.mkdirSync(out, { recursive: true });

async function load(rel, hash) {
  const errors = [];
  const vc = new VirtualConsole();
  vc.on("jsdomError", (e) => errors.push(String(e && e.message || e)));
  vc.on("error", (e) => errors.push(String(e)));
  const file = path.join(docs, rel);
  const url = require("url").pathToFileURL(file).href + (hash || "");
  const dom = await JSDOM.fromFile(file, { url, runScripts: "dangerously", resources: "usable", virtualConsole: vc,
    pretendToBeVisual: true,
    // jsdom does not implement layout; scrollIntoView exists in every browser
    beforeParse(win) { win.Element.prototype.scrollIntoView = function () {}; } });
  await new Promise((res) => { if (dom.window.document.readyState === "complete") res(); else dom.window.addEventListener("load", res); });
  await new Promise((r) => setTimeout(r, 50));
  return { dom, w: dom.window, d: dom.window.document, errors };
}

function fire(w, el, type) { el.dispatchEvent(new w.Event(type, { bubbles: true })); }

(async () => {
  const report = { pages: {} };
  // arm pages
  const arms = fs.readdirSync(path.join(docs, "arms")).filter((f) => f.endsWith(".html") && f !== "index.html");
  for (const f of arms) {
    const arm = f.slice(0, -5);
    const { w, d, errors } = await load("arms/" + f);
    const rows = d.querySelectorAll("#allele-mosaic svg g.row").length;
    const nAll = w.ATLAS_ARM.alleles.length;
    const members = d.querySelectorAll("#member-table tbody tr").length;
    // keyboard selection of the first allele filters the member table
    const g = d.querySelector("#allele-mosaic svg g.row");
    let filtered = null, expected = null;
    if (g) {
      const id = g.getAttribute("data-id");
      g.dispatchEvent(new w.KeyboardEvent("keydown", { key: "Enter", bubbles: true }));
      filtered = d.querySelectorAll("#member-table tbody tr").length;
      expected = w.ATLAS_ARM.members.filter((m) => m[5] === id).length;
    }
    // tier filter: hide tier C
    const cb = d.getElementById("f-tier-C");
    cb.checked = false; fire(w, cb, "change");
    const rowsNoC = d.querySelectorAll("#allele-mosaic svg g.row").length;
    const expNoC = w.ATLAS_ARM.alleles.filter((a) => a.tier !== "C").length;
    cb.checked = true; fire(w, cb, "change");
    // nativity filter
    const nsel = d.getElementById("f-nat");
    let natRows = null, natExp = null;
    if (nsel) {
      nsel.value = "foreign"; fire(w, nsel, "change");
      natRows = d.querySelectorAll("#allele-mosaic svg g.row").length;
      natExp = w.ATLAS_ARM.alleles.filter((a) => a.nc === "foreign").length;
      nsel.value = ""; fire(w, nsel, "change");
    }
    // tooltip content on keyboard focus of the first row
    const g0 = d.querySelector("#allele-mosaic svg g.row");
    let tipOk = null;
    if (g0) {
      g0.dispatchEvent(new w.FocusEvent("focus"));
      const t = d.getElementById("atlas-tooltip");
      const a0 = w.ATLAS_ARM.alleles.find((a) => a.id === g0.getAttribute("data-id"));
      tipOk = !!t && t.innerHTML.indexOf("Nativity class") >= 0 && t.innerHTML.indexOf(a0.id) >= 0 &&
        a0.blk.every((b) => !w.ATLAS_PALETTE.families[b[3]] || t.innerHTML.indexOf(b[3] + " (" + w.ATLAS_PALETTE.families[b[3]].join(", ") + ")") >= 0);
      g0.dispatchEvent(new w.FocusEvent("blur"));
    }
    // member mosaic toggle
    d.getElementById("toggle-members-mosaic").click();
    const mrects = d.querySelectorAll("#member-mosaic svg rect").length;
    const svg = d.querySelector("#allele-mosaic svg");
    if (["20q", "3q", "6p", "8p", "16p", "11p"].includes(arm) && svg) {
      fs.writeFileSync(path.join(out, `arm_${arm}_alleles.svg`), svg.outerHTML.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" '));
      const ms = d.querySelector("#member-mosaic svg");
      if (ms) fs.writeFileSync(path.join(out, `arm_${arm}_members.svg`), ms.outerHTML.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" '));
    }
    report.pages["arms/" + f] = { allele_rows: rows, alleles: nAll, member_rows: members, members: w.ATLAS_ARM.members.length,
      first_allele_filtered_rows: filtered, first_allele_expected: expected, rows_without_tierC: rowsNoC,
      expected_without_tierC: expNoC, member_mosaic_rects: mrects, foreign_rows: natRows, foreign_expected: natExp,
      tooltip_ok: tipOk, errors };
  }
  // hash deep link
  {
    const { d, errors } = await load("arms/20q.html", "#SA-20q-01");
    report.deep_link = { selected: !!d.querySelector('#allele-mosaic g.row.selected[data-id="SA-20q-01"]'),
      member_rows: d.querySelectorAll("#member-table tbody tr").length, errors };
  }
  // alleles page search
  {
    const { w, d, errors } = await load("alleles.html");
    const all = d.querySelectorAll("#allele-table tbody tr").length;
    const q = d.getElementById("q"); q.value = "20q 3q"; fire(w, q, "input");
    const hits = Array.from(d.querySelectorAll("#allele-table tbody tr")).map((tr) => tr.cells[0].textContent);
    const ty = d.getElementById("q-type"); ty.value = "haplogroup"; q.value = ""; fire(w, q, "input"); fire(w, ty, "change");
    const hgRows = d.querySelectorAll("#allele-table tbody tr").length;
    ty.value = ""; fire(w, ty, "change");
    const fam = d.getElementById("q-fam"); const nat = d.getElementById("q-nat");
    const famRows = {};
    for (const f of Object.keys(w.ATLAS_INDEX.families)) {
      fam.value = f; fire(w, fam, "change");
      famRows[f] = [d.querySelectorAll("#allele-table tbody tr").length,
        w.ATLAS_INDEX.rows.filter((r) => (" " + r[6] + " ").indexOf(" " + f + " ") >= 0).length];
    }
    fam.value = ""; nat.value = "within_family"; fire(w, nat, "change");
    const wfRows = [d.querySelectorAll("#allele-table tbody tr").length, w.ATLAS_INDEX.rows.filter((r) => r[5] === "within_family").length];
    nat.value = ""; fire(w, nat, "change");
    q.value = "20q PF04"; fire(w, q, "input");
    const pf04 = Array.from(d.querySelectorAll("#allele-table tbody tr")).map((tr) => tr.cells[0].textContent).slice(0, 5);
    report.pages["alleles.html"] = { rows: all, index_rows: w.ATLAS_INDEX.rows.length, search_20q_3q: hits.slice(0, 8),
      search_hits: hits.length, haplogroup_rows: hgRows, family_filter_rows: famRows, within_family_rows: wfRows,
      search_20q_PF04: pf04, errors };
    if (Object.values(famRows).some((x) => x[0] !== x[1]) || wfRows[0] !== wfRows[1]) { report.alleles_fail = true; }
  }
  // other pages: script errors only
  const others = [];
  (function walk(dir) {
    for (const f of fs.readdirSync(dir)) {
      const p = path.join(dir, f);
      if (fs.statSync(p).isDirectory()) walk(p);
      else if (f.endsWith(".html")) others.push(path.relative(docs, p));
    }
  })(docs);
  for (const rel of others) {
    if (report.pages[rel]) continue;
    const { errors } = await load(rel);
    report.pages[rel] = { errors };
  }
  let nerr = 0, nfail = 0;
  for (const [k, v] of Object.entries(report.pages)) {
    nerr += v.errors.length;
    if (k.startsWith("arms/") && k !== "arms/index.html" && (v.allele_rows !== v.alleles || v.member_rows !== v.members ||
        v.first_allele_filtered_rows !== v.first_allele_expected || v.rows_without_tierC !== v.expected_without_tierC ||
        !(v.member_mosaic_rects > 0) || v.foreign_rows !== v.foreign_expected || v.tooltip_ok !== true)) {
      nfail++; console.log("FAIL", k, JSON.stringify(v)); }
  }
  if (!report.deep_link.selected) { nfail++; console.log("FAIL deep link"); }
  if (report.alleles_fail) { nfail++; console.log("FAIL alleles filters"); }
  report.summary = { pages: Object.keys(report.pages).length, script_errors: nerr, failed_checks: nfail };
  fs.writeFileSync(path.join(out, "render_report.json"), JSON.stringify(report, null, 1));
  console.log(JSON.stringify(report.summary), JSON.stringify(report.pages["alleles.html"]), JSON.stringify(report.deep_link));
  process.exit(nerr || nfail ? 1 : 0);
})();
