/* Arm page: interactive allele mosaic, member-end mosaic and member table.
   Data: window.ATLAS_ARM (data/arms/<arm>.js), window.ATLAS_PALETTE (data/palette.js). */
(function () {
  "use strict";
  var A = window.Atlas, D = window.ATLAS_ARM, P = window.ATLAS_PALETTE;
  if (!A || !D || !P) return;
  var SVGNS = "http://www.w3.org/2000/svg";
  var POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"];
  var KB = 100000;
  // geometry (viewBox units)
  var W = 1000, X0 = 168, X1 = 760, TEL0 = 150, ROWH = 15, TOP = 26;
  function xs(b) { return X0 + (Math.max(0, Math.min(KB, b)) / KB) * (X1 - X0); }
  function s(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    for (var k in attrs) if (Object.prototype.hasOwnProperty.call(attrs, k)) e.setAttribute(k, attrs[k]);
    return e;
  }
  function srcName(i) { return P.arm[i]; }
  function colour(i) { return P.hex[i]; }

  var selected = null;
  var memberTable = null;

  function alleleTip(a) {
    var esc = A.esc;
    var sw = a.sw.length ? a.sw.map(function (v, i) {
      return (v / 1000).toFixed(1) + " kb" + (a.swmin[i] !== a.swmax[i] ?
        " (" + (a.swmin[i] / 1000).toFixed(1) + "\u2013" + (a.swmax[i] / 1000).toFixed(1) + ")" : "");
    }).join("; ") : "none (single block)";
    var fq = "<table><tr>" + POPS.map(function (p) { return "<th>" + p + "</th>"; }).join("") + "<th>all</th></tr><tr>" +
      POPS.map(function (p) { return "<td>" + A.fmt(a.fq[p], 2) + "</td>"; }).join("") + "<td>" + A.fmt(a.f, 3) +
      "</td></tr></table>";
    var c = a.evc;
    var NC = { native: "native: only home (and none) blocks",
      within_family: "within family: non-home blocks all from a paralog family that contains this arm",
      foreign: "foreign: at least one block from another family or arm" };
    var blocks = a.blk.map(function (b) {
      var lab = b[3], fam = P.families[lab];
      var name = lab === "home" ? "home" : lab === "none" ? "none (no CHM13 end within 0.80 identity)" :
        fam ? lab + " (" + fam.join(", ") + ")" : "CHM13 " + lab;
      return "<tr><td>" + (b[0] / 1000).toFixed(1) + "\u2013" + (b[1] / 1000).toFixed(1) + " kb</td><td>" + esc(name) +
        "</td><td>" + (b[4] ? "best end " + esc(b[4]) : "") + "</td></tr>";
    }).join("");
    return "<strong>" + esc(a.id) + "</strong> <span class='tier tier-" + a.tier + "'>" + a.tier + "</span><br>" +
      "Structure (telomere \u2192 centromere): <code>" + esc(a.st) + "</code>" + (a.form ? " \u00b7 form " + esc(a.form) : "") +
      "<br>Nativity class: " + esc(NC[a.nc] || a.nc || "") +
      "<table>" + blocks + "</table>" +
      "<br>Median core switch(es) from the telomere boundary: " + sw +
      "<br>" + a.ne + " ends, " + a.ns + " samples; Release 2 frequency " + A.fmt(a.f, 3) +
      "<br>Frequency by superpopulation:" + fq +
      "Evidence: " + (a.ev ? esc(a.ev.split(",").join(", ")) : "none (tier C)") +
      "<br><span class='muted'>Verkko concordant " + c.vkc + ", discordant " + c.vkd + "; read layout " + c.lay +
      "; HG002 reads " + c.rd + "; pedigree intact " + c.ped + "; T2T members " + c.t2t + "</span>" +
      "<br><span class='muted'>Representative end: " + esc(a.rep) + "</span>";
  }

  function visibleAlleles() {
    var tiers = {};
    ["A", "B", "C"].forEach(function (t) { var cb = document.getElementById("f-tier-" + t); tiers[t] = !cb || cb.checked; });
    var minS = parseInt((document.getElementById("f-minsamples") || {}).value || "1", 10) || 1;
    var nat = (document.getElementById("f-nat") || {}).value || "";
    var order = (document.getElementById("f-order") || {}).value || "members";
    var al = D.alleles.filter(function (a) { return tiers[a.tier] && a.ns >= minS && (!nat || a.nc === nat); });
    if (order === "switch") al.sort(function (a, b) { return (a.sw[0] || 1e9) - (b.sw[0] || 1e9) || b.ne - a.ne; });
    else if (order === "id") al.sort(function (a, b) { return a.id.localeCompare(b.id, undefined, { numeric: true }); });
    else al.sort(function (a, b) { return b.ne - a.ne; });
    return al;
  }

  function drawAxis(svg, y, label) {
    [0, 25, 50, 75, 100].forEach(function (k) {
      var x = xs(k * 1000);
      svg.appendChild(s("line", { x1: x, x2: x, y1: y - 4, y2: y, stroke: "#555", "stroke-width": 1 }));
      var t = s("text", { x: x, y: y - 7, "text-anchor": "middle", "font-size": 10, fill: "#333" });
      t.textContent = k + (k === 100 ? " kb" : ""); svg.appendChild(t);
    });
    var t2 = s("text", { x: TEL0, y: y - 7, "font-size": 10, fill: "#333", "text-anchor": "middle" });
    t2.textContent = "TEL"; svg.appendChild(t2);
    if (label) {
      var t3 = s("text", { x: 790, y: y - 7, "font-size": 10, fill: "#333" });
      t3.textContent = label; svg.appendChild(t3);
    }
  }

  function renderAlleles() {
    var host = document.getElementById("allele-mosaic");
    var al = visibleAlleles();
    var H = TOP + al.length * ROWH + 8;
    var svg = s("svg", { viewBox: "0 0 " + W + " " + H, "class": "mosaic", role: "group",
      "aria-label": "Structural alleles of " + D.arm + ": one row per allele, blocks coloured by the most similar CHM13 end" });
    drawAxis(svg, TOP - 4, "Release 2 frequency \u00b7 samples");
    var fmax = Math.max.apply(null, al.map(function (a) { return a.f || 0; }).concat([0.05]));
    al.forEach(function (a, i) {
      var y = TOP + i * ROWH;
      var g = s("g", { "class": "row" + (selected === a.id ? " selected" : ""), tabindex: 0, role: "button",
        "aria-label": a.id + ", tier " + a.tier + ", " + (a.nc || "").replace("_", " ") + ", structure " + a.st + ", " + a.ns + " samples, frequency " +
          A.fmt(a.f, 3) + ". Press Enter to show members.", "data-id": a.id });
      g.appendChild(s("rect", { "class": "hit", x: 0, y: y, width: W, height: ROWH, fill: "transparent" }));
      var lab = s("text", { x: 4, y: y + 11, "font-size": 10.5, fill: "#1c1c1c" }); lab.textContent = a.id; g.appendChild(lab);
      var tb = s("rect", { x: 118, y: y + 2, width: 16, height: 11, fill: a.tier === "A" ? "#2b2b2b" : a.tier === "B" ? "#6b6b6b" : "#dcdcdc" });
      g.appendChild(tb);
      var tt = s("text", { x: 126, y: y + 11, "font-size": 9, "text-anchor": "middle", fill: a.tier === "C" ? "#1c1c1c" : "#fff", "font-weight": 700 });
      tt.textContent = a.tier; g.appendChild(tt);
      g.appendChild(s("rect", { x: TEL0 - 6, y: y + 2, width: 12, height: ROWH - 4, fill: "#3b3b3b" }));
      a.blk.forEach(function (b) {
        var attrs = { x: xs(b[0]), y: y + 2, width: Math.max(0.6, xs(b[1]) - xs(b[0])), height: ROWH - 4, fill: colour(b[2]) };
        if (P.families[b[3]]) { attrs.stroke = "#1c1c1c"; attrs["stroke-width"] = 1; attrs["stroke-dasharray"] = "3 2"; }
        g.appendChild(s("rect", attrs));
      });
      a.sw.forEach(function (v, k) {
        if (a.swmax[k] > a.swmin[k]) g.appendChild(s("rect", { x: xs(a.swmin[k]), y: y + ROWH - 3, width: Math.max(1, xs(a.swmax[k]) - xs(a.swmin[k])), height: 2, fill: "#1c1c1c", opacity: 0.45 }));
        g.appendChild(s("path", { d: "M" + (xs(v) - 3) + " " + (y) + " h6 l-3 4 z", fill: "#1c1c1c" }));
      });
      var fw = 120 * (a.f || 0) / fmax;
      g.appendChild(s("rect", { x: 790, y: y + 3, width: Math.max(0.8, fw), height: ROWH - 6, fill: "#4d4d4d" }));
      var ft = s("text", { x: 790 + fw + 4, y: y + 11, "font-size": 10, fill: "#1c1c1c" });
      ft.textContent = A.fmt(a.f, 3) + " \u00b7 " + a.ns; g.appendChild(ft);
      var title = s("title", {}); title.textContent = a.id + " " + a.st; g.appendChild(title);
      function show(ev) {
        var r = g.getBoundingClientRect();
        var x = ev && ev.clientX !== undefined ? ev.clientX : r.left + r.width * 0.5;
        var yy = ev && ev.clientY !== undefined ? ev.clientY : r.bottom;
        A.showTip(alleleTip(a), x, yy);
      }
      g.addEventListener("mousemove", show);
      g.addEventListener("focus", function () { show(null); });
      g.addEventListener("mouseleave", A.hideTip);
      g.addEventListener("blur", A.hideTip);
      function pick() { selectAllele(selected === a.id ? null : a.id, true); }
      g.addEventListener("click", pick);
      g.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(); } });
      svg.appendChild(g);
    });
    host.innerHTML = "";
    host.appendChild(svg);
    var cnt = document.getElementById("allele-count");
    if (cnt) cnt.textContent = al.length + " of " + D.alleles.length + " structural alleles shown";
  }

  function renderLegend() {
    var used = {};
    D.alleles.forEach(function (a) { a.blk.forEach(function (b) { used[b[2]] = 1; }); });
    D.members.forEach(function (m) { (m[m.length - 1] || []).forEach(function (b) { used[b[2]] = 1; }); });
    var host = document.getElementById("mosaic-legend");
    if (!host) return;
    host.innerHTML = "";
    Object.keys(used).map(Number).sort(function (a, b) { return a - b; }).forEach(function (i) {
      host.appendChild(A.el("span", {}, [A.el("span", { "class": "sw", style: "background:" + colour(i) }),
        srcName(i) === "none" ? "none (no CHM13 end within 0.80)" :
          "CHM13 " + srcName(i) + (P.family_of[srcName(i)] ? " (" + P.family_of[srcName(i)] + ")" : "")]));
    });
    host.appendChild(A.el("span", { "class": "muted" }, ["Dashed outline = block labelled by paralog family: " +
      Object.keys(P.families).map(function (f) { return f + " = " + P.families[f].join("/"); }).join(", ")]));
  }

  function renderMembersMosaic() {
    var host = document.getElementById("member-mosaic");
    var al = visibleAlleles();
    var idx = {}; D.mcols.forEach(function (c, i) { idx[c] = i; });
    var groups = al.map(function (a) {
      return [a, D.members.filter(function (m) {
        return m[idx.structural_allele_id] === a.id && m[idx.set] !== "pedigree_aux";
      })];
    }).filter(function (g) { return g[1].length; });
    var n = groups.reduce(function (t, g) { return t + g[1].length; }, 0);
    var rh = n > 900 ? 1 : 2, gap = 4;
    var H = TOP + n * rh + groups.length * gap + 8;
    var svg = s("svg", { viewBox: "0 0 " + W + " " + H, "class": "mosaic", role: "img",
      "aria-label": "Mosaic of " + n + " member ends of " + D.arm + " grouped by structural allele" });
    drawAxis(svg, TOP - 4, "allele");
    var y = TOP;
    groups.forEach(function (g) {
      var y0 = y;
      g[1].sort(function (a, b) {
        var ba = a[a.length - 1] || [], bb = b[b.length - 1] || [];
        return (ba.length > 1 ? ba[0][1] : 1e9) - (bb.length > 1 ? bb[0][1] : 1e9);
      });
      g[1].forEach(function (m) {
        svg.appendChild(s("rect", { x: TEL0 - 6, y: y, width: 12, height: rh, fill: "#3b3b3b" }));
        (m[m.length - 1] || []).forEach(function (b) {
          svg.appendChild(s("rect", { x: xs(b[0]), y: y, width: Math.max(0.5, xs(b[1]) - xs(b[0])), height: rh, fill: colour(b[2]) }));
        });
        y += rh;
      });
      if (y - y0 >= 9) {
        var t = s("text", { x: 790, y: (y0 + y) / 2 + 4, "font-size": 10, fill: "#1c1c1c" });
        t.textContent = g[0].id + " (" + g[1].length + ")"; svg.appendChild(t);
      }
      y += gap;
    });
    host.innerHTML = "";
    host.appendChild(svg);
  }

  function buildMemberTable() {
    var idx = {}; D.mcols.forEach(function (c, i) { idx[c] = i; });
    var rows = D.members.map(function (m) {
      var r = {}; D.mcols.forEach(function (c, i) { r[c] = m[i]; }); return r;
    });
    var cols = [
      { key: "end_id", label: "End ID" },
      { key: "sample", label: "Sample" },
      { key: "superpopulation", label: "Superpop." },
      { key: "set", label: "Set" },
      { key: "assembly_accession", label: "Assembly", render: function (r) {
        var v = r.assembly_accession;
        if (!v) return "";
        if (/^GC[AF]_\d+\.\d+$/.test(v)) return A.el("a", { href: "https://www.ncbi.nlm.nih.gov/datasets/genome/" + v + "/", text: v });
        return v;
      } },
      { key: "structural_allele_id", label: "Structural allele" },
      { key: "haplogroup_id", label: "Haplogroup" },
      { key: "nativity_class", label: "Nativity", title: "Class of the end's own core structure" },
      { key: "S20_state", label: "S20 state", title: "Majority label of the first 20 kb of subtelomere" },
      { key: "S20_best_end", label: "S20 best end", title: "Best single CHM13 end over the first 20 kb" },
      { key: "form", label: "Form" },
      { key: "sa_assignment", label: "Assignment" },
      { key: "flags", label: "Flags" }
    ];
    memberTable = new A.AtlasTable({ container: document.getElementById("member-table"), columns: cols, rows: rows,
      countEl: document.getElementById("member-count"), initialSort: { key: "structural_allele_id", dir: 1 } });
    var sel = document.getElementById("m-allele");
    D.alleles.slice().sort(function (a, b) { return a.id.localeCompare(b.id, undefined, { numeric: true }); })
      .forEach(function (a) { sel.appendChild(A.el("option", { value: a.id, text: a.id + " (" + a.ne + ")" })); });
    var hsel = document.getElementById("m-hg");
    D.haplogroups.forEach(function (h) { hsel.appendChild(A.el("option", { value: h.id, text: h.id + " (" + h.ne + ")" })); });
    ["m-search", "m-allele", "m-hg", "m-nat", "m-pop", "m-set"].forEach(function (id) {
      var e = document.getElementById(id); e.addEventListener("input", applyMemberFilter); e.addEventListener("change", applyMemberFilter);
    });
    document.getElementById("m-download").addEventListener("click", function () {
      A.downloadText("atlas_" + D.version + "_" + D.arm + "_members.tsv", memberTable.toTSV());
    });
  }

  function applyMemberFilter() {
    var q = document.getElementById("m-search").value.trim().toLowerCase();
    var al = document.getElementById("m-allele").value, hg = document.getElementById("m-hg").value;
    var pop = document.getElementById("m-pop").value, set = document.getElementById("m-set").value;
    var nat = document.getElementById("m-nat").value;
    memberTable.setFilter(function (r) {
      if (al && r.structural_allele_id !== al) return false;
      if (hg && r.haplogroup_id !== hg) return false;
      if (pop && (r.superpopulation || "NA") !== pop) return false;
      if (set && r.set !== set) return false;
      if (nat && r.nativity_class !== nat) return false;
      if (q) {
        var hay = (r.end_id + " " + r.sample + " " + (r.assembly_accession || "") + " " + (r.flags || "")).toLowerCase();
        if (hay.indexOf(q) < 0) return false;
      }
      return true;
    });
  }

  function selectAllele(id, scroll) {
    selected = id;
    document.getElementById("m-allele").value = id || "";
    applyMemberFilter();
    renderAlleles();
    if (id && history.replaceState) {
      try { history.replaceState(null, "", "#" + id); } catch (e) { /* file:// pages have an opaque origin */ }
    }
    var st = document.getElementById("selection-status");
    if (st) st.textContent = id ? "Member table filtered to " + id + "." : "";
    if (scroll && id) document.getElementById("members").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function init() {
    ["f-tier-A", "f-tier-B", "f-tier-C", "f-nat", "f-minsamples", "f-order"].forEach(function (id) {
      var e = document.getElementById(id);
      if (e) { e.addEventListener("change", function () { renderAlleles(); if (memberShown) renderMembersMosaic(); });
        e.addEventListener("input", function () { renderAlleles(); if (memberShown) renderMembersMosaic(); }); }
    });
    var memberShown = false;
    var btn = document.getElementById("toggle-members-mosaic");
    btn.addEventListener("click", function () {
      memberShown = !memberShown;
      btn.setAttribute("aria-expanded", memberShown ? "true" : "false");
      btn.textContent = memberShown ? "Hide member-end mosaic" : "Show member-end mosaic";
      var host = document.getElementById("member-mosaic");
      if (memberShown) { host.hidden = false; renderMembersMosaic(); } else { host.hidden = true; host.innerHTML = ""; }
    });
    buildMemberTable();
    renderLegend();
    var h = decodeURIComponent((location.hash || "").slice(1));
    if (h && D.alleles.some(function (a) { return a.id === h; })) selectAllele(h, false);
    else if (h && D.haplogroups.some(function (x) { return x.id === h; })) {
      document.getElementById("m-hg").value = h; applyMemberFilter(); renderAlleles();
      var e = document.getElementById(h); if (e) e.scrollIntoView();
    } else renderAlleles();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
