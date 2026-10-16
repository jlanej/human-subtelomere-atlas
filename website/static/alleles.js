/* Alleles page: client-side search over all structural alleles and haplogroups (window.ATLAS_INDEX). */
(function () {
  "use strict";
  var A = window.Atlas, X = window.ATLAS_INDEX;
  if (!A || !X) return;
  function init() {
    var rows = X.rows.map(function (r) {
      var o = {}; X.cols.forEach(function (c, i) { o[c] = r[i]; }); return o;
    });
    var cols = [
      { key: "id", label: "ID", render: function (r) {
        return A.el("a", { href: "arms/" + encodeURIComponent(r.arm) + ".html#" + encodeURIComponent(r.id), text: r.id });
      } },
      { key: "type", label: "Type" },
      { key: "arm", label: "Arm" },
      { key: "structure", label: "Structure / main structural allele" },
      { key: "nativity_class", label: "Nativity", render: function (r) { return (r.nativity_class || "").replace("_", " "); } },
      { key: "families", label: "Families" },
      { key: "tier", label: "Tier", render: function (r) { return A.el("span", { "class": "tier tier-" + r.tier, text: r.tier }); } },
      { key: "evidence", label: "Evidence" },
      { key: "n_ends", label: "Ends", num: true },
      { key: "n_samples", label: "Samples", num: true },
      { key: "r2_freq", label: "R2 freq.", num: true, digits: 3 },
      { key: "AFR", label: "AFR", num: true, digits: 2 },
      { key: "AMR", label: "AMR", num: true, digits: 2 },
      { key: "EAS", label: "EAS", num: true, digits: 2 },
      { key: "EUR", label: "EUR", num: true, digits: 2 },
      { key: "SAS", label: "SAS", num: true, digits: 2 }
    ];
    var t = new A.AtlasTable({ container: document.getElementById("allele-table"), columns: cols, rows: rows,
      countEl: document.getElementById("allele-count"), initialSort: { key: "n_samples", dir: -1 } });
    var arms = {}; rows.forEach(function (r) { arms[r.arm] = 1; });
    var asel = document.getElementById("q-arm");
    Object.keys(arms).sort(function (a, b) { return a.localeCompare(b, undefined, { numeric: true }); })
      .forEach(function (a) { asel.appendChild(A.el("option", { value: a, text: a })); });
    function apply() {
      var q = document.getElementById("q").value.trim().toLowerCase();
      var terms = q ? q.split(/\s+/) : [];
      var ty = document.getElementById("q-type").value, ti = document.getElementById("q-tier").value;
      var arm = asel.value, nat = document.getElementById("q-nat").value, fam = document.getElementById("q-fam").value;
      var minS = parseInt(document.getElementById("q-minsamples").value || "0", 10) || 0;
      t.setFilter(function (r) {
        if (ty && r.type !== ty) return false;
        if (ti && r.tier !== ti) return false;
        if (arm && r.arm !== arm) return false;
        if (nat && r.nativity_class !== nat) return false;
        if (fam && (" " + (r.families || "") + " ").indexOf(" " + fam + " ") < 0) return false;
        if (r.n_samples < minS) return false;
        if (terms.length) {
          var hay = (r.id + " " + r.arm + " " + r.structure + " " + r.form + " " + r.families + " " + r.nativity_class + " " + r.evidence).toLowerCase();
          for (var i = 0; i < terms.length; i++) if (hay.indexOf(terms[i]) < 0) return false;
        }
        return true;
      });
    }
    ["q", "q-type", "q-tier", "q-arm", "q-nat", "q-fam", "q-minsamples"].forEach(function (id) {
      var e = document.getElementById(id); e.addEventListener("input", apply); e.addEventListener("change", apply);
    });
    var params = new URLSearchParams(location.search);
    if (params.get("q")) { document.getElementById("q").value = params.get("q"); }
    apply();
    document.getElementById("q-download").addEventListener("click", function () {
      A.downloadText("atlas_" + X.version + "_alleles_filtered.tsv", t.toTSV());
    });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
