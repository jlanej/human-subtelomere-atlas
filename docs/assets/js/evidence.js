/* Read-evidence claims table (window.EVIDENCE_CLAIMS): search, filters, TSV export. */
(function () {
  "use strict";
  var A = window.Atlas, X = window.EVIDENCE_CLAIMS;
  if (!A || !X) return;
  var COH = { HG002_Q100: "HG002 Q100", HG002_trio: "HG003/HG004", "1KG_ONT": "1KG-ONT" };
  var DEF = { definitive_multi_dataset: "definitive (>= 2 datasets)", definitive: "definitive",
    confirmed_not_definitive: "confirmed, not definitive", not_definitive: "" };
  function init() {
    var rows = X.rows.map(function (r) { var o = {}; X.cols.forEach(function (c, i) { o[c] = r[i]; }); return o; });
    var cols = [
      { key: "claim_id", label: "Claim", render: function (r) {
        return r.page ? A.el("a", { href: "claim/" + encodeURIComponent(r.slug) + ".html", text: r.claim_id }) : r.claim_id; } },
      { key: "cohort", label: "Cohort", render: function (r) { return COH[r.cohort] || r.cohort; } },
      { key: "individual", label: "Individual" },
      { key: "haplotype", label: "Haplotype" },
      { key: "arm", label: "Arm" },
      { key: "allele_id", label: "Allele", render: function (r) {
        return A.el("a", { href: "../arms/" + encodeURIComponent(r.arm) + ".html#" + encodeURIComponent(r.allele_id), text: r.allele_id }); } },
      { key: "donor_D", label: "Donor (CHM13)", render: function (r) { return (r.donor_D || "").replace("chr", ""); } },
      { key: "v1_verdict", label: "v1 verdict" },
      { key: "definitive_verdict", label: "Definitive tier", render: function (r) { return DEF[r.definitive_verdict] || ""; } },
      { key: "n_conclusive_concordant", label: "Conclusive reads", num: true },
      { key: "n_definitive", label: "Definitive reads", num: true },
      { key: "n_definitive_datasets", label: "Datasets", num: true },
      { key: "n_contradicts", label: "Contradicting", num: true }
    ];
    var t = new A.AtlasTable({ container: document.getElementById("e-table"), columns: cols, rows: rows,
      countEl: document.getElementById("e-count"), initialSort: { key: "n_definitive", dir: -1 } });
    function apply() {
      var q = document.getElementById("e-q").value.trim().toLowerCase();
      var terms = q ? q.split(/\s+/) : [];
      var co = document.getElementById("e-cohort").value, v1 = document.getElementById("e-v1").value;
      var de = document.getElementById("e-def").value, ev = document.getElementById("e-evid").checked;
      t.setFilter(function (r) {
        if (co && r.cohort !== co) return false;
        if (v1 && r.v1_verdict !== v1) return false;
        if (de === "definitive" && String(r.definitive_verdict).indexOf("definitive") !== 0) return false;
        if (de && de !== "definitive" && r.definitive_verdict !== de) return false;
        if (ev && !r.page) return false;
        if (terms.length) {
          var hay = (r.claim_id + " " + r.individual + " " + r.arm + " " + r.allele_id + " " + r.donor_D + " " + r.recipient_R + " " +
            r.haplotype + " " + r.nativity_class).toLowerCase();
          for (var i = 0; i < terms.length; i++) if (hay.indexOf(terms[i]) < 0) return false;
        }
        return true;
      });
    }
    ["e-q", "e-cohort", "e-v1", "e-def", "e-evid"].forEach(function (id) {
      var e = document.getElementById(id); e.addEventListener("input", apply); e.addEventListener("change", apply);
    });
    document.getElementById("e-download").addEventListener("click", function () {
      var keep = t.visibleRows();
      var lines = [X.cols.join("\t")].concat(keep.map(function (r) { return X.cols.map(function (c) { return r[c] === null ? "" : r[c]; }).join("\t"); }));
      A.downloadText("read_evidence_claims.tsv", lines.join("\n") + "\n");
    });
    apply();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
