/* human-subtelomere-atlas: shared helpers (vanilla JS, no dependencies, works from file://). */
(function () {
  "use strict";

  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    if (attrs) {
      for (var k in attrs) {
        if (!Object.prototype.hasOwnProperty.call(attrs, k)) continue;
        var v = attrs[k];
        if (v === null || v === undefined || v === false) continue;
        if (k === "text") e.textContent = v;
        else if (k === "html") e.innerHTML = v;
        else if (k.slice(0, 2) === "on") e.addEventListener(k.slice(2), v);
        else e.setAttribute(k, v === true ? "" : v);
      }
    }
    (children || []).forEach(function (c) {
      if (c === null || c === undefined) return;
      e.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
    });
    return e;
  }

  function esc(s) {
    return String(s === null || s === undefined ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function fmt(v, digits) {
    if (v === null || v === undefined || v === "" || (typeof v === "number" && isNaN(v))) return "";
    if (typeof v === "number" && digits !== undefined) return v.toFixed(digits);
    return String(v);
  }

  /* Sortable, filterable table.
     opts: {container, columns: [{key, label, num, digits, render(row) -> Node|string, title}], rows: [obj],
            filter: fn(row) -> bool, caption, countEl, initialSort: {key, dir}} */
  function AtlasTable(opts) {
    this.o = opts;
    this.sortKey = opts.initialSort ? opts.initialSort.key : null;
    this.sortDir = opts.initialSort ? opts.initialSort.dir : 1;
    this.filter = opts.filter || null;
    this.build();
    this.render();
  }
  AtlasTable.prototype.build = function () {
    var self = this, o = this.o;
    var table = el("table", {});
    if (o.caption) table.appendChild(el("caption", { "class": "muted", text: o.caption }));
    var tr = el("tr");
    this.ths = [];
    o.columns.forEach(function (c) {
      var th = el("th", { scope: "col", "class": c.num ? "num" : null, title: c.title || null });
      var b = el("button", { type: "button", "class": "sort", text: c.label,
        "aria-label": "Sort by " + c.label,
        onclick: function () { self.sortBy(c.key); } });
      th.appendChild(b);
      self.ths.push([c.key, th]);
      tr.appendChild(th);
    });
    table.appendChild(el("thead", {}, [tr]));
    this.tbody = el("tbody");
    table.appendChild(this.tbody);
    o.container.innerHTML = "";
    o.container.appendChild(table);
  };
  AtlasTable.prototype.sortBy = function (key) {
    if (this.sortKey === key) this.sortDir = -this.sortDir; else { this.sortKey = key; this.sortDir = 1; }
    this.render();
  };
  AtlasTable.prototype.setFilter = function (fn) { this.filter = fn; this.render(); };
  AtlasTable.prototype.visibleRows = function () {
    var f = this.filter, rows = this.o.rows.filter(function (r) { return !f || f(r); });
    var k = this.sortKey, d = this.sortDir;
    if (k) {
      rows.sort(function (a, b) {
        var x = a[k], y = b[k];
        var xn = (x === null || x === undefined || x === ""), yn = (y === null || y === undefined || y === "");
        if (xn && yn) return 0; if (xn) return 1; if (yn) return -1;
        if (typeof x === "number" && typeof y === "number") return d * (x - y);
        return d * String(x).localeCompare(String(y), undefined, { numeric: true });
      });
    }
    return rows;
  };
  AtlasTable.prototype.render = function () {
    var self = this, o = this.o, rows = this.visibleRows();
    this.ths.forEach(function (p) {
      p[1].setAttribute("aria-sort", p[0] === self.sortKey ? (self.sortDir > 0 ? "ascending" : "descending") : "none");
    });
    var frag = document.createDocumentFragment();
    rows.forEach(function (r) {
      var tr = document.createElement("tr");
      if (r.__id) tr.id = r.__id;
      o.columns.forEach(function (c) {
        var td = document.createElement("td");
        if (c.num) td.className = "num";
        var v = c.render ? c.render(r) : fmt(r[c.key], c.digits);
        if (v && v.nodeType) td.appendChild(v); else td.textContent = v === undefined ? "" : v;
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    });
    this.tbody.innerHTML = "";
    this.tbody.appendChild(frag);
    if (o.countEl) o.countEl.textContent = rows.length + " of " + o.rows.length + " rows shown";
    this.lastRows = rows;
  };
  AtlasTable.prototype.toTSV = function () {
    var cols = this.o.columns;
    var lines = [cols.map(function (c) { return c.key; }).join("\t")];
    (this.lastRows || []).forEach(function (r) {
      lines.push(cols.map(function (c) {
        var v = r[c.key]; return v === null || v === undefined ? "" : String(v).replace(/[\t\n]/g, " ");
      }).join("\t"));
    });
    return lines.join("\n") + "\n";
  };

  function downloadText(name, text) {
    var blob = new Blob([text], { type: "text/tab-separated-values" });
    var a = el("a", { href: URL.createObjectURL(blob), download: name });
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }

  /* Tooltip shared by SVG rows */
  var tip = null;
  function showTip(html, x, y) {
    if (!tip) { tip = el("div", { "class": "tooltip", role: "tooltip", id: "atlas-tooltip" }); document.body.appendChild(tip); }
    tip.innerHTML = html;
    tip.style.display = "block";
    var w = tip.offsetWidth, h = tip.offsetHeight, vw = window.innerWidth, vh = window.innerHeight;
    var left = x + 14, top = y + 14;
    if (left + w > vw - 8) left = Math.max(8, x - w - 14);
    if (top + h > vh - 8) top = Math.max(8, vh - h - 8);
    tip.style.left = left + "px"; tip.style.top = top + "px";
  }
  function hideTip() { if (tip) tip.style.display = "none"; }
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") hideTip(); });

  window.Atlas = { el: el, esc: esc, fmt: fmt, AtlasTable: AtlasTable, downloadText: downloadText,
    showTip: showTip, hideTip: hideTip };
})();
