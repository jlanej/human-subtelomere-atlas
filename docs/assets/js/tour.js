/* Guided read tour: chapters from tour/tour.json (website/tour_pages.py), one igv.js browser at a time. */
(function () {
  "use strict";
  var host = document.getElementById("igv"), status = document.getElementById("tour-status");
  var $ = function (id) { return document.getElementById(id); };
  var chapters = [], cur = -1, gi = 0, browser = null, token = 0;

  function setStatus(msg, isErr) { status.textContent = msg || ""; status.className = isErr ? "tour-error" : "muted"; }

  function steps() {
    var ol = $("tour-steps"); ol.innerHTML = "";
    chapters.forEach(function (c, i) {
      var li = document.createElement("li"), b = document.createElement("button");
      b.type = "button"; b.textContent = (i + 1) + ". " + c.short; b.title = c.title;
      b.addEventListener("click", function () { show(i); });
      if (i === cur) b.setAttribute("aria-current", "step");
      li.appendChild(b); ol.appendChild(li);
    });
    $("tour-prev").disabled = cur <= 0; $("tour-next").disabled = cur >= chapters.length - 1;
  }

  function viewButtons(views) {
    var box = $("tour-views"); box.innerHTML = "";
    if (!views || views.length < 2) return;
    views.forEach(function (v, i) {
      var b = document.createElement("button"); b.type = "button"; b.textContent = v.label;
      if (i === 0) b.setAttribute("aria-pressed", "true");
      b.addEventListener("click", function () {
        Array.prototype.forEach.call(box.children, function (x) { x.removeAttribute("aria-pressed"); });
        b.setAttribute("aria-pressed", "true");
        if (browser) browser.search(v.locus.join(" "));
      });
      box.appendChild(b);
    });
  }

  async function load(cfg) {
    var my = ++token;
    if (browser) { igv.removeBrowser(browser); browser = null; }
    host.innerHTML = "";
    setStatus("Loading molecules and reference ...");
    try {
      var b = await igv.createBrowser(host, {
        reference: cfg.ref, locus: cfg.views[0].locus, tracks: cfg.tracks,
        showCenterGuide: false, showCursorGuide: false, showSampleNames: false,
        showControls: true, showNavigation: true, showRuler: true, queryParametersSupported: false
      });
      if (my !== token) { igv.removeBrowser(b); return; }
      browser = b; setStatus("");
    } catch (e) { setStatus("Could not load this view: " + (e && e.message ? e.message : e), true); }
  }

  function texts(c, extra) {
    $("tour-title").textContent = c.title;
    $("tour-lead").innerHTML = c.lead || "";
    var ul = $("tour-notes"); ul.innerHTML = "";
    (c.notes || []).concat(extra || []).forEach(function (n) { var li = document.createElement("li"); li.innerHTML = n; ul.appendChild(li); });
  }

  function showGallery(c, j) {
    gi = Math.max(0, Math.min(c.gallery.length - 1, j));
    var g = c.gallery[gi]; $("tour-claim").value = String(gi);
    $("gal-page").href = g.page;
    texts(c, ["<b>" + g.label + "</b>: " + g.n_definitive + " definitive molecules (" + g.verdict.replace(/_/g, " ") +
              "); " + g.n_shown + " shown. Donor end " + g.donor + ", recipient end " + g.recipient + ". " + (g.note || "")]);
    viewButtons(g.views); history.replaceState(null, "", "#gallery/" + encodeURIComponent(g.slug));
    load(g);
  }

  function show(i, slug) {
    cur = Math.max(0, Math.min(chapters.length - 1, i)); var c = chapters[cur]; steps();
    var gal = $("tour-gallery");
    if (c.gallery) {
      gal.hidden = false; var sel = $("tour-claim");
      if (!sel.options.length) {
        var groups = {};
        c.gallery.forEach(function (g, k) {
          if (!groups[g.cohort]) { groups[g.cohort] = document.createElement("optgroup"); groups[g.cohort].label = g.cohort.replace(/_/g, " "); sel.appendChild(groups[g.cohort]); }
          var o = document.createElement("option"); o.value = String(k); o.textContent = g.label; groups[g.cohort].appendChild(o);
        });
      }
      var k = slug ? c.gallery.findIndex(function (g) { return g.slug === slug; }) : gi;
      showGallery(c, k < 0 ? 0 : k);
    } else {
      gal.hidden = true; texts(c); viewButtons(c.views);
      history.replaceState(null, "", "#" + c.id); load(c);
    }
  }

  function fromHash() {
    var h = decodeURIComponent(location.hash.replace(/^#/, ""));
    if (h.indexOf("gallery/") === 0) return [chapters.findIndex(function (c) { return c.gallery; }), h.slice(8)];
    var i = chapters.findIndex(function (c) { return c.id === h; });
    return [i < 0 ? 0 : i, null];
  }

  document.addEventListener("DOMContentLoaded", function () {
    if (typeof igv === "undefined") { setStatus("The genome browser (igv.js) did not load.", true); return; }
    fetch("tour/tour.json").then(function (r) { return r.json(); }).then(function (j) {
      chapters = j.chapters;
      $("tour-prev").addEventListener("click", function () { show(cur - 1); });
      $("tour-next").addEventListener("click", function () { show(cur + 1); });
      $("tour-claim").addEventListener("change", function (e) { showGallery(chapters[cur], parseInt(e.target.value, 10)); });
      $("gal-prev").addEventListener("click", function () { showGallery(chapters[cur], gi - 1); });
      $("gal-next").addEventListener("click", function () { showGallery(chapters[cur], gi + 1); });
      document.addEventListener("keydown", function (e) {
        var t = e.target && e.target.tagName;
        if (t === "INPUT" || t === "SELECT" || t === "TEXTAREA" || e.altKey || e.ctrlKey || e.metaKey) return;
        if (e.key === "ArrowRight") { show(cur + 1); } else if (e.key === "ArrowLeft") { show(cur - 1); }
      });
      var h = fromHash(); show(h[0], h[1]);
    }).catch(function (e) { setStatus("Could not load the tour: " + e, true); });
  });
})();
