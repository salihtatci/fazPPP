/* FazPPP arayuzu */
(() => {
  "use strict";
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));

  // ------------------------------------------------------------------ dil
  let lang = "tr";
  try { lang = localStorage.getItem("fazppp-lang") || (navigator.language || "tr").slice(0, 2); } catch (e) { /* yok */ }
  if (!I18N[lang]) lang = "tr";
  const t = (k) => (I18N[lang][k] ?? I18N.tr[k] ?? k);

  function applyLang() {
    document.documentElement.lang = lang;
    document.title = t("doc.title");
    $$("[data-i18n]").forEach((el) => { el.textContent = t(el.dataset.i18n); });
    $$(".lang button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.lang === lang)));
    renderSlots();
    renderMethod();
    syncKf();
    if (skyData) drawHeroSky(false);
    if (state.result) renderResults(state.result, state.activeTab);
  }
  $$(".lang button").forEach((b) => b.addEventListener("click", () => {
    lang = b.dataset.lang;
    try { localStorage.setItem("fazppp-lang", lang); } catch (e) { /* yok */ }
    applyLang();
  }));

  // ------------------------------------------------------------------ sayi bicimleri
  const dec = () => (lang === "tr" ? "," : ".");
  function num(v, d = 3) {
    if (v === null || v === undefined || Number.isNaN(v)) return "–";
    const s = Number(v).toFixed(d).replace("-", "−");
    return lang === "tr" ? s.replace(".", ",") : s;
  }
  function coord(v, d = 4) {
    if (v === null || v === undefined) return "–";
    const neg = v < 0;
    const [i, f] = Math.abs(v).toFixed(d).split(".");
    const g = i.replace(/\B(?=(\d{3})+(?!\d))/g, " ");
    return (neg ? "−" : "") + g + (f ? dec() + f : "");
  }
  const sgn = (v, d = 3) => (v === null || v === undefined ? "–" : (v > 0 ? "+" : v < 0 ? "−" : "") + num(Math.abs(v), d));
  const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  // ic uydu numarasi: GPS 1-32, Galileo 100 + PRN
  const svn = (p) => (p >= 100 ? "E" + String(p - 100).padStart(2, "0") : "G" + String(p).padStart(2, "0"));

  // ------------------------------------------------------------------ dosya yuvalari
  const SLOTS_MAIN = [
    { id: "obs", req: true, ext: ".o,.obs,.rnx,.gz,.*o" },
    { id: "nav", req: false, ext: ".n,.nav,.rnx,.gz,.*n" },
    { id: "sp3", req: true, multi: true, ext: ".sp3,.SP3,.gz,.eph" },
    { id: "clk", req: false, multi: true, ext: ".clk,.CLK,.gz" },
    { id: "met", req: false, ext: ".m,.met,.rnx,.gz,.*m" },
  ];
  const SLOTS_MORE = [
    { id: "atx", req: false, ext: ".atx,.gz" },
    { id: "dcb", req: false, ext: ".DCB,.dcb,.gz" },
  ];
  const state = { files: {}, result: null, activeTab: "summary", running: false, csv: {}, charts: {} };

  function renderSlots() {
    for (const [host, list] of [["#slots-main", SLOTS_MAIN], ["#slots-more", SLOTS_MORE]]) {
      const el = $(host);
      el.innerHTML = list.map((s) => {
        const f = state.files[s.id];
        const names = f ? f.map((x) => esc(x.name)).join(", ") : "";
        return `<div class="slot${f ? " filled" : ""}" data-slot="${s.id}">
          <div class="slot-head"><span class="slot-name">${esc(t("up." + s.id))}</span>
          <span class="slot-tag${s.req ? " req" : ""}">${esc(t(s.req ? "up.required" : "up.optional"))}</span></div>
          <p class="slot-hint">${esc(t("up." + s.id + ".hint"))}</p>
          <input type="file" ${s.multi ? "multiple" : ""} aria-label="${esc(t("up." + s.id))}">
          <div class="slot-file">${f ? `<span class="fname">${names}</span><button type="button" data-clear="${s.id}">${esc(t("up.clear"))}</button>`
            : `<span class="drop">${esc(t("up.drop"))}</span>`}</div></div>`;
      }).join("");
    }
    $$(".slot").forEach((slot) => {
      const id = slot.dataset.slot;
      const inp = $("input[type=file]", slot);
      inp.addEventListener("change", () => takeFiles(id, inp.files));
      slot.addEventListener("dragover", (e) => { e.preventDefault(); slot.classList.add("over"); });
      slot.addEventListener("dragleave", () => slot.classList.remove("over"));
      slot.addEventListener("drop", (e) => { e.preventDefault(); slot.classList.remove("over"); takeFiles(id, e.dataTransfer.files); });
    });
    $$("[data-clear]").forEach((b) => b.addEventListener("click", (e) => {
      e.preventDefault(); e.stopPropagation(); delete state.files[b.dataset.clear]; renderSlots();
    }));
  }

  async function readFileText(file) {
    const name = file.name;
    if (/\.Z$/.test(name)) throw new Error(t("err.Z"));
    let buf = await file.arrayBuffer();
    if (/\.gz$/i.test(name)) {
      const ds = new DecompressionStream("gzip");
      buf = await new Response(new Blob([buf]).stream().pipeThrough(ds)).arrayBuffer();
    }
    return new TextDecoder("iso-8859-1").decode(buf);
  }

  async function takeFiles(id, fileList) {
    const arr = Array.from(fileList || []);
    if (!arr.length) return;
    ensureWorker();
    const def = [...SLOTS_MAIN, ...SLOTS_MORE].find((s) => s.id === id);
    const use = def.multi ? arr : arr.slice(0, 1);
    try {
      const out = [];
      for (const f of use) out.push({ name: f.name, text: await readFileText(f) });
      state.files[id] = out;
      hideError();
    } catch (err) {
      showError("read", `${t("err.read")} ${err.message}`);
    }
    renderSlots();
  }

  // ------------------------------------------------------------------ worker
  let worker = null;
  let workerReady = false;
  function ensureWorker() {
    if (worker) return worker;
    worker = new Worker("worker.js");
    worker.onmessage = (ev) => {
      const m = ev.data;
      if (m.type === "boot") setStatus(t("boot." + m.stage) + "…");
      else if (m.type === "ready") { workerReady = true; if (!state.running) setStatus(t("boot.ready")); }
      else if (m.type === "progress") markStep(m.step);
      else if (m.type === "result") onResult(m.payload);
    };
    worker.onerror = (e) => { state.running = false; setBusy(false); showError("worker", e.message || ""); };
    worker.postMessage({ type: "init" });
    return worker;
  }
  const setStatus = (s) => { $("#status").textContent = s; };

  const STEP_KEYS = ["read", "preprocess", "orbit", "atmosphere", "corrections", "l3", "tables", "iono", "single", "multipath", "spp", "plots"];
  function renderSteps() {
    const ol = $("#steps");
    ol.innerHTML = STEP_KEYS.map((k) => `<li data-step="${k}">${esc(t("step." + k))}</li>`).join("");
    ol.hidden = false;
  }
  function markStep(key) {
    if (key === "done") { $$("#steps li").forEach((li) => { li.classList.remove("active"); li.classList.add("done"); }); return; }
    let seen = false;
    $$("#steps li").forEach((li) => {
      if (li.dataset.step === key) { li.classList.add("active"); li.classList.remove("done"); seen = true; }
      else if (!seen) { li.classList.remove("active"); li.classList.add("done"); }
    });
    setStatus(t("run.running") + ": " + t("step." + key));
  }
  function setBusy(b) { $("#run").disabled = b; $("#demo-btn").disabled = b; }

  function params() {
    const v = (id) => $(id).value.trim();
    const n = (id) => { const x = parseFloat(v(id).replace(",", ".")); return Number.isFinite(x) ? x : null; };
    const ref = [n("#p-rx"), n("#p-ry"), n("#p-rz")];
    return {
      epoch: v("#p-epoch") || null,
      cutoff: n("#p-cutoff") ?? 10,
      ref_xyz: ref.every((x) => x !== null) ? ref : null,
      sig_phase: n("#p-sp") ?? 0.01, sig_code: n("#p-sc") ?? 1.0,
      zwd_block_h: n("#p-zwd") ?? 2, spp_step: n("#p-spp") ?? 300,
      mode: $("#p-mode").value, systems: $("#p-sys").value,
      rate: n("#p-rate"), min_arc_min: n("#p-minarc") ?? 10,
      kf_dyn: $("#p-kfdyn").value, kf_qh: n("#p-kfqh"), kf_qv: n("#p-kfqv"), kf_smooth: $("#p-kfsm").checked,
    };
  }

  function run() {
    if (state.running) return;
    if (!state.files.obs || !state.files.sp3) { setStatus(t("run.needObs")); return; }
    hideError();
    state.running = true;
    setBusy(true);
    renderSteps();
    setStatus(workerReady ? t("run.running") + "…" : t("boot.pyodide") + "…");
    const one = (k) => (state.files[k] ? state.files[k][0].text : null);
    state.lastRun = { files: Object.fromEntries(Object.entries(state.files).map(([k, v]) => [k, v.map((f) => f.name).join(", ")])), params: params() };
    ensureWorker().postMessage({
      type: "run",
      files: { obs: one("obs"), sp3: state.files.sp3.map((f) => f.text), clk: state.files.clk ? state.files.clk.map((f) => f.text) : null,
        nav: one("nav"), met: one("met"), atx: one("atx"), dcb: one("dcb") },
      params: params(),
    });
  }
  $("#run").addEventListener("click", run);
  const syncKf = () => {
    $("#kf-box").hidden = $("#p-mode").value !== "kalman";
    const pv = $("#p-kfdyn").value === "pv";
    $("#p-kfqh").placeholder = (pv ? 0.5 : 1).toLocaleString(lang === "tr" ? "tr-TR" : "en-GB");
    $("#p-kfqv").placeholder = (pv ? 0.01 : 0.05).toLocaleString(lang === "tr" ? "tr-TR" : "en-GB");
  };
  $("#p-mode").addEventListener("change", syncKf);
  $("#p-kfdyn").addEventListener("change", syncKf);

  function onResult(payload) {
    state.running = false;
    setBusy(false);
    let r;
    try { r = JSON.parse(payload); } catch (e) { showError("internal", String(e)); return; }
    if (!r.ok) { showError(r.code, r.message, r.trace); setStatus(""); return; }
    markStep("done");
    setStatus("");
    state.result = r.result;
    renderResults(r.result, "summary");
    $("#results").hidden = false;
    $("#results").scrollIntoView({ behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth" });
  }

  function showError(code, msg, trace) {
    const box = $("#error");
    const known = I18N[lang]["err." + code];
    box.innerHTML = `<h3>${esc(t("err.title"))}</h3><p>${esc(known || msg || t("err.internal"))}</p>` +
      ((trace || (!known && msg)) ? `<details><summary>${esc(t("err.details"))}</summary><pre>${esc(trace || msg)}</pre></details>` : "");
    box.hidden = false;
  }
  const hideError = () => { $("#error").hidden = true; };

  // ------------------------------------------------------------------ ornek veri
  const DEMO = { dir: "examples/ista_2014_240/", obs: "ISTA2400.14o", nav: "ISTA2400.14n", sp3: "igs18074.sp3", met: "ISTA2400.14m",
    epoch: "13:15:00", ref: ["4208830.1293", "2334850.4702", "4171267.3316"] };
  $("#demo-btn").addEventListener("click", async () => {
    setBusy(true);
    ensureWorker();
    try {
      for (const k of ["obs", "nav", "sp3", "met"]) {
        const txt = await (await fetch(DEMO.dir + DEMO[k])).text();
        state.files[k] = [{ name: DEMO[k], text: txt }];
      }
      delete state.files.atx; delete state.files.dcb; delete state.files.clk;
      $("#p-mode").value = "static"; $("#p-sys").value = "GE"; $("#p-rate").value = ""; syncKf();
      $("#p-epoch").value = DEMO.epoch;
      [$("#p-rx").value, $("#p-ry").value, $("#p-rz").value] = DEMO.ref;
      renderSlots();
      $("#compute").scrollIntoView({ behavior: "smooth" });
      setBusy(false);
      run();
    } catch (e) { setBusy(false); showError("read", `${t("err.read")} ${e.message}`); }
  });

  // ------------------------------------------------------------------ giris gokyuzu grafigi
  let skyData = null;
  async function drawHeroSky(animate) {
    const svg = $("#sky");
    if (!skyData) {
      try { skyData = await (await fetch("assets/demo_sky.json")).json(); } catch (e) { return; }
    }
    const P = (az, el) => { const r = (90 - el) / 90 * 100; const a = az * Math.PI / 180; return [r * Math.sin(a), -r * Math.cos(a)]; };
    const card = lang === "tr" ? ["K", "D", "G", "B"] : ["N", "E", "S", "W"];
    let h = "";
    for (const el of [0, 30, 60]) h += `<circle class="sky-ring${el === 0 ? " outer" : ""}" r="${(90 - el) / 90 * 100}"/>`;
    for (let a = 0; a < 360; a += 30) { const [x, y] = P(a, 0); h += `<line class="sky-spoke" x1="0" y1="0" x2="${x.toFixed(1)}" y2="${y.toFixed(1)}"/>`; }
    [[0, card[0]], [90, card[1]], [180, card[2]], [270, card[3]]].forEach(([a, l]) => { const [x, y] = P(a, -9); h += `<text class="sky-lab" x="${x}" y="${y}">${l}</text>`; });
    h += `<text class="sky-el" x="2" y="${-(60 / 90 * 100) + 6}">30°</text><text class="sky-el" x="2" y="${-(30 / 90 * 100) + 6}">60°</text>`;
    const epochSv = new Set(skyData.epoch.map((e) => e.sv));
    let i = 0;
    for (const [sv, tr] of Object.entries(skyData.sky)) {
      let d = ""; let pen = false; let prev = null;
      for (let k = 0; k < tr.az.length; k++) {
        const az = tr.az[k], el = tr.el[k];
        if (az === null || el === null || el < 0) { pen = false; continue; }
        if (prev !== null) { const da = Math.abs(az - prev[0]); if (Math.min(da, 360 - da) > 30 || Math.abs(el - prev[1]) > 10) pen = false; }
        const [x, y] = P(az, el);
        d += (pen ? "L" : "M") + x.toFixed(1) + " " + y.toFixed(1);
        pen = true; prev = [az, el];
      }
      if (d) h += `<path class="sky-track${epochSv.has(sv) ? " lit" : ""}" d="${d}" style="--d:${(i++ * 0.035).toFixed(2)}s"/>`;
    }
    for (const e of skyData.epoch) {
      const [x, y] = P(e.az, Math.max(e.el, 0));
      h += `<circle class="sky-sat${e.el < 10 ? " low" : ""}" cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3.6"/>`;
      h += `<text class="sky-prn" x="${(x + 5).toFixed(1)}" y="${(y - 4).toFixed(1)}">${e.sv}</text>`;
    }
    svg.innerHTML = h;
    if (animate && !matchMedia("(prefers-reduced-motion: reduce)").matches) {
      $$(".sky-track", svg).forEach((p) => p.style.setProperty("--len", Math.ceil(p.getTotalLength())));
      svg.classList.add("sky-draw");
    }
  }

  function renderMethod() {
    const keys = ["read", "preprocess", "orbit", "atmosphere", "corrections", "l3", "iono", "single", "multipath", "spp"];
    $("#method-list").innerHTML = keys.map((k) => `<li>${esc(t("step." + k))}</li>`).join("");
  }

  // ------------------------------------------------------------------ tablolar
  function table(id, title, cols, rows, opts = {}) {
    state.csv[id] = { cols, rows };
    const head = cols.map((c) => `<th scope="col">${esc(t(c.l))}</th>`).join("");
    let shown = rows;
    if (opts.maxRows && rows.length > opts.maxRows) {
      const step = rows.length / opts.maxRows;
      shown = Array.from({ length: opts.maxRows }, (_, i) => rows[Math.floor(i * step)]);
    }
    const body = shown.map((r) => {
      const off = opts.offKey && r[opts.offKey] === false;
      return `<tr${off ? ' class="off"' : ""}>` + cols.map((c) => {
        const v = r[c.k];
        let s;
        if (c.f) s = c.f(v, r);
        else if (typeof v === "number") s = c.c ? coord(v, c.d ?? 3) : (c.sg ? sgn(v, c.d ?? 3) : num(v, c.d ?? 3));
        else if (typeof v === "boolean") s = t(v ? "col.yes" : "col.no");
        else s = esc(v ?? "–");
        return `<td${typeof v === "string" ? ' class="txt"' : ""}>${s}</td>`;
      }).join("") + "</tr>";
    }).join("");
    return `<div class="tbl-head"><h4>${esc(title)}</h4><button type="button" data-csv="${id}">${esc(t("csv"))}</button></div>
      <div class="tbl-wrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
  }
  document.addEventListener("click", (e) => {
    const b = e.target.closest("[data-csv]");
    if (!b) return;
    const { cols, rows } = state.csv[b.dataset.csv];
    const line = (a) => a.map((x) => { const s = String(x ?? ""); return /[",;\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; }).join(",");
    const txt = [line(cols.map((c) => t(c.l))), ...rows.map((r) => line(cols.map((c) => r[c.k])))].join("\n");
    download(`${b.dataset.csv}.csv`, txt, "text/csv");
  });
  function download(name, text, type) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type }));
    a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }
  $("#dl-json").addEventListener("click", () => {
    if (!state.result) return;
    const st = state.result.station.marker || "station";
    download(`fazppp_${st}_${state.result.epoch.used.replace(/:/g, "")}.json`, JSON.stringify(state.result, null, 1), "application/json");
  });

  // ------------------------------------------------------------------ PDF
  const JSPDF_URL = "https://cdn.jsdelivr.net/npm/jspdf@2.5.2/dist/jspdf.umd.min.js";
  const AUTOTABLE_URL = "https://cdn.jsdelivr.net/npm/jspdf-autotable@3.8.4/dist/jspdf.plugin.autotable.min.js";
  const scriptCache = {};
  function loadScript(src) {
    if (!scriptCache[src]) scriptCache[src] = new Promise((res, rej) => {
      const el = document.createElement("script"); el.src = src; el.onload = res; el.onerror = () => rej(new Error(src)); document.head.appendChild(el);
    });
    return scriptCache[src];
  }
  async function fontB64(url) {
    const buf = new Uint8Array(await (await fetch(url)).arrayBuffer());
    let bin = ""; for (let i = 0; i < buf.length; i += 0x8000) bin += String.fromCharCode.apply(null, buf.subarray(i, i + 0x8000));
    return btoa(bin);
  }
  const INK = [27, 39, 44], MUTED = [100, 117, 122], LINE = [200, 209, 204], PAPER = [238, 241, 236], ACCENT = [224, 165, 14];

  async function exportPdf() {
    const R = state.result;
    if (!R) return;
    const btn = $("#dl-pdf");
    const label = btn.textContent;
    btn.disabled = true; btn.textContent = t("res.pdfBusy");
    try {
      await loadScript(JSPDF_URL);
      await loadScript(AUTOTABLE_URL);
      if (!window.Plotly) await new Promise((r) => setTimeout(r, 1000));
      const { jsPDF } = window.jspdf;
      const doc = new jsPDF({ unit: "pt", format: "a4", compress: true });
      doc.addFileToVFS("DejaVuSans.ttf", await fontB64("assets/fonts/DejaVuSans.subset.ttf"));
      doc.addFont("DejaVuSans.ttf", "DejaVu", "normal");
      doc.addFileToVFS("DejaVuSans-Bold.ttf", await fontB64("assets/fonts/DejaVuSans-Bold.subset.ttf"));
      doc.addFont("DejaVuSans-Bold.ttf", "DejaVu", "bold");
      const W = doc.internal.pageSize.getWidth(), H = doc.internal.pageSize.getHeight(), M = 42, CW = W - 2 * M;
      let y = M;
      const ensure = (h) => { if (y + h > H - M - 14) { doc.addPage(); y = M; } };
      const text = (str, { size = 9, bold = false, color = INK, gap = 4, indent = 0 } = {}) => {
        if (!str) return;
        doc.setFont("DejaVu", bold ? "bold" : "normal"); doc.setFontSize(size); doc.setTextColor(...color);
        const lines = doc.splitTextToSize(String(str), CW - indent);
        const lh = size * 1.35;
        for (const ln of lines) { ensure(lh); doc.text(ln, M + indent, y + size); y += lh; }
        y += gap;
      };
      const heading = (str, size = 14) => {
        ensure(size * 2.6);
        y += size * 0.5;
        text(str, { size, bold: true, gap: 2 });
        doc.setDrawColor(...ACCENT); doc.setLineWidth(1.4); doc.line(M, y, M + 46, y); y += 10;
      };
      const tableOpts = (extra) => Object.assign({
        startY: y, margin: { left: M, right: M, bottom: M + 14 }, theme: "grid",
        styles: { font: "DejaVu", fontSize: 7.2, cellPadding: 2.6, textColor: INK, lineColor: LINE, lineWidth: 0.4, overflow: "linebreak" },
        headStyles: { fillColor: PAPER, textColor: INK, fontStyle: "bold" },
        didDrawPage: () => {},
      }, extra);
      const afterTable = () => { y = doc.lastAutoTable.finalY + 10; };
      const kvTable = (pairs) => {
        doc.autoTable(tableOpts({ body: pairs, theme: "plain", styles: { font: "DejaVu", fontSize: 8.4, cellPadding: 2.2, textColor: INK },
          columnStyles: { 0: { textColor: MUTED, cellWidth: CW * 0.42 }, 1: { fontStyle: "bold", halign: "right" } } }));
        afterTable();
      };
      const lightLayout = (lay) => {
        const L = JSON.parse(JSON.stringify(lay));
        L.paper_bgcolor = "#ffffff"; L.plot_bgcolor = "#ffffff";
        L.font = Object.assign({}, L.font, { color: "#1b272c", family: "DejaVu Sans, Arial, sans-serif" });
        for (const ax of ["xaxis", "yaxis"]) if (L[ax]) Object.assign(L[ax], { gridcolor: "#e3e8e5", linecolor: "#b9c4bf", zerolinecolor: "#b9c4bf", tickcolor: "#b9c4bf" });
        if (L.polar) { L.polar.bgcolor = "#ffffff"; ["radialaxis", "angularaxis"].forEach((k) => Object.assign(L.polar[k], { gridcolor: "#e3e8e5", linecolor: "#b9c4bf" })); }
        return L;
      };
      const chartImage = async (id, title, tall) => {
        const p = pending[id];
        if (!p || !window.Plotly) return;
        const w = 1000, h = tall ? 560 : 400;
        const L = lightLayout(p.layout);
        if (p.traces.filter((tr) => tr.showlegend !== false).length > 10) L.showlegend = false;
        L.width = w; L.height = h;
        const url = await Plotly.toImage({ data: p.traces, layout: L }, { format: "png", width: w, height: h, scale: 2 });
        const ih = CW * h / w;
        ensure(ih + 22);
        if (title) text(title, { size: 8.6, bold: true, color: [60, 75, 80], gap: 2 });
        doc.addImage(url, "PNG", M, y, CW, ih, undefined, "FAST");
        y += ih + 10;
      };

      // ---- kapak bilgisi
      const S = R.station, E = R.epoch, L3 = R.l3, run = state.lastRun || { files: {}, params: {} };
      doc.setFillColor(...ACCENT); doc.circle(M + 7, y + 8, 5, "F");
      doc.setFont("DejaVu", "bold"); doc.setFontSize(20); doc.setTextColor(...INK); doc.text(t("pdf.title"), M + 20, y + 15); y += 32;
      text(`${S.marker || "–"}  ·  ${E.used_full} GPST`, { size: 11, bold: true, gap: 2 });
      text(`${t("pdf.created")}: ${new Date().toLocaleString(lang === "tr" ? "tr-TR" : "en-GB")}   ·   https://salihtatci.github.io/fazPPP/`, { size: 8, color: MUTED, gap: 8 });
      const P = R.meta.params;
      kvTable([
        [t("pdf.files"), Object.entries(run.files).map(([k, v]) => `${t("up." + k)}: ${v}`).join("\n") || "–"],
        [t("pdf.params"), `${t("set.mode")}: ${t("set.mode." + (P.mode || "static"))}\n${t("set.systems")}: ${t("set.systems." + (P.systems || "G"))}\n${t("set.rate")}: ${P.rate ?? "–"}${P.kf_dyn ? `\n${t("set.kf.dyn")}: ${t("set.kf.dyn." + P.kf_dyn)}, q ${P.kf_qh} / ${P.kf_qv}, RTS ${P.kf_smooth ? "✓" : "–"}` : ""}\n${t("set.cutoff")}: ${P.cutoff}°\nσ faz/kod: ${P.sig_phase} / ${P.sig_code} m\n${t("set.zwd")}: ${P.zwd_block_h}\n${t("set.spp")}: ${P.spp_step}`],
        [t("pdf.ref"), R.meta.ref_source === "user" ? R.meta.ref_xyz.map((v) => coord(v, 4)).join("; ") + " m" : t("pdf.none")],
        [t("res.runtime"), num(R.meta.runtime_s, 1) + " s"],
      ]);

      // ---- sekmeler
      for (const k of TABS) {
        if (k !== "summary") { doc.addPage(); y = M; }
        heading(t("tab." + k), 15);
        const panel = $("#pn-" + k);
        if (!panel) continue;
        if (k === "summary") {
          const xyz = (X, s) => ["X", "Y", "Z"].map((a, i) => [a, coord(X[i], 4) + " m", s ? "± " + num(s[i], 1) + " mm" : ""]);
          text(R.kin ? `${S.marker}: ${t(L3.mode === "kalman" ? "sum.kfTrack" : "sum.kinTrack")}, ${E.used} (${t("sum.kinAt")})` : `${S.marker}: ${t("sum.static")}`, { size: 10, bold: true, gap: 2 });
          doc.autoTable(tableOpts({ body: xyz(L3.X, L3.sX_mm), theme: "plain", styles: { font: "DejaVu", fontSize: 11, cellPadding: 2.5, textColor: INK },
            columnStyles: { 0: { cellWidth: 24, textColor: MUTED }, 1: { fontStyle: "bold" }, 2: { textColor: MUTED, halign: "right" } } }));
          afterTable();
          text(`φ ${S.lat}    λ ${S.lon}    h ${num(S.h, 4)} m`, { size: 9, gap: 8 });
          if (L3.Xk) {
            text(`${t("sum.kin")}, ${E.used}`, { size: 10, bold: true, gap: 2 });
            doc.autoTable(tableOpts({ body: xyz(L3.Xk, L3.sXk_mm), theme: "plain", styles: { font: "DejaVu", fontSize: 10, cellPadding: 2.2, textColor: INK },
              columnStyles: { 0: { cellWidth: 24, textColor: MUTED }, 1: { fontStyle: "bold" }, 2: { textColor: MUTED, halign: "right" } } }));
            afterTable();
            text(`${t("l3.kinMinusStatic")}: ${neuStr(L3.kin_minus_static_neu_mm)}`, { size: 8.6, color: MUTED, gap: 8 });
          }
          kvTable([
            [`${t("sum.clock")}, ${E.used}`, `${num(L3.clock_ns, 1)} ns  (c·dt = ${num(L3.clock_m, 3)} ± ${num(L3.sclock_m, 3)} m)`],
            [`${t("sum.ztd")}, ${E.used}`, `${num(L3.ztd, 4)} m  (Saastamoinen ${num(L3.ztd_model, 4)} m)`],
            [t("sum.dop") + (L3.n_sat_kin ? ` (${L3.n_sat_kin} ${t("sum.sats")})` : ""), L3.dop && L3.dop.PDOP ? ["GDOP", "PDOP", "HDOP", "VDOP", "TDOP"].map((d) => `${d} ${num(L3.dop[d], 2)}`).join("  ") : "–"],
            ...(L3.isb_ns !== null && L3.isb_ns !== undefined ? [[t("sum.isb"), `${num(L3.isb_ns, 2)} ± ${num(L3.sisb_ns, 2)} ns`]] : []),
          ]);
        }
        const els = $$("p.lead2, dl.kv, h4, .tbl-wrap table, p.chart-title, .chart, p.note, .warns", panel)
          .filter((el) => !el.closest(".plate") && !(el.closest(".warns") && !el.classList.contains("warns")));
        let pendingTitle = null;
        for (const el of els) {
          if (el.matches("p.lead2")) text(el.textContent, { size: 8.8, color: [60, 75, 80], gap: 8 });
          else if (el.matches("p.note")) text(el.textContent, { size: 7.6, color: MUTED, gap: 6 });
          else if (el.matches("h4")) { ensure(60); text(el.textContent, { size: 10, bold: true, gap: 3 }); }
          else if (el.matches("dl.kv")) kvTable($$("dt", el).map((dt) => [dt.textContent, dt.nextElementSibling ? dt.nextElementSibling.textContent : ""]));
          else if (el.matches("table")) {
            const ncol = el.querySelectorAll("thead th").length;
            doc.autoTable(tableOpts({ html: el, styles: { font: "DejaVu", fontSize: ncol > 9 ? 6.2 : ncol > 6 ? 6.9 : 7.6, cellPadding: 2.4, textColor: INK, lineColor: LINE, lineWidth: 0.4, halign: "right" },
              columnStyles: { 0: { halign: "left" } },
              didParseCell: (d) => {
                if (d.column.index === 0) d.cell.styles.halign = "left";
                if (d.row.raw && d.row.raw.classList && d.row.raw.classList.contains("off")) d.cell.styles.textColor = MUTED;
              } }));
            afterTable();
          }
          else if (el.matches("p.chart-title")) pendingTitle = el.textContent;
          else if (el.matches(".chart")) { await chartImage(el.id.slice(3), pendingTitle, el.classList.contains("tall")); pendingTitle = null; }
          else if (el.matches(".warns")) {
            text(el.querySelector("h3").textContent, { size: 10, bold: true, gap: 2 });
            $$("li", el).forEach((li) => text("• " + li.textContent, { size: 8.4, color: [156, 74, 28], gap: 1, indent: 6 }));
            y += 6;
          }
        }
      }

      // ---- sayfa altliklari
      const n = doc.getNumberOfPages();
      for (let i = 1; i <= n; i++) {
        doc.setPage(i);
        doc.setFont("DejaVu", "normal"); doc.setFontSize(7); doc.setTextColor(...MUTED);
        doc.setDrawColor(...LINE); doc.setLineWidth(0.4); doc.line(M, H - M + 2, W - M, H - M + 2);
        doc.text(`FazPPP · ${S.marker || ""} · ${E.used_full}`, M, H - M + 13);
        doc.text(`${t("pdf.page")} ${i} / ${n}`, W - M, H - M + 13, { align: "right" });
      }
      doc.save(`fazppp_${S.marker || "station"}_${E.used.replace(/:/g, "")}.pdf`);
    } catch (err) {
      showError("internal", `${t("res.pdfFail")} ${err && err.message || err}`);
    } finally {
      btn.disabled = false; btn.textContent = label;
    }
  }
  $("#dl-pdf").addEventListener("click", exportPdf);

  const kv = (pairs) => `<dl class="kv">${pairs.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join("")}</dl>`;
  const neuStr = (a, d = 1, u = "mm") => (a ? a.map((x) => sgn(x, d)).join(" / ") + " " + u : "–");

  // ------------------------------------------------------------------ grafikler
  const PALETTE = ["#1b6f80", "#8a4777", "#b98500", "#3d7a3a", "#5b5fa8", "#a8532c", "#2f8f8a", "#7a6a2f", "#a33b5b", "#4d6b84", "#6f8f2f", "#8f5f9f"];
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  function baseLayout(extra = {}) {
    const ink = css("--ink-2"), grid = css("--line-2"), line = css("--line");
    return Object.assign({
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Archivo, Arial, sans-serif", size: 12, color: ink },
      margin: { l: 58, r: 16, t: 14, b: 44 },
      xaxis: { gridcolor: grid, linecolor: line, zerolinecolor: line, ticks: "outside", tickcolor: line },
      yaxis: { gridcolor: grid, linecolor: line, zerolinecolor: line, ticks: "outside", tickcolor: line },
      legend: { orientation: "h", y: -0.22, font: { size: 11 } },
      separators: lang === "tr" ? ", " : ". ",
      hovermode: "closest",
    }, extra);
  }
  const PCONF = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"] };
  const pending = {};
  function chart(id, title, traces, layout, tall) {
    pending[id] = { traces, layout };
    return `<p class="chart-title">${esc(title)}</p><div class="chart${tall ? " tall" : ""}" id="ch-${id}"></div>`;
  }
  function drawPending(panel) {
    if (!window.Plotly) { setTimeout(() => drawPending(panel), 300); return; }
    $$(".chart", panel).forEach((el) => {
      const id = el.id.slice(3);
      const p = pending[id];
      if (!p || el.dataset.drawn) return;
      if (p.traces.filter((tr) => tr.showlegend !== false).length > 10) p.layout.showlegend = false;
      Plotly.newPlot(el, p.traces, p.layout, PCONF);
      el.dataset.drawn = "1";
    });
  }
  const hoursAxis = (lab) => ({ title: { text: lab }, range: null });
  const tLab = () => (lang === "tr" ? "GPS zamanı [saat]" : "GPS time [h]");
  function bySat(tArr, yArr, prnArr, mode = "markers", size = 3) {
    const g = {};
    for (let i = 0; i < tArr.length; i++) {
      if (yArr[i] === null) continue;
      const k = svn(prnArr[i]);
      (g[k] = g[k] || { x: [], y: [] }).x.push(tArr[i]); g[k].y.push(yArr[i]);
    }
    return Object.entries(g).sort().map(([k, v], i) => ({ type: "scatter", mode, name: k, x: v.x, y: v.y,
      marker: { size, color: PALETTE[i % PALETTE.length] }, line: { width: 1.3, color: PALETTE[i % PALETTE.length] } }));
  }
  const epochLine = (h) => ({ type: "line", x0: h, x1: h, yref: "paper", y0: 0, y1: 1, line: { color: css("--accent"), width: 2, dash: "dot" } });

  // ------------------------------------------------------------------ sonuclar
  const TABS_BASE = ["summary", "pre", "orbit", "atm", "corr", "l3", "iono", "single", "mp", "spp"];
  let TABS = TABS_BASE.slice();
  function renderResults(R, active) {
    Object.keys(pending).forEach((k) => delete pending[k]);
    TABS = R.kin ? ["summary", "kin", ...TABS_BASE.slice(1)] : TABS_BASE.slice();
    if (!TABS.includes(active)) active = "summary";
    $("#runtime").textContent = `${t("res.runtime")}: ${num(R.meta.runtime_s, 1)} s`;
    $("#tabs").innerHTML = TABS.map((k) => `<button type="button" role="tab" id="tab-${k}" aria-controls="pn-${k}" aria-selected="${k === active}" data-tab="${k}">${esc(t("tab." + k))}</button>`).join("");
    const html = { summary: pSummary(R), pre: pPre(R), orbit: pOrbit(R), atm: pAtm(R), corr: pCorr(R), l3: pL3(R),
      iono: pIono(R), single: pSingle(R), mp: pMp(R), spp: pSpp(R), kin: R.kin ? pKin(R) : "" };
    $("#panels").innerHTML = TABS.map((k) => `<section class="panel" role="tabpanel" id="pn-${k}" aria-labelledby="tab-${k}" ${k === active ? "" : "hidden"}>${html[k]}</section>`).join("");
    $$("#tabs button").forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
    $("#tabs").addEventListener("keydown", (e) => {
      if (!["ArrowRight", "ArrowLeft"].includes(e.key)) return;
      const i = TABS.indexOf(state.activeTab) + (e.key === "ArrowRight" ? 1 : -1);
      const k = TABS[(i + TABS.length) % TABS.length]; showTab(k); $("#tab-" + k).focus();
    });
    showTab(active);
  }
  function showTab(k) {
    state.activeTab = k;
    $$("#tabs button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.tab === k)));
    $$(".panel").forEach((p) => { p.hidden = p.id !== "pn-" + k; });
    drawPending($("#pn-" + k));
  }

  function solName(n) { return t("sol." + n); }

  function pSummary(R) {
    const L = R.l3, S = R.station, E = R.epoch;
    const xyz = (X, s) => ["X", "Y", "Z"].map((a, i) => `<span class="ax">${a}</span><span class="v">${coord(X[i], 4)}</span><span class="s">${s ? "± " + num(s[i], 1) + " mm" : ""}</span>`).join("");
    const kinH = L.Xk ? `<div class="kin"><h3>${esc(t("sum.kin"))}, ${esc(E.used)}</h3><div class="coord">${xyz(L.Xk, L.sXk_mm)}</div>
      <p class="geo">${esc(t("l3.kinMinusStatic"))}: ${neuStr(L.kin_minus_static_neu_mm)}</p></div>` : "";
    const dop = L.dop && L.dop.PDOP ? `<div class="dops">${["GDOP", "PDOP", "HDOP", "VDOP", "TDOP"].map((k) => `<span>${k} <b>${num(L.dop[k], 2)}</b></span>`).join("")}</div>
      <p class="note">${L.n_sat_kin} ${esc(t("sum.sats"))}</p>` : "–";
    const rows = R.solutions.filter((s) => s.X).map((s) => ({ name: solName(s.name), X: s.X[0], Y: s.X[1], Z: s.X[2],
      N: s.neu?.N ?? null, E: s.neu?.E ?? null, U: s.neu?.U ?? null, H: s.neu?.H ?? null, D3: s.neu?.D3 ?? null }));
    const cols = [{ k: "name", l: "col.solution" }, { k: "X", l: "col.X", d: 4 }, { k: "Y", l: "col.Y", d: 4 }, { k: "Z", l: "col.Z", d: 4 },
      { k: "N", l: "col.N", d: 3, sg: 1 }, { k: "E", l: "col.E", d: 3, sg: 1 }, { k: "U", l: "col.U", d: 3, sg: 1 }, { k: "H", l: "col.H", d: 3 }, { k: "D3", l: "col.D3", d: 3 }];
    const warns = (R.meta.warnings || []).map((w) => {
      const [k, extra] = w.split(":");
      return `<li>${esc(t("warn." + k))}${extra ? " " + esc(extra) : ""}</li>`;
    }).join("");
    return `
      <div class="plate">
        <div class="plate-main">
          <h3>${esc(S.marker || t("sum.station"))}: ${R.kin ? esc(t(L.mode === "kalman" ? "sum.kfTrack" : "sum.kinTrack")) + ", " + esc(E.used) + " (" + esc(t("sum.kinAt")) + ")" : esc(t("sum.static"))}</h3>
          <div class="coord">${xyz(L.X, L.sX_mm)}</div>
          <p class="geo">φ ${esc(S.lat)} &nbsp; λ ${esc(S.lon)} &nbsp; h ${num(S.h, 4)} m</p>
          ${kinH}
        </div>
        <div class="plate-side">
          <div><h3>${esc(t("sum.clock"))}, ${esc(E.used)}</h3><div class="big-num">${num(L.clock_ns, 1)} ns</div>
            <p class="note">c·dt = ${num(L.clock_m, 3)} m ± ${num(L.sclock_m, 3)} m</p></div>
          <div><h3>${esc(t("sum.ztd"))}, ${esc(E.used)}</h3><div class="big-num">${num(L.ztd, 4)} m</div>
            <p class="note">Saastamoinen: ${num(L.ztd_model, 4)} m</p></div>
          <div><h3>${esc(t("sum.dop"))}</h3>${dop}</div>
          ${L.isb_ns !== null && L.isb_ns !== undefined ? `<div><h3>${esc(t("sum.isb"))}</h3><div class="big-num">${num(L.isb_ns, 2)} ns</div>
            <p class="note">± ${num(L.sisb_ns, 2)} ns · GPS ${L.n_gps_used} / Galileo ${L.n_gal_used} ${esc(t("sum.sats"))}</p></div>` : ""}
        </div>
      </div>
      ${table("solutions", t("sum.compare"), cols, rows)}
      <p class="note">${esc(t("sum.compare.lead." + (R.meta.ref_source || "static")))}</p>
      ${warns ? `<div class="warns"><h3>${esc(t("sum.warnings"))}</h3><ul>${warns}</ul></div>` : ""}`;
  }

  function pPre(R) {
    const S = R.station, P = R.pre, pl = R.plots;
    const info = kv([
      [t("pre.receiver"), esc(S.receiver)], [t("pre.antenna"), esc(S.antenna)],
      [t("pre.approx"), S.approx.map((v) => coord(v, 4)).join("; ") + " m"],
      [t("pre.arp"), S.delta_hen.map((v) => num(v, 4)).join(" / ") + " m"],
      [t("pre.version"), esc(String(S.rinex_version))], [t("pre.types"), esc(S.obs_types.join(" "))],
      [t("pre.codes"), esc(Object.entries(S.used_codes).filter(([, v]) => v).map(([k, v]) => `${k}=${v}`).join(", "))],
      ...(S.used_codes_E ? [[t("pre.codesE"), esc(Object.entries(S.used_codes_E).filter(([, v]) => v).map(([k, v]) => `${k}=${v}`).join(", "))]] : []),
      [t("pre.interval"), num(S.interval, 0) + " s"], [t("pre.span"), esc(S.first + " – " + S.last.slice(11))],
      [t("pre.counts"), `${S.n_epochs} / ${S.n_sats}`],
    ]);
    const st = kv([
      [t("pre.records"), P.n_records], [t("pre.frozen"), P.n_frozen], [t("pre.incomplete"), P.n_incomplete],
      [t("pre.arcs"), P.n_arcs], [t("pre.usableArcs"), P.n_usable_arcs], [t("pre.outliers"), P.n_outlier],
      [t("pre.usable"), P.n_usable], [t("pre.thr"), `GF ${num(P.gf_thr, 2)} m, MW ${num(P.mw_thr, 0)}`],
      ...(P.n_win_slips ? [[t("pre.winSlips"), P.n_win_slips]] : []),
    ]);
    const jt = P.jumps.length
      ? table("jumps", t("pre.jumps"), [{ k: "t", l: "col.time" }, { k: "median_m", l: "col.median", d: 3 }, { k: "nsat", l: "col.nsat", d: 0 }, { k: "k_ms", l: "col.k", d: 0 }], P.jumps)
        + `<p class="note">${esc(t("pre.jumps.lead"))}</p>`
      : `<h4>${esc(t("pre.jumps"))}</h4><p class="note">${esc(t("pre.jumps.none"))}</p>`;
    const at = table("arcs", t("pre.arcsEpoch"), [{ k: "sv", l: "col.sv" }, { k: "start", l: "col.start" }, { k: "end", l: "col.end" },
      { k: "n", l: "col.n", d: 0 }, { k: "min", l: "col.min", d: 1 }, { k: "mw_std", l: "col.mwstd", d: 3 }, { k: "usable", l: "col.used" }], P.arcs_epoch, { offKey: "usable" });
    let charts = "";
    if (pl.clockjump) {
      const c = pl.clockjump;
      charts += `<div class="grid2"><div>${chart("cj1", `${c.sv}: ${lang === "tr" ? "ham" : "raw"} λ₁L1 − P1 [km]`,
        [{ type: "scatter", mode: "markers", x: c.t, y: c.raw_km, marker: { size: 3, color: css("--l2") }, name: "raw" }],
        baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, shapes: [epochLine(pl.epoch_h)], showlegend: false }))}</div>
        <div>${chart("cj2", `${c.sv}: ${lang === "tr" ? "düzeltilmiş" : "corrected"} λ₁L1 − P1 [m]`,
        [{ type: "scatter", mode: "markers", x: c.t, y: c.cor_m, marker: { size: 3, color: css("--l1") }, name: "cor" }],
        baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, shapes: [epochLine(pl.epoch_h)], showlegend: false }))}</div></div>`;
    }
    const gf = Object.entries(pl.gfmw || {});
    if (gf.length) {
      const tr = (key) => gf.map(([sv, d], i) => ({ type: "scatter", mode: "lines", name: sv, x: d.t, y: d[key], line: { width: 1.3, color: PALETTE[i % PALETTE.length] } }));
      charts += `<div class="grid2"><div>${chart("gf", "GF − GF₀ [m]", tr("gf"), baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, shapes: [epochLine(pl.epoch_h)] }))}</div>
        <div>${chart("mw", lang === "tr" ? "MW − ortalama [devir]" : "MW − mean [cycles]", tr("mw"), baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, shapes: [epochLine(pl.epoch_h)] }))}</div></div>`;
    }
    return `<div class="grid2"><div><h4>${esc(t("pre.station"))}</h4>${info}</div><div><h4>${esc(t("pre.stats"))}</h4>${st}<p class="note">${esc(t("pre.method"))}</p></div></div>
      ${jt}${at}${charts}`;
  }

  function pOrbit(R) {
    const O = R.orbit, pl = R.plots;
    const cols = [{ k: "sv", l: "col.sv" }, { k: "X", l: "col.X", d: 3 }, { k: "Y", l: "col.Y", d: 3 }, { k: "Z", l: "col.Z", d: 3 },
      { k: "tau_ms", l: "col.tau", d: 3 }, { k: "dts_us", l: "col.dts", d: 3 }, { k: "rel_ns", l: "col.rel", d: 3 }, { k: "sagnac", l: "col.sagnac", d: 3 },
      { k: "el", l: "col.el", d: 2 }, { k: "az", l: "col.az", d: 2 }, { k: "used", l: "col.used" }];
    const ck = O.check ? kv([[t("orb.spread") + " P1 − ρ + c·dtˢ", num(O.check.spread_P1, 2) + " m"],
      ["P3 − ρ + c·dtˢ", num(O.check.spread_P3, 2) + " m"], ["P3 − ρ + c·dtˢ − T", num(O.check.spread_P3T, 2) + " m"],
      [t("orb.clockRough"), num(O.check.clock_rough_ms, 4) + " ms"]]) : "";
    const sky = Object.entries(pl.sky).map(([sv, d], i) => {
      const r = [], th = [];
      for (let k = 0; k < d.az.length; k++) {
        if (k > 0 && d.az[k] !== null && d.az[k - 1] !== null) {
          const da = Math.abs(d.az[k] - d.az[k - 1]); const dd = Math.min(da, 360 - da);
          if (dd > 30 || Math.abs(d.el[k] - d.el[k - 1]) > 10) { r.push(null); th.push(null); }
        }
        r.push(d.el[k] === null ? null : 90 - d.el[k]); th.push(d.az[k]);
      }
      return { type: "scatterpolar", mode: "lines", name: sv, r, theta: th, connectgaps: false,
        line: { width: 1, color: PALETTE[i % PALETTE.length] }, opacity: 0.55, hoverinfo: "name" };
    });
    sky.push({ type: "scatterpolar", mode: "markers+text", name: R.epoch.used, r: pl.sky_epoch.map((e) => 90 - e.el), theta: pl.sky_epoch.map((e) => e.az),
      text: pl.sky_epoch.map((e) => e.sv), textposition: "top right", marker: { size: 10, color: pl.sky_epoch.map((e) => (e.used ? css("--accent") : css("--sheet"))), line: { color: css("--ink"), width: 1.2 } } });
    const pol = baseLayout({ showlegend: false, margin: { l: 30, r: 30, t: 30, b: 30 },
      polar: { bgcolor: "rgba(0,0,0,0)", radialaxis: { range: [0, 90], tickvals: [30, 60, 90], ticktext: ["60°", "30°", "0°"], gridcolor: css("--line-2"), linecolor: css("--line") },
        angularaxis: { direction: "clockwise", rotation: 90, gridcolor: css("--line-2"), linecolor: css("--line") } } });
    return `<p class="lead2">${esc(t("orb.lead"))}</p>
      ${table("orbit", `${R.epoch.used}`, cols, O.sats, { offKey: "used" })}
      <div class="grid2"><div><h4>${esc(t("orb.check"))}</h4>${ck}<p class="note">${esc(t("orb.check.lead"))}</p></div>
      <div>${chart("sky", t("orb.chart.sky"), sky, pol, true)}</div></div>`;
  }

  function pAtm(R) {
    const D = R.data, A = R.atm, pl = R.plots;
    const met = kv([[t("atm.met"), `${num(D.met.P, 1)} mbar, ${num(D.met.T, 1)} °C, ${num(D.met.RH, 1)} % (${esc(t("atm.met.source." + D.met.source))})`],
      ["ZHD", num(D.zhd, 4) + " m"], ["ZWD", num(D.zwd, 4) + " m"], ["ZTD", num(D.zhd + D.zwd, 4) + " m"]]);
    const cols = [{ k: "sv", l: "col.sv" }, { k: "el", l: "col.el", d: 2 }, { k: "mh", l: "col.mh", d: 3 }, { k: "mw", l: "col.mw", d: 3 },
      { k: "T", l: "col.T", d: 3 }, { k: "I1_code", l: "col.I1code", d: 3 }, { k: "used", l: "col.used" }];
    const tr = [{ type: "scatter", mode: "markers", x: pl.trop.el, y: pl.trop.T, marker: { size: 3, color: css("--l1"), opacity: .6 }, name: "T" }];
    const tide = [{ type: "scatter", mode: "lines", x: pl.tide.t, y: pl.tide.up_cm, line: { color: css("--l2"), width: 2 }, name: "Up" }];
    return `<p class="lead2">${esc(t("atm.lead"))}</p>${met}
      ${table("atm", R.epoch.used, cols, A.rows, { offKey: "used" })}
      <div class="grid2"><div>${chart("trop", t("atm.chart.trop"), tr, baseLayout({ showlegend: false, xaxis: { ...baseLayout().xaxis, title: { text: t("col.el") } }, yaxis: { ...baseLayout().yaxis, title: { text: "m" } } }))}</div>
      <div>${chart("tide", t("atm.chart.tide"), tide, baseLayout({ showlegend: false, xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, shapes: [epochLine(pl.epoch_h)] }))}</div></div>`;
  }

  function pCorr(R) {
    const C = R.corr, L = R.l3, A = R.atm;
    const cols = [{ k: "sv", l: "col.sv" }, { k: "el", l: "col.el", d: 1 }, { k: "nadir", l: "col.nadir", d: 1 }, { k: "rec", l: "col.rec", d: 1 },
      { k: "sat", l: "col.sat", d: 1 }, { k: "tide", l: "col.tide", d: 1 }, { k: "shapiro", l: "col.shapiro", d: 1 }, { k: "windup", l: "col.windup", d: 1 },
      { k: "total", l: "col.total", d: 1 }, { k: "used", l: "col.used" }];
    const info = kv([
      [t("corr.recPco"), R.meta.rec_ant_found ? neuStr(A.rec_pco_l3_mm, 2) : esc(t("corr.recMissing"))],
      [t("corr.tide"), sgn(A.tide_up_cm, 1) + " cm (" + num(A.tide_range_cm[0], 1) + " … " + num(A.tide_range_cm[1], 1) + ")"],
      [t("corr.effect"), neuStr(L.no4b_static_minus_full_neu_mm)],
      [t("l3.rms") + " (" + (lang === "tr" ? "düzeltmesiz" : "without") + ")", num(L.no4b_rms_phase_mm, 1) + " mm → " + num(L.rms_phase_mm, 1) + " mm"],
    ]);
    return `<p class="lead2">${esc(t("corr.lead"))}</p>${info}${table("corr", R.epoch.used + " [mm]", cols, C.rows, { offKey: "used" })}`;
  }

  function pL3(R) {
    const L = R.l3, pl = R.plots;
    const D = R.data;
    const info = kv([
      [t("l3.mode"), esc(t("set.mode." + (L.mode || "static"))) + (R.meta.params.kf_dyn ? " · " + esc(t("set.kf.dyn." + R.meta.params.kf_dyn)) : "")],
      [t("l3.clockSrc"), `${esc(t("l3.clockSrc." + (D.sat_clock || "sp3")))}, ${num(D.rate, 0)} s`],
      [t("l3.sats"), `${L.n_gps_used ?? "–"} / ${L.n_gal_used ?? 0}`],
      ...(L.isb_ns !== null && L.isb_ns !== undefined ? [[t("l3.isb"), `${num(L.isb_ns, 2)} ± ${num(L.sisb_ns, 2)} ns`]] : []),
      [t("l3.obs"), `${L.n_obs} / ${L.n_par}`], [t("l3.epochs"), `${L.n_epochs}, ${L.n_arcs}, ${L.n_zwd}`],
      [t("l3.s0"), num(L.s0, 3)], [t("l3.rms"), `${num(L.rms_phase_mm, 1)} mm / ${num(L.rms_code_m, 2)} m`],
      ...(L.iter_mm && L.iter_mm.length ? [[t("l3.iter"), L.iter_mm.map((x) => num(x, 1)).join(", ") + " mm"]] : []),
      [t("l3.kinMinusStatic"), neuStr(L.kin_minus_static_neu_mm)],
    ]);
    const res = bySat(pl.l3.t, pl.l3.v_mm, pl.l3.prn, "markers", 4);
    const clk = [{ type: "scatter", mode: "lines+markers", x: pl.l3.clk_t, y: pl.l3.clk_us, line: { color: css("--l1"), width: 1.5 }, marker: { size: 4 }, name: "clock" }];
    const zwd = [{ type: "scatter", mode: "lines", line: { shape: "hv", color: css("--l2"), width: 2 }, x: [...pl.l3.zwd_t, pl.l3.zwd_t[pl.l3.zwd_t.length - 1] + R.meta.params.zwd_block_h], y: [...pl.l3.zwd_mm, pl.l3.zwd_mm[pl.l3.zwd_mm.length - 1]], name: "ZWD" }];
    const xa = { ...baseLayout().xaxis, title: { text: tLab() } };
    return `<p class="lead2">${esc(t(L.mode === "kalman" ? "l3.lead.kalman" : "l3.lead"))}</p>${info}
      ${chart("l3res", t("l3.chart.res"), res, baseLayout({ xaxis: xa, shapes: [epochLine(pl.epoch_h)] }))}
      <div class="grid2"><div>${chart("l3clk", t("l3.chart.clk"), clk, baseLayout({ xaxis: xa, showlegend: false, shapes: [epochLine(pl.epoch_h)] }))}</div>
      <div>${chart("l3zwd", t("l3.chart.zwd"), zwd, baseLayout({ xaxis: xa, showlegend: false, shapes: [epochLine(pl.epoch_h)] }))}</div></div>`;
  }

  function pIono(R) {
    const I = R.iono, pl = R.plots;
    if (!I || !I.available) return `<p class="lead2">${esc(t("iono.lead"))}</p><p>${esc(t("iono.unavailable"))}</p>`;
    const info = kv([
      [t("iono.station"), num(I.vtec_station, 2) + " TECU"], [t("iono.rate"), sgn(I.vtec_rate, 2) + " TECU/h"],
      [t("iono.dcbr"), num(I.dcb_rec_ns, 2) + " ns"], [t("iono.dcbs"), esc(t("iono.dcbs." + (I.sat_dcb === "file" ? "file" : "estimated")))],
      [t("iono.coderms"), num(I.code_rms_m, 3) + " m"], [t("iono.modelrms"), num(I.model_rms_tecu, 2) + " TECU"],
      ...(I.klob_minus ? [[t("iono.klob"), sgn(I.klob_minus.mean, 2) + " / " + num(I.klob_minus.rms, 2) + " m"]] : []),
    ]);
    const cols = [{ k: "sv", l: "col.sv" }, { k: "el", l: "col.el", d: 2 }, { k: "ipp_lat", l: "col.ipplat", d: 3 }, { k: "ipp_lon", l: "col.ipplon", d: 3 },
      { k: "stec", l: "col.stec", d: 2 }, { k: "vtec", l: "col.vtec", d: 2 }, { k: "I1", l: "col.I1", d: 3 }, { k: "I2", l: "col.I2", d: 3 }];
    if (I.rows.length && I.rows[0].I1_klob !== undefined) cols.push({ k: "I1_klob", l: "col.I1klob", d: 3 });
    const xa = { ...baseLayout().xaxis, title: { text: tLab() } };
    let charts = "";
    if (pl.iono) {
      const v = bySat(pl.iono.t, pl.iono.vtec, pl.iono.prn, "markers", 3);
      const p4 = [{ type: "scatter", mode: "markers", x: pl.iono.t, y: pl.iono.p4, marker: { size: 2, color: css("--line") }, name: "P2 − P1", showlegend: false },
        ...bySat(pl.iono.t, pl.iono.lev, pl.iono.prn, "markers", 2).map((tr) => ({ ...tr, showlegend: false }))];
      charts = `<div class="grid2"><div>${chart("vtec", t("iono.chart.vtec"), v, baseLayout({ xaxis: xa, yaxis: { ...baseLayout().yaxis, title: { text: "TECU" } }, shapes: [epochLine(pl.epoch_h)] }))}</div>
        <div>${chart("lev", t("iono.chart.lev"), p4, baseLayout({ xaxis: xa, yaxis: { ...baseLayout().yaxis, title: { text: "m" } }, shapes: [epochLine(pl.epoch_h)] }))}</div></div>`;
    }
    const a5cols = [{ k: "sv", l: "col.sv" }, { k: "P1", l: "col.P1", d: 3 }, { k: "P2", l: "col.P2", d: 3 }, { k: "P3", l: "col.P3", d: 3 },
      { k: "L1", l: "col.L1", d: 3 }, { k: "L2", l: "col.L2", d: 3 }, { k: "L3", l: "col.L3", d: 3 }, { k: "used", l: "col.used" }];
    return `<p class="lead2">${esc(t("iono.lead"))}</p>${info}${table("iono", R.epoch.used, cols, I.rows, { offKey: "used" })}
      ${charts}${table("a5", "P3, L3 – " + R.epoch.used, a5cols, R.a5 || [], { offKey: "used" })}`;
  }

  function pSingle(R) {
    const rows = R.single.map((s) => ({ freq: s.freq, model: t("model." + s.model), rms: s.rms_phase_mm,
      N: s.neu_kin?.N ?? null, E: s.neu_kin?.E ?? null, U: s.neu_kin?.U ?? null, kin: s.neu_kin?.D3 ?? null, stat: s.neu_static?.D3 ?? null,
      clock: s.clock_ns, ztd: s.ztd }));
    const cols = [{ k: "freq", l: "col.freq" }, { k: "model", l: "col.model" }, { k: "rms", l: "col.rms", d: 1 },
      { k: "N", l: "col.N", d: 3, sg: 1 }, { k: "E", l: "col.E", d: 3, sg: 1 }, { k: "U", l: "col.U", d: 3, sg: 1 },
      { k: "kin", l: "col.kin3d", d: 3 }, { k: "stat", l: "col.stat3d", d: 3 }, { k: "clock", l: "col.clock", d: 1 }, { k: "ztd", l: "col.ztd", d: 4 }];
    const noK = !R.data.klobuchar ? `<p class="note">${esc(t("single.noKlob"))}</p>` : "";
    return `<p class="lead2">${esc(t("single.lead"))}</p>${noK}${table("single", R.epoch.used, cols, rows)}
      <p class="note">${esc(t("sum.compare.lead." + (R.meta.ref_source || "static")))}</p>`;
  }

  function pMp(R) {
    const M = R.mp, pl = R.plots;
    const info = kv([[t("mp.epochRms"), `MP1 ${num(M.epoch_rms[0], 3)} m, MP2 ${num(M.epoch_rms[1], 3)} m`],
      [t("mp.cutRms"), `MP1 ${num(M.cut_rms[0], 3)} m, MP2 ${num(M.cut_rms[1], 3)} m`],
      [t("mp.allRms"), `MP1 ${num(M.all_rms[0], 3)} m, MP2 ${num(M.all_rms[1], 3)} m`]]);
    const cols = [{ k: "sv", l: "col.sv" }, { k: "el", l: "col.el", d: 2 }, { k: "az", l: "col.az", d: 2 }, { k: "MP1", l: "col.MP1", d: 3, sg: 1 }, { k: "MP2", l: "col.MP2", d: 3, sg: 1 }];
    const bcols = [{ k: "bin", l: "col.bin" }, { k: "MP1", l: "col.MP1", d: 3 }, { k: "MP2", l: "col.MP2", d: 3 }, { k: "n", l: "col.count", d: 0 }];
    const sc = [{ type: "scatter", mode: "markers", name: "MP1", x: pl.mp.el, y: pl.mp.mp1, marker: { size: 3, color: css("--l1"), opacity: .5 } },
      { type: "scatter", mode: "markers", name: "MP2", x: pl.mp.el, y: pl.mp.mp2, marker: { size: 3, color: css("--l2"), opacity: .5 } }];
    const bar = [{ type: "bar", name: "MP1", x: M.per_sat.map((r) => r.sv), y: M.per_sat.map((r) => r.MP1), marker: { color: css("--l1") } },
      { type: "bar", name: "MP2", x: M.per_sat.map((r) => r.sv), y: M.per_sat.map((r) => r.MP2), marker: { color: css("--l2") } }];
    return `<p class="lead2">${esc(t("mp.lead"))}</p>${info}
      <div class="grid2"><div>${table("mp", R.epoch.used, cols, M.rows)}</div><div>${table("mpbins", t("mp.bins"), bcols, M.bins)}</div></div>
      <div class="grid2"><div>${chart("mpsc", t("mp.chart.scatter"), sc, baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: t("col.el") } }, yaxis: { ...baseLayout().yaxis, range: [-3, 3], title: { text: "m" } } }))}</div>
      <div>${chart("mpbar", t("mp.chart.persat"), bar, baseLayout({ barmode: "group", yaxis: { ...baseLayout().yaxis, title: { text: "RMS [m]" } } }))}</div></div>`;
  }


  function pKin(R) {
    const K = R.kin, S = K.summary, P = K.plot, pl = R.plots;
    const info = kv([
      [t("kin.epochs"), `${S.n_epochs}, ${num(S.rate, 0)} s (${esc(S.first)} – ${esc(S.last)})`],
      [t("kin.span"), `${num(S.span_N, 2)} / ${num(S.span_E, 2)} / ${num(S.span_U, 2)} m`],
      [t("kin.path"), num(S.path_m, 1) + " m"],
      [t("kin.speed"), `${num(S.v_mean_moving, 2)} / ${num(S.v_max, 2)} m/s`],
      [t("kin.moving"), num((S.moving_frac ?? 0) * 100, 0) + " %"],
      [t("kin.sig"), (S.sig_med_mm || []).map((x) => num(x, 1)).join(" / ") + " mm"],
      [t("kin.origin"), R.meta.ref_xyz.map((v) => coord(v, 3)).join("; ") + " m"],
    ]);
    const xa = { ...baseLayout().xaxis, title: { text: tLab() } };
    const track = [{ type: "scatter", mode: "lines+markers", x: P.E, y: P.N, text: P.t.map((h) => { const s = Math.round(h * 3600); return `${String(Math.floor(s / 3600)).padStart(2, "0")}:${String(Math.floor(s % 3600 / 60)).padStart(2, "0")}:${String(s % 60).padStart(2, "0")}`; }),
      hovertemplate: "%{text}<br>ΔE %{x:.3f} m<br>ΔN %{y:.3f} m<extra></extra>",
      line: { width: 1, color: css("--line") }, marker: { size: 5, color: P.t, colorscale: "Viridis", showscale: true, colorbar: { title: { text: lang === "tr" ? "saat" : "h" }, thickness: 12 } }, name: K.kfplot ? t("kin.smooth") : "track" },
      { type: "scatter", mode: "markers", x: [P.E[S.j_epoch]], y: [P.N[S.j_epoch]], marker: { size: 13, color: css("--accent"), line: { color: css("--ink"), width: 1.5 } }, name: R.epoch.used }];
    const trLay = baseLayout({ showlegend: !!K.kfplot, xaxis: { ...baseLayout().xaxis, title: { text: "ΔE [m]" } }, yaxis: { ...baseLayout().yaxis, title: { text: "ΔN [m]" }, scaleanchor: "x", scaleratio: 1 } });
    const neuTr = ["N", "E", "U"].map((k, i) => ({ type: "scatter", mode: "lines", name: "Δ" + k, x: P.t, y: P[k], line: { width: 1.4, color: [css("--l1"), css("--l2"), PALETTE[2]][i] } }));
    const sigTr = ["sN", "sE", "sU"].map((k, i) => ({ type: "scatter", mode: "lines", name: "σ" + k.slice(1), x: P.t, y: P[k], line: { width: 1.4, color: [css("--l1"), css("--l2"), PALETTE[2]][i] } }));
    const nTr = [{ type: "scatter", mode: "lines", name: lang === "tr" ? "uydu" : "sats", x: P.t, y: P.n, line: { shape: "hv", width: 1.6, color: css("--l1") } },
      { type: "scatter", mode: "lines", name: "PDOP", x: P.t, y: P.pdop, yaxis: "y2", line: { width: 1.2, color: css("--l2") } }];
    const nLay = baseLayout({ xaxis: xa, yaxis: { ...baseLayout().yaxis, title: { text: lang === "tr" ? "uydu" : "satellites" } },
      yaxis2: { overlaying: "y", side: "right", title: { text: "PDOP" }, showgrid: false, linecolor: css("--line") } });
    const cols = [{ k: "t", l: "col.time" }, { k: "X", l: "col.X", d: 4 }, { k: "Y", l: "col.Y", d: 4 }, { k: "Z", l: "col.Z", d: 4 },
      { k: "N", l: "col.N", d: 3, sg: 1 }, { k: "E", l: "col.E", d: 3, sg: 1 }, { k: "U", l: "col.U", d: 3, sg: 1 },
      { k: "sN", l: "col.sN", d: 1 }, { k: "sE", l: "col.sE", d: 1 }, { k: "sU", l: "col.sU", d: 1 },
      { k: "n", l: "col.nsat", d: 0 }, { k: "pdop", l: "col.pdop", d: 2 }, { k: "v", l: "col.v", d: 2 }];
    let kfH = "";
    const F = S.kf, FP = K.kfplot;
    if (F && FP) {
      cols.push({ k: "fN", l: "col.fN", d: 3, sg: 1 }, { k: "fE", l: "col.fE", d: 3, sg: 1 }, { k: "fU", l: "col.fU", d: 3, sg: 1 });
      track.splice(1, 0, { type: "scatter", mode: "lines", x: FP.fE, y: FP.fN, line: { width: 1, color: css("--l2"), dash: "dot" }, name: t("kin.fwd"), hoverinfo: "skip" });
      const finfo = kv([
        [t("kin.kf.model"), `${esc(t("set.kf.dyn." + F.dyn))}, ${num(F.q_h, 3)} / ${num(F.q_v, 3)}`],
        [t("kin.kf.smoothed"), esc(t(F.smoothed ? "col.yes" : "col.no"))],
        [t("kin.kf.state"), F.n_state],
        [t("kin.kf.rej"), `${F.n_reject} / ${F.n_reset}`],
        [t("kin.kf.conv"), F.conv_min !== null ? num(F.conv_min, 1) + " " + esc(t("kin.kf.min")) : "–"],
        [t("kin.kf.diff"), (F.fwd_minus_smooth_rms_mm || []).map((x) => num(x, 1)).join(" / ") + " mm"],
        [t("kin.kf.diffLate"), (F.fwd_minus_smooth_rms_late_mm || []).map((x) => num(x, 1)).join(" / ") + " mm"],
      ]);
      const sigTr2 = [{ type: "scatter", mode: "lines", name: t("kin.fwd"), x: P.t, y: FP.fsH, line: { width: 1.6, color: css("--l2") } },
        { type: "scatter", mode: "lines", name: t("kin.smooth"), x: P.t, y: FP.sH, line: { width: 1.6, color: css("--l1") } }];
      const dTr = ["N", "E", "U"].map((k, i) => ({ type: "scatter", mode: "lines", name: "Δ" + k, x: P.t, y: FP["d" + k], line: { width: 1.3, color: [css("--l1"), css("--l2"), PALETTE[2]][i] } }));
      kfH = `<h4>${esc(t("kin.kf.title"))}</h4><p class="lead2">${esc(t("kin.kf.lead"))}</p>${finfo}
        <div class="grid2"><div>${chart("kfsig", t("kin.chart.kfsig"), sigTr2, baseLayout({ xaxis: xa, yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "mm" } } }))}</div>
        <div>${chart("kfdiff", t("kin.chart.kfdiff"), dTr, baseLayout({ xaxis: xa, yaxis: { ...baseLayout().yaxis, title: { text: "m" } } }))}</div></div>`;
    }
    return `<p class="lead2">${esc(t(F ? "l3.lead.kalman" : "kin.lead"))}</p>${info}
      ${chart("kintrack", t("kin.chart.track"), track, trLay, false).replace('class="chart"', 'class="chart square"')}
      ${chart("kinneu", t("kin.chart.neu"), neuTr, baseLayout({ xaxis: xa, shapes: [epochLine(pl.epoch_h)] }))}
      <div class="grid2"><div>${chart("kinsig", t("kin.chart.sig"), sigTr, baseLayout({ xaxis: xa }))}</div>
      <div>${chart("kinn", t("kin.chart.nsat"), nTr, nLay)}</div></div>
      ${kfH}
      ${table("kin", t("kin.table"), cols, K.rows, { maxRows: 400 })}
      <p class="note">${esc(t("kin.tableNote"))}</p>`;
  }

  function pSpp(R) {
    const S = R.spp, pl = R.plots;
    const ep = S.epoch;
    const info = kv([[t("spp.n"), `${S.n} (${esc(t("spp.source." + S.source))})`],
      [t("spp.epoch"), ep ? ep.X.map((v) => coord(v, 3)).join("; ") + " m" : "–"],
      [t("sum.clock"), ep ? num(ep.clock_ns, 1) + " ns" : "–"],
      [t("spp.mean"), S.mean ? S.mean.map((v) => coord(v, 3)).join("; ") + " m" : "–"],
      [t("spp.rms"), num(S.rms3d, 2) + " m"]]);
    let ch = "";
    if (pl.spp) {
      const tr = ["N", "E", "U"].map((k, i) => ({ type: "scatter", mode: "lines+markers", name: "Δ" + k, x: pl.spp.t, y: pl.spp[k], marker: { size: 3 }, line: { width: 1, color: [css("--l1"), css("--l2"), PALETTE[2]][i] } }));
      ch = chart("spp", t("spp.chart"), tr, baseLayout({ xaxis: { ...baseLayout().xaxis, title: { text: tLab() } }, yaxis: { ...baseLayout().yaxis, ...(R.kin ? {} : { range: [-10, 10] }), title: { text: "m" } }, shapes: [epochLine(pl.epoch_h)] }), true);
    }
    return `<p class="lead2">${esc(t("spp.lead"))}</p>${info}${ch}`;
  }

  // GitHub Pages'te yayinlaniyorsa kaynak kodu baglantisini repoya yonlendir
  (() => {
    const h = location.hostname, seg = location.pathname.split("/").filter(Boolean);
    if (h.endsWith(".github.io")) {
      const user = h.split(".")[0];
      $("#src-link").href = `https://github.com/${user}/${seg[0] && !seg[0].includes(".") ? seg[0] : user + ".github.io"}`;
    }
  })();

  // ------------------------------------------------------------------ baslat
  applyLang();
  drawHeroSky(true);
})();
