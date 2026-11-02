/* Hesaplari ayri bir is parcaciginda (Web Worker) Pyodide ile calistirir. */
const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
importScripts(PYODIDE_URL + "pyodide.js");

let py = null;
let ready = null;

async function boot() {
  post({ type: "boot", stage: "pyodide" });
  py = await loadPyodide({ indexURL: PYODIDE_URL });
  post({ type: "boot", stage: "numpy" });
  await py.loadPackage("numpy");
  post({ type: "boot", stage: "code" });
  const man = await (await fetch("py/manifest.json", { cache: "no-cache" })).json();
  py.FS.mkdirTree("/app/gnssfaz");
  py.FS.mkdirTree("/app/data");
  for (const f of man.files) {
    const txt = await (await fetch("py/" + f, { cache: "no-cache" })).text();
    py.FS.writeFile("/app/" + f, txt);
  }
  for (const name of ["igs08.atx.gz", "igs20.atx.gz"]) {
    const atx = await (await fetch("data/" + name)).arrayBuffer();
    py.FS.writeFile("/app/data/" + name, new Uint8Array(atx));
  }
  py.runPython("import sys; sys.path.insert(0, '/app'); import web_entry");
  post({ type: "ready" });
}

function post(msg) { self.postMessage(msg); }

self.onmessage = async (ev) => {
  const msg = ev.data;
  try {
    if (!ready) ready = boot();
    await ready;
    if (msg.type === "init") return;
    if (msg.type === "run") {
      const progress = (step, frac) => post({ type: "progress", step, frac });
      const mod = py.globals.get("web_entry");
      const run = mod.run_json;
      const atxPaths = JSON.stringify({ igs08: "/app/data/igs08.atx.gz", igs20: "/app/data/igs20.atx.gz" });
      const out = run(JSON.stringify(msg.files), JSON.stringify(msg.params), progress, atxPaths);
      run.destroy();
      mod.destroy();
      post({ type: "result", payload: out });
    }
  } catch (err) {
    post({ type: "result", payload: JSON.stringify({ ok: false, code: "worker", message: String(err && err.message || err) }) });
  }
};
