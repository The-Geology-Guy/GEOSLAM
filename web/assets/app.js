/* GEOSLAM web app: UI, data fetching, map. The SLAM maths runs in Python (worker.js). */

const COMCAT = "https://earthquake.usgs.gov/fdsnws/event/1/query";
const THREEDEP = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage";
const QFAULTS = "https://earthquake.usgs.gov/arcgis/rest/services/haz/Qfaults/MapServer/21/query";
const MAX_CELLS = 800;

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

/* ------------------------------------------------------------------ worker */
const worker = new Worker("assets/worker.js");
let seq = 0;
const pending = new Map();
worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === "progress") { busyText(m.data); return; }
  const p = pending.get(m.id);
  if (!p) return;
  pending.delete(m.id);
  m.ok ? p.resolve(m.result) : p.reject(Object.assign(new Error(m.error), { detail: m.detail }));
};
function rpc(op, args = [], transfer = []) {
  return new Promise((resolve, reject) => {
    const id = ++seq;
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, op, args }, transfer);
  });
}
const pycall = (fn, ...args) => rpc("call", [fn, ...args]);

const pyReady = (async () => {
  const man = await fetch("py/manifest.json").then((r) => r.json());
  return rpc("init", ["../py/" + man.wheel]);
})();
pyReady.catch((e) => setStatus("The calculator could not start: " + e.message, "err"));

/* ------------------------------------------------------------------- state */
const S = {
  kind: null,          // "comcat" | "manual"
  caseObj: null,
  mech: 0,
  su: null, du: null,
  plan: null,
  demKey: null, demInfo: null,
  upload: null,        // {name, kind, data}
  results: null,
  layers: { NP1: true, NP2: true },
  lightPlane: 0,
  faultsKey: null, faults: null,
  quakesKey: null, quakes: null,
};

/* ---------------------------------------------------------------------- UI */
function setStatus(msg, kind = "") {
  const el = $("#status");
  el.textContent = msg || "";
  el.className = "status " + kind;
}
function busy(on, text = "") { $("#busy").hidden = !on; busyText(text); }
function busyText(t) { if (t) $("#busy-text").textContent = t; }

$$(".tab").forEach((t) => t.addEventListener("click", () => {
  $$(".tab").forEach((x) => { x.classList.toggle("is-on", x === t); x.setAttribute("aria-selected", x === t); });
  $$(".tabpane").forEach((p) => { p.hidden = p.dataset.pane !== t.dataset.tab; });
}));
$$("[data-open]").forEach((b) => b.addEventListener("click", () => $("#dlg-" + b.dataset.open).showModal()));
$$(".example").forEach((b) => b.addEventListener("click", () => { $("#event-id").value = b.dataset.id; loadEvent(b.dataset.id); }));
$("#load-event").addEventListener("click", () => loadEvent($("#event-id").value.trim()));
$("#event-id").addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); loadEvent(e.target.value.trim()); } });
$("#do-search").addEventListener("click", searchArea);
$("#use-manual").addEventListener("click", useManual);
$$('input[name="dem-src"]').forEach((r) => r.addEventListener("change", () => {
  $("#upload-box").hidden = demSource() !== "upload";
  S.demKey = null;
  refreshPlan();
}));
$("#dem-file").addEventListener("change", onFile);
$("#flow").addEventListener("submit", (e) => { e.preventDefault(); runAll(); });

for (const [id, out, fmt] of [["su", "su-out", (v) => `±${v}°`], ["du", "du-out", (v) => `±${v}°`], ["mult", "mult-out", (v) => `×${(+v).toFixed(1)}`]]) {
  const el = $("#" + id);
  const show = () => { $("#" + out).textContent = fmt(el.value); };
  el.addEventListener("input", () => { show(); if (id !== "mult") { S[id] = +el.value; } liveRerun(); });
  show();
}
$$('input[name="method"], input[name="datum"]').forEach((r) => r.addEventListener("change", liveRerun));
["#a-eh", "#a-ez"].forEach((s) => $(s).addEventListener("change", liveRerun));
["#half-km", "#cell-m"].forEach((s) => $(s).addEventListener("change", () => { S.demKey = null; refreshPlan(); }));
$("#show-bounds").addEventListener("change", () => setVis("bounds", $("#show-bounds").checked));
$("#show-hs").addEventListener("change", () => setVis("hillshade", $("#show-hs").checked));

const today = new Date();
$("#s-end").value = today.toISOString().slice(0, 10);
$("#s-start").value = new Date(today.getTime() - 15 * 365.25 * 864e5).toISOString().slice(0, 10);

const demSource = () => $('input[name="dem-src"]:checked').value;
const method = () => $('input[name="method"]:checked').value;
const datum = () => $('input[name="datum"]:checked').value;
const num = (s) => { const v = parseFloat($(s).value); return Number.isFinite(v) ? v : null; };
const fmtDeg = (v) => (v == null ? "?" : Math.round(v));
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

/* -------------------------------------------------------------------- map */
const BASEMAPS = {
  topo: { tiles: ["https://basemap.nationalmap.gov/arcgis/rest/services/USGSTopo/MapServer/tile/{z}/{y}/{x}"], max: 16, attr: "USGS The National Map" },
  relief: { tiles: ["https://basemap.nationalmap.gov/arcgis/rest/services/USGSShadedReliefOnly/MapServer/tile/{z}/{y}/{x}"], max: 13, attr: "USGS The National Map" },
  imagery: { tiles: ["https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}"], max: 16, attr: "USGS The National Map" },
  osm: { tiles: ["https://tile.openstreetmap.org/{z}/{x}/{y}.png"], max: 19, attr: "© OpenStreetMap contributors" },
};
const style = { version: 8, sources: {}, layers: [{ id: "bg", type: "background", paint: { "background-color": "#dfe4e8" } }] };
for (const [k, b] of Object.entries(BASEMAPS)) {
  style.sources["bm-" + k] = { type: "raster", tiles: b.tiles, tileSize: 256, maxzoom: b.max, attribution: b.attr };
  style.layers.push({ id: "bm-" + k, type: "raster", source: "bm-" + k, layout: { visibility: k === "topo" ? "visible" : "none" } });
}
const map = new maplibregl.Map({
  container: "map", style, center: [-98.5, 39], zoom: 3.4, preserveDrawingBuffer: true,
  attributionControl: { compact: true },
});
map.addControl(new maplibregl.NavigationControl({ visualizePitch: false }), "top-right");
map.addControl(new maplibregl.ScaleControl({ unit: "metric" }), "bottom-right");
$("#basemap").addEventListener("change", (e) => {
  for (const k of Object.keys(BASEMAPS)) map.setLayoutProperty("bm-" + k, "visibility", k === e.target.value ? "visible" : "none");
});

const mapReady = new Promise((res) => map.on("load", res));
const EMPTY = { type: "FeatureCollection", features: [] };

mapReady.then(() => {
  const c = { np1: css("--np1"), np2: css("--np2"), fault: css("--fault"), quake: css("--quake"), ink: "#1e2733" };
  for (const id of ["swaths", "traces", "bounds", "faults", "quakes", "epi", "search"]) map.addSource(id, { type: "geojson", data: EMPTY });
  const planeColor = ["match", ["get", "plane"], "NP2", c.np2, c.np1];
  map.addLayer({ id: "swath-fill", type: "fill", source: "swaths", paint: { "fill-color": planeColor, "fill-opacity": 0.22 } });
  map.addLayer({ id: "swath-line", type: "line", source: "swaths", paint: { "line-color": planeColor, "line-width": 2.2 } });
  map.addLayer({ id: "bounds", type: "line", source: "bounds", paint: { "line-color": planeColor, "line-width": 0.9, "line-opacity": 0.75 } });
  map.addLayer({ id: "faults", type: "line", source: "faults", paint: { "line-color": c.fault, "line-width": 1.6 } });
  map.addLayer({ id: "traces", type: "line", source: "traces", paint: { "line-color": c.ink, "line-width": 2 } });
  map.addLayer({ id: "quakes", type: "circle", source: "quakes", paint: {
    "circle-radius": ["interpolate", ["linear"], ["get", "mag"], 1, 2, 3, 4, 5, 9, 7, 16],
    "circle-color": "rgba(0,0,0,0)", "circle-stroke-color": c.quake, "circle-stroke-width": 1.2 } });
  map.addLayer({ id: "search", type: "circle", source: "search", paint: {
    "circle-radius": ["interpolate", ["linear"], ["get", "mag"], 3, 5, 7, 13],
    "circle-color": c.np1, "circle-opacity": 0.5, "circle-stroke-color": "#fff", "circle-stroke-width": 1 } });
  map.addLayer({ id: "epi", type: "symbol", source: "epi", layout: {
    "icon-image": "beachball", "icon-size": 0.42, "icon-allow-overlap": true, "icon-offset": [62, -62] } });
  map.addLayer({ id: "epi-star", type: "circle", source: "epi", paint: {
    "circle-radius": 6, "circle-color": "#ffd400", "circle-stroke-color": "#1e2733", "circle-stroke-width": 1.5 } });
  map.on("click", "search", (e) => { const id = e.features[0].properties.id; $("#event-id").value = id; loadEvent(id); });
  map.on("mouseenter", "search", () => { map.getCanvas().style.cursor = "pointer"; });
  map.on("mouseleave", "search", () => { map.getCanvas().style.cursor = ""; });
  map.on("click", "faults", (e) => {
    const p = e.features[0].properties;
    new maplibregl.Popup().setLngLat(e.lngLat).setHTML(`<strong>${esc(p.fault_name)}</strong><br>${esc(p.age || "")}${p.slip_sense ? ", " + esc(p.slip_sense) : ""}`).addTo(map);
  });
  readHash();
});

function setData(id, data) { const s = map.getSource(id); if (s) s.setData(data); }
function setVis(id, on) {
  if (id === "hillshade") { if (map.getLayer("hillshade")) map.setLayoutProperty("hillshade", "visibility", on ? "visible" : "none"); return; }
  if (map.getLayer(id)) map.setLayoutProperty(id, "visibility", on ? "visible" : "none");
}

/* ---------------------------------------------------------------- loading */
async function fetchJSON(url, what) {
  const r = await fetch(url);
  if (r.status === 404 || r.status === 204) throw new Error(`${what} not found.`);
  if (!r.ok) throw new Error(`${what} request failed (HTTP ${r.status}).`);
  return r.json();
}

async function loadEvent(id) {
  if (!id) { setStatus("Enter a ComCat event ID, such as nc72282711.", "err"); return; }
  setStatus(`Loading ${id} from the USGS catalog…`);
  try {
    const detail = await fetchJSON(`${COMCAT}?eventid=${encodeURIComponent(id)}&format=geojson`, `Event ${id}`);
    const products = detail.properties?.products || {};
    const urls = [];
    for (const t of ["focal-mechanism", "moment-tensor"]) for (const p of products[t] || []) {
      const u = p.contents?.["quakeml.xml"]?.url; if (u) urls.push(u);
    }
    const qml = {};
    await Promise.all(urls.map(async (u) => { try { qml[u] = await (await fetch(u)).text(); } catch { /* skip */ } }));
    await pyReady;
    const c = JSON.parse(await pycall("case_from_comcat", JSON.stringify(detail), JSON.stringify(qml)));
    if (!c.mechanisms.length) {
      prefillManual(c.event);
      setStatus(`${id} has no focal mechanism in ComCat. Its location is copied into “Enter values” so you can add a mechanism from another source.`, "err");
      return;
    }
    setCase("comcat", c);
    setStatus("");
  } catch (e) {
    setStatus(e.message, "err");
  }
}

function prefillManual(ev) {
  $$(".tab").find((t) => t.dataset.tab === "manual").click();
  $("#m-lat").value = ev.lat.toFixed(5); $("#m-lon").value = ev.lon.toFixed(5);
  $("#m-depth").value = ev.depth_km; if (ev.mag != null) $("#m-mag").value = ev.mag;
}

async function useManual() {
  const v = Object.fromEntries(["lat", "lon", "depth", "strike", "dip", "rake", "su", "du", "mag", "eh1", "eh1az", "eh2", "ez"].map((k) => [k, num("#m-" + k)]));
  const missing = ["lat", "lon", "depth", "strike", "dip", "rake"].filter((k) => v[k] == null);
  if (missing.length) { setStatus("Fill in latitude, longitude, depth, strike, dip and rake.", "err"); return; }
  if (v.dip <= 0 || v.dip > 90) { setStatus("Dip must be between 0 and 90 degrees.", "err"); return; }
  await pyReady;
  const np1 = { strike: v.strike, dip: v.dip, rake: v.rake, su: v.su || 0, du: v.du || 0, ru: 0 };
  const np2 = JSON.parse(await pycall("auxiliary", JSON.stringify(np1)));
  const unc = v.eh1 != null && v.ez != null ? { eh1: v.eh1, eh1_az: v.eh1az || 0, eh2: v.eh2 ?? v.eh1, ez: v.ez } : null;
  setCase("manual", {
    event: { id: "entered", lon: v.lon, lat: v.lat, depth_km: v.depth, mag: v.mag, mag_type: "", time: "", title: "Entered hypocentre" },
    uncertainty: unc,
    mechanisms: [{ label: "Entered mechanism", source: "user", product: "manual", method: "", has_uncertainty: !!(v.su || v.du), np1, np2 }],
  });
}

async function setCase(kind, c, keep = {}) {
  S.kind = kind; S.caseObj = c; S.mech = keep.mech ?? 0; S.results = null; S.demKey = null;
  S.faultsKey = S.quakesKey = null;
  const ev = c.event;
  const when = ev.time ? new Date(ev.time).toISOString().replace("T", " ").slice(0, 16) + " UTC" : "";
  $("#event-card").hidden = false;
  $("#event-card").innerHTML = `<div class="title">${esc(ev.title || ev.id)}</div>
    <div class="meta">${ev.mag != null ? `M${(+ev.mag).toFixed(1)}${ev.mag_type ? " " + esc(ev.mag_type) : ""}, ` : ""}${when ? when + ", " : ""}depth ${(+ev.depth_km).toFixed(1)} km, ${ev.lat.toFixed(4)}°, ${ev.lon.toFixed(4)}°${ev.url ? `, <a href="${esc(ev.url)}" target="_blank" rel="noopener">USGS page</a>` : ""}</div>`;
  $("#step-event").classList.add("is-done");
  renderMechs(keep);
  await renderUncertainty();
  ["#step-mech", "#step-unc", "#step-dem"].forEach((s) => { $(s).hidden = false; });
  $("#run").disabled = false;
  setData("search", EMPTY);
  setData("epi", { type: "FeatureCollection", features: [{ type: "Feature", geometry: { type: "Point", coordinates: [ev.lon, ev.lat] }, properties: {} }] });
  await refreshPlan();
  if (S.plan) map.fitBounds([[S.plan.lonlat_bounds[0], S.plan.lonlat_bounds[1]], [S.plan.lonlat_bounds[2], S.plan.lonlat_bounds[3]]], { padding: 30, duration: 600 });
  updateBeachball();
}

function renderMechs(keep = {}) {
  const box = $("#mech-list");
  box.innerHTML = "";
  S.caseObj.mechanisms.forEach((m, i) => {
    const p1 = m.np1, p2 = m.np2;
    const unc = m.has_uncertainty
      ? `<span class="tag">Uncertainty reported: strike ±${fmtDeg(p1.su)}°, dip ±${fmtDeg(p1.du)}°</span>`
      : `<span class="tag none">No orientation uncertainty reported</span>`;
    const lab = document.createElement("label");
    lab.className = "mech";
    lab.innerHTML = `<input type="radio" name="mech" value="${i}" ${i === S.mech ? "checked" : ""}>
      <div><div><strong>${esc(m.label)}</strong></div>
      <div class="planes">NP1 ${fmtDeg(p1.strike)}/${fmtDeg(p1.dip)}/${fmtDeg(p1.rake)}, NP2 ${fmtDeg(p2.strike)}/${fmtDeg(p2.dip)}/${fmtDeg(p2.rake)} <span class="muted">(strike/dip/rake)</span></div>
      ${unc}</div>`;
    lab.querySelector("input").addEventListener("change", () => { S.mech = i; applyMechUncertainty(); updateBeachball(); liveRerun(); });
    box.appendChild(lab);
  });
  applyMechUncertainty(keep);
}

function applyMechUncertainty(keep = {}) {
  const m = S.caseObj.mechanisms[S.mech];
  const assumed = !m.has_uncertainty;
  S.su = keep.su ?? (assumed ? 10 : Math.round(m.np1.su));
  S.du = keep.du ?? (assumed ? 10 : Math.round(m.np1.du));
  $("#su").value = S.su; $("#du").value = S.du;
  $("#su").dispatchEvent(new Event("input")); $("#du").dispatchEvent(new Event("input"));
  $("#mech-note").textContent = assumed
    ? "This solution has no published orientation uncertainty, so ±10° is assumed (the convention used by Rasaka, 2016). Adjust the sliders to test other values."
    : "Starting values are the catalog's reported uncertainties. Move the sliders to see how the swath responds.";
}

async function renderUncertainty() {
  const u = S.caseObj.uncertainty;
  $("#unc-assume").hidden = !!u;
  if (!u) {
    $("#unc-text").textContent = "The catalog reports no location uncertainty for this hypocentre. Enter assumed values; Cronin's students used about ±2 km for poorly located historical events.";
    return;
  }
  const d = JSON.parse(await pycall("describe_uncertainty", JSON.stringify(u)));
  const conf = d.confidence ? `${Math.round(d.confidence * 100)}% confidence` : "confidence level not stated";
  const h = d.horizontal_m[0] === d.horizontal_m[1] ? `±${d.horizontal_m[0]} m` : `±${d.horizontal_m[0]}–${d.horizontal_m[1]} m`;
  $("#unc-text").innerHTML = `Horizontal ${h}, vertical ±${d.vertical_m} m <span class="muted">(${conf}; ${esc(prettySource(d.source))})</span>.`;
}
const prettySource = (s) => s.startsWith("comcat:error-ellipse") ? "3-D error ellipsoid from ComCat"
  : s.startsWith("comcat:") ? "ComCat horizontal and vertical error" : s === "entered" ? "entered values" : s;

async function refreshPlan() {
  if (!S.caseObj) return;
  await pyReady;
  const half = num("#half-km"), cell = num("#cell-m");
  S.plan = JSON.parse(await pycall("dem_plan", JSON.stringify(caseForRun()), half, cell, MAX_CELLS));
  const km = ((S.plan.bbox[2] - S.plan.bbox[0]) / 1000).toFixed(0);
  const zone = `UTM zone ${S.plan.zone}${S.plan.north ? "N" : "S"}`;
  $("#dem-plan").textContent = demSource() === "3dep"
    ? `Will fetch a ${km} × ${km} km area at ${S.plan.cellsize} m resolution (${S.plan.size[0]} × ${S.plan.size[1]} cells), projected to ${zone}.`
    : `Analysis area ${km} × ${km} km around the epicentre, in ${zone}.`;
}

/* case with slider overrides applied (sent to Python) */
function caseForRun() {
  const c = structuredClone(S.caseObj);
  if (S.su != null) for (const m of c.mechanisms) for (const k of ["np1", "np2"]) { m[k].su = S.su; m[k].du = S.du; }
  return c;
}

/* ------------------------------------------------------------------ search */
async function searchArea() {
  const b = map.getBounds();
  const params = new URLSearchParams({
    format: "geojson", starttime: $("#s-start").value, endtime: $("#s-end").value,
    minmagnitude: $("#s-minmag").value || "3", minlatitude: b.getSouth().toFixed(4), maxlatitude: b.getNorth().toFixed(4),
    minlongitude: b.getWest().toFixed(4), maxlongitude: b.getEast().toFixed(4), orderby: "magnitude", limit: "200",
  });
  setStatus("Searching the USGS catalog…");
  try {
    const sets = await Promise.all(["focal-mechanism", "moment-tensor"].map((pt) =>
      fetchJSON(`${COMCAT}?${params}&producttype=${pt}`, "Search").catch(() => EMPTY)));
    const byId = new Map();
    for (const s of sets) for (const f of s.features || []) byId.set(f.id, f);
    const feats = [...byId.values()].sort((a, b) => b.properties.mag - a.properties.mag).slice(0, 150);
    const list = $("#search-results");
    list.innerHTML = "";
    for (const f of feats) {
      const li = document.createElement("li");
      const d = new Date(f.properties.time).toISOString().slice(0, 10);
      li.innerHTML = `<button type="button"><strong>M${(+f.properties.mag).toFixed(1)}</strong> ${esc(f.properties.place || f.id)}<div class="meta">${d}, depth ${(+f.geometry.coordinates[2]).toFixed(1)} km, ${esc(f.id)}</div></button>`;
      li.querySelector("button").addEventListener("click", () => { $("#event-id").value = f.id; loadEvent(f.id); });
      list.appendChild(li);
    }
    setData("search", { type: "FeatureCollection", features: feats.map((f) => ({ type: "Feature", geometry: { type: "Point", coordinates: f.geometry.coordinates.slice(0, 2) }, properties: { id: f.id, mag: f.properties.mag } })) });
    setStatus(feats.length ? `${feats.length} earthquakes with focal mechanisms. Pick one from the list or the map.` : "No earthquakes with focal mechanisms here. Zoom out, widen the dates or lower the magnitude.", feats.length ? "" : "err");
  } catch (e) { setStatus(e.message, "err"); }
}

/* -------------------------------------------------------------------- DEM */
async function onFile(e) {
  const f = e.target.files[0];
  S.upload = null; S.demKey = null;
  if (!f) return;
  const isAscii = /\.(asc|txt|dat)$/i.test(f.name);
  $("#ascii-zone").hidden = !isAscii;
  if (isAscii && S.plan && !$("#ascii-zone-n").value) $("#ascii-zone-n").value = S.plan.zone;
  S.upload = { name: f.name, kind: isAscii ? "ascii" : "tiff", data: isAscii ? await f.text() : await f.arrayBuffer() };
  setStatus(`${f.name} ready (${(f.size / 1e6).toFixed(1)} MB).`, "ok");
}

async function ensureDem() {
  const src = demSource();
  const key = src === "3dep" ? `3dep:${S.plan.epsg}:${S.plan.bbox}:${S.plan.cellsize}` : `up:${S.upload?.name}:${S.plan.bbox}:${$("#ascii-zone-n").value}${$("#ascii-hemi").value}`;
  if (S.demKey === key) return;
  if (src === "3dep") {
    busy(true, "Fetching terrain from USGS 3DEP");
    const [n, m] = S.plan.size;
    const p = new URLSearchParams({ bbox: S.plan.bbox.join(","), bboxSR: S.plan.epsg, imageSR: S.plan.epsg, size: `${n},${m}`,
      format: "tiff", pixelType: "F32", interpolation: "RSP_BilinearInterpolation", noDataInterpretation: "esriNoDataMatchAny", f: "image" });
    const r = await fetch(`${THREEDEP}?${p}`);
    if (!r.ok) throw new Error(`USGS 3DEP did not return terrain (HTTP ${r.status}). Try again, or upload your own DEM.`);
    const { data, w, h, nodata } = await readTiff(await r.arrayBuffer());
    let valid = 0;
    for (let i = 0; i < data.length; i++) { const v = data[i]; if (v === nodata || v < -1e4 || v > 1e5 || Number.isNaN(v)) data[i] = NaN; else valid++; }
    if (valid < 0.5 * data.length) throw new Error("USGS 3DEP has little or no elevation data here (it covers the United States). Choose “My own DEM” and upload a GeoTIFF.");
    const cs = (S.plan.bbox[2] - S.plan.bbox[0]) / w;
    busyText("Preparing terrain");
    S.demInfo = JSON.parse(await rpc("setDemUtm", [data, h, w, S.plan.bbox[0], S.plan.bbox[3], cs, S.plan.epsg, "USGS 3DEP"], [data.buffer]));
  } else {
    if (!S.upload) throw new Error("Choose a DEM file to upload, or switch to USGS 3DEP.");
    busy(true, "Reading your DEM");
    if (S.upload.kind === "ascii") {
      const z = num("#ascii-zone-n");
      if (!z) throw new Error("Enter the UTM zone your ASCII grid is in.");
      await rpc("setDemAscii", [S.upload.data, z, $("#ascii-hemi").value === "N"]);
    } else {
      const { data, w, h, nodata, origin, res, epsg, geographic } = await readTiff(S.upload.data.slice(0));
      for (let i = 0; i < data.length; i++) if (data[i] === nodata) data[i] = NaN;
      if (geographic) {
        await rpc("setDemGeo", [data, h, w, origin[0], origin[1], res[0], -res[1], S.caseObj.event.lon, S.caseObj.event.lat, num("#cell-m"), JSON.stringify(S.plan.bbox), S.upload.name], [data.buffer]);
      } else if (epsg && ((epsg >= 32601 && epsg <= 32760) || (epsg >= 26901 && epsg <= 26923))) {
        if (Math.abs(res[0] + res[1]) > 1e-6 * Math.abs(res[0])) throw new Error("The DEM's cells are not square. Resample it to square cells first.");
        await rpc("setDemUtm", [data, h, w, origin[0], origin[1], res[0], epsg, S.upload.name], [data.buffer]);
      } else {
        throw new Error(`This GeoTIFF uses EPSG:${epsg || "unknown"}. Use a lat/lon (WGS84/NAD83) or UTM GeoTIFF.`);
      }
    }
    await pycall("crop_dem", JSON.stringify(S.plan.bbox));
    S.demInfo = JSON.parse(await pycall("thin_dem", MAX_CELLS));
    if (S.demInfo.valid_fraction < 0.3) throw new Error("Your DEM barely overlaps the area around this earthquake. Check its location and coordinate system.");
  }
  S.demKey = key;
  await drawHillshade();
}

async function readTiff(buf) {
  const tiff = await GeoTIFF.fromArrayBuffer(buf);
  const img = await tiff.getImage();
  const w = img.getWidth(), h = img.getHeight();
  const ras = await img.readRasters({ samples: [0], interleave: true });
  const data = new Float32Array(ras);  // own copy, safe to transfer
  const nd = img.getGDALNoData();
  const keys = img.getGeoKeys() || {};
  const epsg = keys.ProjectedCSTypeGeoKey || keys.GeographicTypeGeoKey || null;
  const geographic = !keys.ProjectedCSTypeGeoKey && (keys.GTModelTypeGeoKey === 2 || !!keys.GeographicTypeGeoKey);
  let origin = null, res = null;
  try { origin = img.getOrigin(); res = img.getResolution(); } catch { /* 3DEP path doesn't need it */ }
  return { data, w, h, nodata: nd == null ? NaN : nd, origin, res, epsg, geographic };
}

async function drawHillshade() {
  const r = S.results?.planes?.[S.lightPlane];
  const az = r ? r.light_azimuth : 315;
  const rgba = await rpc("hillshade", [az, 15]);
  const [h, w] = S.demInfo.shape;
  const cv = document.createElement("canvas");
  cv.width = w; cv.height = h;
  cv.getContext("2d").putImageData(new ImageData(rgba, w, h), 0, 0);
  const url = cv.toDataURL("image/png");
  const coords = S.demInfo.corners;
  if (map.getSource("hillshade")) map.getSource("hillshade").updateImage({ url, coordinates: coords });
  else {
    map.addSource("hillshade", { type: "image", url, coordinates: coords });
    map.addLayer({ id: "hillshade", type: "raster", source: "hillshade", paint: { "raster-opacity": 0.55, "raster-fade-duration": 0 } }, "swath-fill");
  }
  setVis("hillshade", $("#show-hs").checked);
}

/* ------------------------------------------------------------- context data */
async function loadContext() {
  const [w, s, e, n] = S.plan.lonlat_bounds;
  const fk = S.plan.lonlat_bounds.join();
  const tasks = [];
  if (S.faultsKey !== fk) tasks.push((async () => {
    const feats = [];
    for (let off = 0; off < 10000; off += 2000) {
      const p = new URLSearchParams({ geometry: `${w},${s},${e},${n}`, geometryType: "esriGeometryEnvelope", inSR: "4326",
        spatialRel: "esriSpatialRelIntersects", outFields: "fault_name,section_name,age,slip_sense,slip_rate", outSR: "4326",
        maxAllowableOffset: "0.0001", geometryPrecision: "6", f: "geojson", resultOffset: off, resultRecordCount: 2000 });
      const js = await fetchJSON(`${QFAULTS}?${p}`, "Fault database").catch(() => EMPTY);
      feats.push(...(js.features || []));
      if (!js.exceededTransferLimit && (js.features || []).length < 2000) break;
    }
    S.faults = { type: "FeatureCollection", features: feats }; S.faultsKey = fk;
    setData("faults", S.faults);
  })());
  const ev = S.caseObj.event;
  const qk = `${fk}:${$("#seis-days").value}:${$("#seis-mag").value}`;
  if (ev.time && S.quakesKey !== qk) tasks.push((async () => {
    const t0 = new Date(ev.time);
    const t1 = new Date(t0.getTime() + (+$("#seis-days").value || 0) * 864e5 + 6e4);
    const p = new URLSearchParams({ format: "geojson", starttime: t0.toISOString(), endtime: t1.toISOString(),
      minlatitude: s, maxlatitude: n, minlongitude: w, maxlongitude: e, minmagnitude: $("#seis-mag").value || "1.5", limit: "5000" });
    const js = await fetchJSON(`${COMCAT}?${p}`, "Seismicity").catch(() => EMPTY);
    S.quakes = { type: "FeatureCollection", features: (js.features || []).filter((f) => f.id !== ev.id).map((f) => ({ type: "Feature", geometry: { type: "Point", coordinates: f.geometry.coordinates.slice(0, 2) }, properties: { mag: f.properties.mag } })) };
    S.quakesKey = qk;
    setData("quakes", S.quakes);
  })());
  else if (!ev.time) setData("quakes", EMPTY);
  await Promise.all(tasks);
}

/* --------------------------------------------------------------------- run */
let running = false, rerunQueued = false;
async function runAll() {
  if (!S.caseObj) return;
  if (running) { rerunQueued = true; return; }
  running = true;
  $("#run").disabled = true;
  try {
    setStatus("");
    busy(true, "Starting the calculator (first run takes a few seconds)");
    await pyReady;
    await refreshPlan();
    const ctx = loadContext();
    await ensureDem();
    busyText("Projecting the nodal planes onto the terrain");
    const t0 = performance.now();
    const opts = {
      method: method(), depth_datum: datum(), multiplier: +$("#mult").value, chunk: 20000,
      assumed_location_m: S.caseObj.uncertainty ? null : [num("#a-eh") || 2000, num("#a-ez") || 2000],
      planes: ["np1", "np2"].map((k, i) => ({ label: "NP" + (i + 1), plane: caseForRun().mechanisms[S.mech][k] })),
    };
    S.results = JSON.parse(await pycall("run", JSON.stringify(caseForRun()), JSON.stringify(opts)));
    S.results.seconds = (performance.now() - t0) / 1000;
    S.results.opts = opts;
    drawResults();
    await drawHillshade();
    await ctx;
    renderResults();
    writeHash();
    setStatus(`Done in ${S.results.seconds.toFixed(1)} s.`, "ok");
  } catch (e) {
    console.error(e.detail || e);
    setStatus(e.message, "err");
  } finally {
    busy(false);
    running = false;
    $("#run").disabled = false;
    if (rerunQueued) { rerunQueued = false; runAll(); }
  }
}
let liveTimer = null;
function liveRerun() {
  if (!S.results) return;
  clearTimeout(liveTimer);
  liveTimer = setTimeout(runAll, 350);
}

function drawResults() {
  const sw = [], tr = [], bd = [];
  for (const p of S.results.planes) {
    if (!S.layers[p.label]) continue;
    const props = { plane: p.label, strike: p.strike, dip: p.dip, rake: p.rake, area_km2: p.area_km2 };
    if (p.polygons.length) sw.push({ type: "Feature", properties: props, geometry: { type: "MultiPolygon", coordinates: p.polygons } });
    for (const [name, lines] of Object.entries(p.traces)) {
      const f = { type: "Feature", properties: { plane: p.label, trace: name }, geometry: { type: "MultiLineString", coordinates: lines } };
      (name === "Mean" ? tr : bd).push(f);
    }
  }
  setData("swaths", { type: "FeatureCollection", features: sw });
  setData("traces", { type: "FeatureCollection", features: tr });
  setData("bounds", { type: "FeatureCollection", features: bd });
  setVis("bounds", $("#show-bounds").checked);
  $("#legend").hidden = false;
}

function renderResults() {
  const box = $("#plane-results");
  box.innerHTML = "";
  S.results.planes.forEach((p, i) => {
    const div = document.createElement("div");
    div.className = "plane-res " + p.label.toLowerCase();
    const rows = p.planes.filter((q) => q.name !== "MThin").map((q) =>
      `<tr><td>${q.name}</td><td>${q.dip.toFixed(1)}°</td><td>${q.dip_azimuth.toFixed(1)}°</td><td>${Math.round(q.half_width_m)} m</td></tr>`).join("");
    div.innerHTML = `<h3><span>Nodal plane ${i + 1}</span>
        <label class="toggle"><input type="checkbox" ${S.layers[p.label] ? "checked" : ""}> Show</label></h3>
      <div class="nums">Strike ${p.strike.toFixed(0)}°, dip ${p.dip.toFixed(0)}°, rake ${p.rake == null ? "?" : p.rake.toFixed(0)}°, uncertainty ±${p.su.toFixed(0)}° / ±${p.du.toFixed(0)}°</div>
      <div class="nums">Seismo-lineament area ${p.area_km2.toFixed(1)} km² within the ${((S.plan.bbox[2] - S.plan.bbox[0]) / 1000).toFixed(0)} km map</div>
      <label class="toggle"><input type="radio" name="light" value="${i}" ${S.lightPlane === i ? "checked" : ""}> Light the terrain across this plane</label>
      <details><summary>Bounding planes (grid north, convergence ${p.gamma.toFixed(2)}°)</summary>
        <table><tr><th>Plane</th><th>Dip</th><th>Dip dir.</th><th>Half-width</th></tr>${rows}</table></details>`;
    div.querySelector('input[type="checkbox"]').addEventListener("change", (e) => { S.layers[p.label] = e.target.checked; drawResults(); });
    div.querySelector('input[type="radio"]').addEventListener("change", () => { S.lightPlane = i; drawHillshade(); });
    box.appendChild(div);
  });
  $("#results").hidden = false;
}

/* --------------------------------------------------------------- beachball */
async function updateBeachball() {
  const m = S.caseObj?.mechanisms?.[S.mech];
  if (!m || m.np1.rake == null) return;
  await pyReady;
  const n = 161;
  const rgba = await rpc("beachball", [m.np1.strike, m.np1.dip, m.np1.rake, n]);
  if (map.hasImage("beachball")) map.removeImage("beachball");
  map.addImage("beachball", { width: n, height: n, data: new Uint8Array(rgba.buffer) });
}

/* ----------------------------------------------------------------- exports */
function download(name, blob) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
const baseName = () => `geoslam_${S.caseObj?.event?.id || "result"}`;
function resultsGeoJSON() {
  const feats = [];
  const meta = { event_id: S.caseObj.event.id, method: S.results.opts.method, multiplier: S.results.opts.multiplier, depth_datum: S.results.opts.depth_datum };
  for (const p of S.results.planes) {
    const props = { ...meta, plane: p.label, strike: p.strike, dip: p.dip, rake: p.rake, strike_unc: p.su, dip_unc: p.du, area_km2: p.area_km2 };
    feats.push({ type: "Feature", properties: { ...props, kind: "seismo-lineament" }, geometry: { type: "MultiPolygon", coordinates: p.polygons } });
    for (const [name, lines] of Object.entries(p.traces))
      feats.push({ type: "Feature", properties: { ...props, kind: name === "Mean" ? "mean_trace" : "bounding_plane_trace", trace: name }, geometry: { type: "MultiLineString", coordinates: lines } });
  }
  const ev = S.caseObj.event;
  feats.push({ type: "Feature", properties: { ...meta, kind: "epicenter", depth_km: ev.depth_km, mag: ev.mag }, geometry: { type: "Point", coordinates: [ev.lon, ev.lat] } });
  return { type: "FeatureCollection", features: feats };
}
$("#dl-geojson").addEventListener("click", () => download(baseName() + ".geojson", new Blob([JSON.stringify(resultsGeoJSON())], { type: "application/geo+json" })));
$("#dl-kml").addEventListener("click", () => {
  const coords = (ring) => ring.map((c) => `${c[0]},${c[1]},0`).join(" ");
  const pm = [];
  for (const f of resultsGeoJSON().features) {
    const p = f.properties, nm = `${p.plane || ""} ${p.kind}${p.trace && p.trace !== "Mean" ? " " + p.trace : ""}`.trim();
    const st = p.plane === "NP2" ? "np2" : "np1";
    let g = "";
    if (f.geometry.type === "MultiPolygon") g = `<MultiGeometry>${f.geometry.coordinates.map((poly) => `<Polygon><outerBoundaryIs><LinearRing><coordinates>${coords(poly[0])}</coordinates></LinearRing></outerBoundaryIs>${poly.slice(1).map((h) => `<innerBoundaryIs><LinearRing><coordinates>${coords(h)}</coordinates></LinearRing></innerBoundaryIs>`).join("")}</Polygon>`).join("")}</MultiGeometry>`;
    else if (f.geometry.type === "MultiLineString") g = `<MultiGeometry>${f.geometry.coordinates.map((l) => `<LineString><coordinates>${coords(l)}</coordinates></LineString>`).join("")}</MultiGeometry>`;
    else g = `<Point><coordinates>${f.geometry.coordinates.join(",")},0</coordinates></Point>`;
    pm.push(`<Placemark><name>${esc(nm)}</name><styleUrl>#${p.kind === "epicenter" ? "epi" : st}</styleUrl>${g}</Placemark>`);
  }
  const kml = `<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>${esc(baseName())}</name>
<Style id="np1"><LineStyle><color>ffa86224</color><width>2</width></LineStyle><PolyStyle><color>40a86224</color></PolyStyle></Style>
<Style id="np2"><LineStyle><color>ff9e4c7a</color><width>2</width></LineStyle><PolyStyle><color>409e4c7a</color></PolyStyle></Style>
<Style id="epi"><IconStyle><color>ff00d4ff</color></IconStyle></Style>${pm.join("\n")}</Document></kml>`;
  download(baseName() + ".kml", new Blob([kml], { type: "application/vnd.google-earth.kml+xml" }));
});
$("#dl-png").addEventListener("click", () => map.getCanvas().toBlob((b) => b && download(baseName() + "_map.png", b)));
$("#dl-txt").addEventListener("click", () => {
  const ev = S.caseObj.event, o = S.results.opts;
  const txt = [`GEOSLAM seismo-lineament analysis`, `Event: ${ev.title || ev.id} (${ev.id}), ${ev.time || ""}`,
    `Hypocentre: ${ev.lat}, ${ev.lon}, depth ${ev.depth_km} km (${o.depth_datum.replace("_", " ")})`,
    `Mechanism: ${S.caseObj.mechanisms[S.mech].label}`, `Method: ${o.method}, width multiplier ${o.multiplier}`,
    `Terrain: ${S.demInfo.source}, ${S.demInfo.cellsize} m cells, EPSG:${S.demInfo.epsg}`, "",
    ...S.results.planes.map((p) => p.summary + "\n"), `Share link: ${location.href}`, "",
    "Method: Cronin, V.S., Millard, M.A., Seidman, L.E., and Bayliss, B.G., 2008, Environmental & Engineering Geoscience 14(3):199-219.",
    "Software: GEOSLAM, https://geoslam.org"].join("\n");
  download(baseName() + "_summary.txt", new Blob([txt], { type: "text/plain" }));
});
$("#share").addEventListener("click", async () => {
  writeHash();
  try { await navigator.clipboard.writeText(location.href); setStatus("Share link copied. Anyone who opens it recomputes this map.", "ok"); }
  catch { setStatus("Copy the link from the address bar to share this map.", "ok"); }
});

/* ------------------------------------------------------------- share links */
function writeHash() {
  if (!S.caseObj) return;
  const p = new URLSearchParams();
  if (S.kind === "comcat") { p.set("e", S.caseObj.event.id); p.set("m", S.mech); }
  else {
    const ev = S.caseObj.event, np = S.caseObj.mechanisms[0].np1, u = S.caseObj.uncertainty;
    p.set("q", [ev.lat, ev.lon, ev.depth_km, np.strike, np.dip, np.rake, ev.mag ?? ""].join(","));
    if (u) p.set("u", [u.eh1, u.eh1_az, u.eh2, u.ez].join(","));
  }
  p.set("su", S.su); p.set("du", S.du); p.set("x", $("#mult").value);
  if (method() !== "envelope") p.set("meth", method());
  if (datum() !== "sea_level") p.set("datum", datum());
  if (num("#half-km")) p.set("half", num("#half-km"));
  if (num("#cell-m")) p.set("cell", num("#cell-m"));
  if (!S.caseObj.uncertainty) p.set("a", `${num("#a-eh") || 2000},${num("#a-ez") || 2000}`);
  history.replaceState(null, "", "#" + p.toString());
}

async function readHash() {
  if (location.hash.length < 3) return;
  const p = new URLSearchParams(location.hash.slice(1));
  if (p.get("meth")) $(`input[name="method"][value="${p.get("meth")}"]`).checked = true;
  if (p.get("datum")) $(`input[name="datum"][value="${p.get("datum")}"]`).checked = true;
  if (p.get("half")) $("#half-km").value = p.get("half");
  if (p.get("cell")) $("#cell-m").value = p.get("cell");
  if (p.get("x")) { $("#mult").value = p.get("x"); $("#mult").dispatchEvent(new Event("input")); }
  if (p.get("a")) { const [h, z] = p.get("a").split(","); $("#a-eh").value = h; $("#a-ez").value = z; }
  const keep = { mech: +(p.get("m") || 0), su: p.has("su") ? +p.get("su") : undefined, du: p.has("du") ? +p.get("du") : undefined };
  try {
    if (p.get("e")) {
      $("#event-id").value = p.get("e");
      setStatus(`Loading shared map for ${p.get("e")}…`);
      const id = p.get("e");
      const detail = await fetchJSON(`${COMCAT}?eventid=${encodeURIComponent(id)}&format=geojson`, `Event ${id}`);
      const urls = [];
      for (const t of ["focal-mechanism", "moment-tensor"]) for (const pr of detail.properties?.products?.[t] || []) { const u = pr.contents?.["quakeml.xml"]?.url; if (u) urls.push(u); }
      const qml = {};
      await Promise.all(urls.map(async (u) => { try { qml[u] = await (await fetch(u)).text(); } catch { /* skip */ } }));
      await pyReady;
      const c = JSON.parse(await pycall("case_from_comcat", JSON.stringify(detail), JSON.stringify(qml)));
      if (!c.mechanisms.length) throw new Error(`${id} has no focal mechanism.`);
      keep.mech = Math.min(keep.mech, c.mechanisms.length - 1);
      await setCase("comcat", c, keep);
    } else if (p.get("q")) {
      const [lat, lon, dep, s, d, r, mag] = p.get("q").split(",").map(parseFloat);
      Object.entries({ lat, lon, depth: dep, strike: s, dip: d, rake: r, mag }).forEach(([k, v]) => { if (Number.isFinite(v)) $("#m-" + k).value = v; });
      if (p.get("u")) { const [a, b, c2, z] = p.get("u").split(",").map(parseFloat); $("#m-eh1").value = a; $("#m-eh1az").value = b; $("#m-eh2").value = c2; $("#m-ez").value = z; }
      else { ["#m-eh1", "#m-eh2", "#m-ez"].forEach((s2) => { $(s2).value = ""; }); }
      $("#m-su").value = keep.su ?? 10; $("#m-du").value = keep.du ?? 10;
      $$(".tab").find((t) => t.dataset.tab === "manual").click();
      await useManual();
      if (keep.su != null) { S.su = keep.su; S.du = keep.du; $("#su").value = S.su; $("#du").value = S.du; $("#su").dispatchEvent(new Event("input")); $("#du").dispatchEvent(new Event("input")); }
    } else return;
    await runAll();
  } catch (e) { setStatus(e.message, "err"); }
}
