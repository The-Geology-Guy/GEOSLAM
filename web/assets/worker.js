/* Runs the GEOSLAM Python package in a Web Worker via Pyodide. */
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/";
importScripts(PYODIDE + "pyodide.js");

let py = null;

const GLUE = `
import json
import numpy as np
from geoslam import web

def _f32(buf):
    return np.asarray(buf.to_py(), dtype=np.float32)

def _set_dem_utm(buf, nrows, ncols, x_min, y_max, cellsize, epsg, source):
    return web.set_dem_utm(_f32(buf), int(nrows), int(ncols), x_min, y_max, cellsize, epsg, source)

def _set_dem_geo(buf, nrows, ncols, lon_min, lat_max, dlon, dlat, clon, clat, cellsize, box_json, source):
    box = json.loads(box_json) if box_json else None
    return web.set_dem_geographic(_f32(buf), int(nrows), int(ncols), lon_min, lat_max, dlon, dlat,
                                  clon, clat, cellsize, box, source)

def _hillshade(az, alt):
    b, w, h = web.hillshade(az, alt)
    return np.frombuffer(b, np.uint8)

def _beachball(s, d, r, n):
    b, n = web.beachball(s, d, r, int(n))
    return np.frombuffer(b, np.uint8)
`;

async function init(wheelUrl) {
  py = await loadPyodide({ indexURL: PYODIDE });
  post("progress", "Loading numerical libraries");
  await py.loadPackage(["numpy", "scipy", "contourpy", "micropip"]);
  post("progress", "Installing GEOSLAM");
  const micropip = py.pyimport("micropip");
  await micropip.install.callKwargs(new URL(wheelUrl, self.location.href).href, { deps: false });
  py.runPython(GLUE);
  return py.runPython("import geoslam; geoslam.__version__");
}

function post(type, data) { self.postMessage({ type, data }); }

function bytesOf(proxy) {
  const b = proxy.getBuffer("u8");
  const out = new Uint8ClampedArray(b.data);
  b.release();
  proxy.destroy();
  return out;
}

const handlers = {
  init: ([wheelUrl]) => init(wheelUrl),
  // generic: geoslam.web.<fn>(*args) with JSON-friendly args; returns str/number
  call: ([fn, ...args]) => {
    const f = py.pyimport("geoslam.web")[fn];
    const r = f(...args);
    if (r && typeof r === "object" && r.toJs) {
      const js = r.toJs({ dict_converter: Object.fromEntries, create_pyproxies: false });
      r.destroy();
      return js;
    }
    return r;
  },
  setDemUtm: ([buf, ...rest]) => py.globals.get("_set_dem_utm")(buf, ...rest),
  setDemGeo: ([buf, ...rest]) => py.globals.get("_set_dem_geo")(buf, ...rest),
  setDemAscii: ([text, zone, north]) => py.pyimport("geoslam.web").set_dem_ascii(text, zone, north),
  hillshade: ([az, alt]) => bytesOf(py.globals.get("_hillshade")(az, alt)),
  beachball: ([s, d, r, n]) => bytesOf(py.globals.get("_beachball")(s, d, r, n)),
};

self.onmessage = async (ev) => {
  const { id, op, args } = ev.data;
  try {
    const result = await handlers[op](args || []);
    const transfer = result instanceof Uint8ClampedArray ? [result.buffer] : [];
    self.postMessage({ id, ok: true, result }, transfer);
  } catch (err) {
    let msg = String(err && err.message ? err.message : err);
    const lines = msg.trim().split("\n");
    self.postMessage({ id, ok: false, error: lines[lines.length - 1], detail: msg });
  }
};
