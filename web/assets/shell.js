/* GEOSLAM "98" shell: windows, menus, taskbar, Start menu, clock and map status bar.
   It only drives the interface; all analysis stays in app.js. */

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const desktop = $("#desktop");
const mobile = () => matchMedia("(max-width: 760px)").matches;
const resizeMap = () => setTimeout(() => window.geoslamMap?.resize(), 0);

/* ----------------------------------------------------------------- windows */
const WIN = { "win-app": "min-app", "win-map": "min-map" };

function setActive(win) {
  $$(".desktop > .window").forEach((w) => w.classList.toggle("is-active", w === win));
}
function syncTasks() {
  for (const [id, cls] of Object.entries(WIN)) {
    const minimized = desktop.classList.contains(cls);
    const task = $(`.task[data-task="${id}"]`);
    task.classList.toggle("is-on", !minimized && $("#" + id).classList.contains("is-active"));
    task.setAttribute("aria-pressed", String(!minimized));
  }
}
function minimize(id) {
  const win = $("#" + id);
  win.classList.remove("is-max");
  desktop.classList.add(WIN[id]);
  const other = Object.keys(WIN).find((k) => k !== id);
  if (!desktop.classList.contains(WIN[other])) setActive($("#" + other));
  syncTasks(); resizeMap();
}
function restore(id) {
  desktop.classList.remove(WIN[id]);
  setActive($("#" + id));
  syncTasks(); resizeMap();
}
function toggleMax(id) {
  const win = $("#" + id);
  const on = !win.classList.contains("is-max");
  $$(".desktop > .window").forEach((w) => w.classList.remove("is-max"));
  win.classList.toggle("is-max", on);
  win.querySelector(".tb-max")?.setAttribute("aria-label", on ? "Restore" : "Maximize");
  setActive(win); syncTasks(); resizeMap();
}

for (const id of Object.keys(WIN)) {
  const win = $("#" + id);
  win.addEventListener("pointerdown", () => { setActive(win); syncTasks(); });
  win.addEventListener("focusin", () => { setActive(win); syncTasks(); });
  win.querySelectorAll(".titlebar [data-win]").forEach((b) => b.addEventListener("click", (e) => {
    e.stopPropagation();
    if (b.dataset.win === "max") toggleMax(id); else minimize(id);
  }));
  win.querySelector(".titlebar").addEventListener("dblclick", (e) => { if (!e.target.closest("button") && !mobile()) toggleMax(id); });
}
$$(".task").forEach((t) => t.addEventListener("click", () => {
  const id = t.dataset.task;
  if (desktop.classList.contains(WIN[id])) restore(id);
  else if ($("#" + id).classList.contains("is-active")) minimize(id);
  else { setActive($("#" + id)); syncTasks(); }
}));
syncTasks();

/* ------------------------------------------------------------------- menus */
const hasResults = () => !$("#results").hidden;
let openMenu = null;

function closeMenus() {
  $$(".menu-pop").forEach((p) => { p.hidden = true; });
  $$(".menu-btn").forEach((b) => b.setAttribute("aria-expanded", "false"));
  $("#start-menu").hidden = true;
  $("#start").setAttribute("aria-expanded", "false");
  openMenu = null;
}
function refreshItems(root) {
  root.querySelectorAll("[data-needs='results']").forEach((b) => { b.disabled = !hasResults(); });
  root.querySelectorAll("[data-basemap]").forEach((b) => b.setAttribute("aria-checked", String($("#basemap").value === b.dataset.basemap)));
  root.querySelectorAll("[data-toggle]").forEach((b) => b.setAttribute("aria-checked", String($(b.dataset.toggle).checked)));
}
function openPop(btn, pop) {
  closeMenus();
  refreshItems(pop);
  pop.hidden = false;
  btn.setAttribute("aria-expanded", "true");
  openMenu = pop;
}

$$(".menu-btn").forEach((btn) => {
  const pop = btn.nextElementSibling;
  btn.addEventListener("click", (e) => {
    e.stopPropagation();
    if (!pop.hidden) { closeMenus(); return; }
    openPop(btn, pop);
    pop.querySelector("button:not(:disabled), a")?.focus();
  });
  // Windows behaviour: once one menu is open, hovering another opens it
  btn.addEventListener("pointerenter", () => { if (openMenu && openMenu !== pop && openMenu.classList.contains("menu-pop")) openPop(btn, pop); });
});

const startBtn = $("#start"), startMenu = $("#start-menu");
startBtn.addEventListener("click", (e) => {
  e.stopPropagation();
  if (!startMenu.hidden) { closeMenus(); return; }
  closeMenus();
  refreshItems(startMenu);
  startMenu.hidden = false;
  startBtn.setAttribute("aria-expanded", "true");
  openMenu = startMenu;
  startMenu.querySelector("button:not(:disabled), a")?.focus();
});

// menu commands; [data-open] items are wired by app.js and only need the menu closed
document.addEventListener("click", (e) => {
  const item = e.target.closest(".menu-pop [role^='menuitem'], .start-items [role='menuitem']");
  if (item) {
    if (item.dataset.click) $(item.dataset.click)?.click();
    if (item.dataset.basemap) { const s = $("#basemap"); s.value = item.dataset.basemap; s.dispatchEvent(new Event("change")); }
    if (item.dataset.toggle) { const c = $(item.dataset.toggle); c.checked = !c.checked; c.dispatchEvent(new Event("change")); }
    if (item.dataset.action === "restart") { location.href = location.pathname; }
    closeMenus();
    return;
  }
  if (openMenu && !e.target.closest(".menu-pop, .start-menu")) closeMenus();
});
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && openMenu) {
    const owner = openMenu === startMenu ? startBtn : openMenu.previousElementSibling;
    closeMenus(); owner?.focus();
  }
  if (!openMenu || !["ArrowDown", "ArrowUp"].includes(e.key)) return;
  const items = [...openMenu.querySelectorAll("button:not(:disabled), a")];
  const i = items.indexOf(document.activeElement);
  const next = items[(i + (e.key === "ArrowDown" ? 1 : -1) + items.length) % items.length];
  next?.focus(); e.preventDefault();
});

/* ------------------------------------------------------------ map window bits */
const BASE_NAMES = { topo: "USGS topographic", relief: "USGS shaded relief", imagery: "USGS imagery", osm: "OpenStreetMap" };
$("#basemap").addEventListener("change", (e) => { $("#t-map").textContent = "Map - " + BASE_NAMES[e.target.value]; });

function fmtLonLat(ll) {
  const la = `${Math.abs(ll.lat).toFixed(4)}°${ll.lat >= 0 ? "N" : "S"}`;
  const lo = `${Math.abs(ll.lng).toFixed(4)}°${ll.lng >= 0 ? "E" : "W"}`;
  return `${la}, ${lo}`;
}
(function hookMap(tries = 0) {
  const map = window.geoslamMap;
  if (!map) { if (tries < 100) setTimeout(() => hookMap(tries + 1), 100); return; }
  const zoom = () => { $("#sb-zoom").textContent = "Zoom " + map.getZoom().toFixed(1); };
  let hovering = false;
  const centre = () => { if (!hovering) $("#sb-coords").textContent = fmtLonLat(map.getCenter()) + " (centre)"; };
  map.on("mousemove", (e) => { hovering = e.originalEvent?.pointerType !== "touch"; if (hovering) $("#sb-coords").textContent = fmtLonLat(e.lngLat); });
  map.on("mouseout", () => { hovering = false; centre(); });
  map.on("move", centre);
  map.on("zoom", zoom);
  zoom();
  $("#sb-coords").textContent = fmtLonLat(map.getCenter()) + " (centre)";
  new ResizeObserver(() => map.resize()).observe($(".map-wrap"));
})();

/* ------------------------------------------------------- collapsible legend */
const legend = $("#legend"), legendBtn = legend.querySelector(".legend-title");
function setLegend(open) {
  legend.classList.toggle("is-collapsed", !open);
  legendBtn.setAttribute("aria-expanded", String(open));
  legendBtn.title = open ? "Hide the legend" : "Show the legend";
}
legendBtn.addEventListener("click", () => setLegend(legend.classList.contains("is-collapsed")));
setLegend(!mobile());   // phones start with it folded so the map stays visible

/* --------------------------------------------- status bar mirrors the busy box */
const sbApp = $("#sb-app");
new MutationObserver(() => {
  const busy = !$("#busy").hidden;
  sbApp.textContent = busy ? ($("#busy-text").textContent || "Working…") : (hasResults() ? "Done" : "Ready");
}).observe($("#busy"), { attributes: true, subtree: true, childList: true, characterData: true });

/* ------------------------------------------------------------------- clock */
const clock = $("#clock");
function tick() {
  const now = new Date();
  clock.textContent = now.toLocaleTimeString([], { hour: "numeric", minute: "2-digit" });
  clock.dateTime = now.toISOString();
  clock.title = now.toLocaleDateString([], { weekday: "long", year: "numeric", month: "long", day: "numeric" });
}
tick();
setInterval(tick, 15000);
