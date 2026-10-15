# GEOSLAM — Seismo-Lineament Analysis Method in Python

GEOSLAM is a Python implementation of **Vince Cronin's Seismo-Lineament Analysis Method (SLAM)**. SLAM projects the nodal planes of a shallow earthquake to the ground surface, along with the uncertainty in the hypocentre location and in the plane orientation. The result outlines the swath (the *seismo-lineament*) that should contain the surface trace of the fault that generated the earthquake.

> SLAM is the intellectual property and work product of Vince Cronin (Baylor University), developed with many students and collaborators — see https://croninprojects.org/Vince/SLAM/. GEOSLAM is a Python translation by Luke Pajer.

**Version 1.0 (2026) is a rebuild.** The 2020 release (0.0.x) is preserved in `legacy/`. Version 1.0 differs from it in four ways:

- **The math follows Cronin's current notebook line by line.** The source is `SLAMcode20210103.nb`, revised in 2022. The details are in [`docs/METHOD.md`](docs/METHOD.md), and tests check the code against his formulas.
- **Every DEM is reprojected to a true UTM grid.** The 2020 code approximated cell sizes in lat/lon.
- **The seismo-lineament is a proper envelope.** It covers every plane within the dip and strike uncertainty. This replaces the manual `fill_*`/`corners_fill` heuristics, and Cronin's exact 7-slab method is still available.
- **The data connectors work with today's services.** They handle the 3-D ComCat error ellipsoid and its units, and they read nodal-plane uncertainties from per-product QuakeML. Supported sources are NCEDC, SCEDC, TexNet, USGS 3DEP, OpenTopography (API key), Copernicus DEM, and the USGS Qfaults REST service.

## Install

```bash
pip install "git+https://github.com/The-Geology-Guy/GEOSLAM@rebuild-2026"
# optional extras: [gis] (rasterio, pyproj), [tables] (pandas, openpyxl), [test]
```

The core needs numpy, scipy, matplotlib, contourpy, requests and tifffile. No GDAL is required.

## Quick start

```bash
geoslam run nc72282711 --out napa/          # 2014 South Napa M6.0, both nodal planes
```

```python
from geoslam.workflow import analyze_event
res = analyze_event("nc72282711", cellsize=30)   # ComCat id -> DEM -> both planes -> faults
fig = res.plot(); res.save("napa/")
for r in res.results:
    print(r.summary())
```

### Bring your own data (Cronin's workflow)

```python
from geoslam import DEMGrid, UTMZone, Hypocenter, HypocenterUncertainty, NodalPlane, seismo_lineament
from geoslam.plotting import plot_slam

dem = DEMGrid.from_arc_ascii("DogValley.asc", UTMZone(10))           # Esri ASCII grid, UTM metres
hyp = Hypocenter(-120.2063, 39.4122, 8.727,
                 HypocenterUncertainty.from_cronin(eh1=300, eh1_azimuth=45, eh2=300, ez=800))
plane = NodalPlane.from_dip_azimuth(300, 80, rake=20, dip_azimuth_uncertainty=8, dip_uncertainty=13)
r = seismo_lineament(dem, hyp, plane)            # method="cronin" reproduces the 7-slab union
plot_slam(r)
```

Cronin's input spreadsheet layout can be read with `geoslam.io.read_cronin_table`.

### Report figures

`geoslam.report.render` draws publication-style maps in the manner of Cronin's SLAM
figures (terrain filling the frame, latitude/longitude on all edges, dashed lineament
edges, black mean trace, labelled beachballs A1, A2, B1…) as PNG, PDF or SVG:

```python
from geoslam.report import render, ReportItem, ReportOptions
png = render(dem, [ReportItem("A1", "A", r)], ReportOptions(fmt="pdf", legend=True), faults=faults)
```

The web app at https://geoslam.org exposes the same figures under **Report view & download**,
along with the terrain image (PNG + world file) and the DEM as an Esri ASCII grid.

## Data sources and API keys

| Data | Default source | Key |
|---|---|---|
| Events and mechanisms | USGS ComCat (also NCEDC, SCEDC, TexNet, ISC, GEOFON via `geoslam.sources`) | none |
| DEM (US) | USGS 3DEP ImageServer | none |
| DEM (global) | OpenTopography (`COP30`, `SRTMGL1`, …) | `OPENTOPOGRAPHY_API_KEY` |
| DEM (global, no key) | Copernicus GLO-30 on AWS | none |
| Faults | USGS Quaternary Fault & Fold DB (REST); GEM Global Active Faults | none |

Please cite the data providers. Citations are listed in `docs/METHOD.md` §8.

## Tests

```bash
pytest                                   # offline: math, parsers, South Napa regression
GEOSLAM_NETWORK_TESTS=1 pytest -m network   # live ComCat / NCEDC / 3DEP / Qfaults
```

## References and credits

- Cronin, V.S., Millard, M.A., Seidman, L.E., Bayliss, B.G., 2008, The Seismo-Lineament Analysis Method (SLAM): a reconnaissance tool to help find seismogenic faults: *Environmental & Engineering Geoscience* 14(3):199–219.
- Theses used for validation:
  - Worrell (2016, South Napa)
  - Ashburn (2015) and Strasser (2017) (Dog Valley fault, Truckee)
  - Rasaka (2016, central Oklahoma)
  - Bayliss (2007, Santa Monica Mountains)

Contacts:
- GEOSLAM: Luke Pajer — ORCID 0000-0002-5218-7650
- SLAM: Vince Cronin — ORCID 0000-0002-3069-6470

MIT License. Contributions welcome under the [code of conduct](CODE_OF_CONDUCT.md).
