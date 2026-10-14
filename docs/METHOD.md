# The Seismo-Lineament Analysis Method (SLAM): mathematics, geophysics and implementation

*GEOSLAM 1.0 technical reference — compiled October 2026.*

SLAM was developed by **Vincent S. Cronin** (Baylor University) beginning in 2004. Its peer-reviewed description is Cronin, Millard, Seidman & Bayliss (2008), *Environmental & Engineering Geoscience* 14(3):199–219.

This document reconstructs the method from four kinds of source:

- **Cronin's current Mathematica notebook.** `SLAMcode20210103.nb` was begun on 5 Aug 2004, revised 3 Jan 2021, with minor changes on 26 Aug 2022. It was read directly from croninprojects.org and is the authoritative reference for every formula below.
- **Theses by his students.** Bayliss 2007, Worrell 2016, Rasaka 2016, Strasser 2017 and Ashburn 2015.
- **Cronin's abstracts and conference papers.** GSA 2007/2012/2014, INQUA PATA 2014, NAGT 2014.
- **The 2020 GEOSLAM port.** It is kept in `legacy/SLAM_2020.py`.

---

## 1. What SLAM does

**Assumptions.** A shallow (upper-crustal, roughly < 30 km) earthquake ruptures a fault. Assume the fault is approximately planar between the hypocentre and the ground surface. Then the fault's surface trace must lie where the fault plane intersects the topography.

**The two nodal planes.** A double-couple focal mechanism gives two orthogonal candidate planes. The radiation pattern alone cannot tell which one slipped, so SLAM evaluates both.

**Seismo-lineament.** Each candidate plane passes through the hypocentre. Projected to the ground surface (a DEM), it traces a line. Neither the hypocentre position nor the plane orientation is known exactly, so the line becomes a swath. Cronin calls this swath a **seismo-lineament**: the region of the ground surface that should contain the trace of the seismogenic fault, if that fault is emergent.

**Characteristic shapes.**

| Uncertainty included | Shape of the seismo-lineament |
|---|---|
| None | A single line |
| Location only | A straight band of constant width |
| Location and orientation | A "bow tie" with its narrowest point near the epicentre |

**Where SLAM fits.** SLAM is a reconnaissance tool that narrows the search area. It is followed by three further steps:

1. Geomorphic analysis: lineaments on low-sun hillshades, plus the Quaternary fault maps.
2. Comparison with aftershocks and other events: overlapping seismo-lineaments, and best-fit planes through the seismicity.
3. Fieldwork: fault-zone fabric, orientation and slip compared statistically with the nodal plane.

**Correlation test.** Cronin calls a mapped fault *correlated* with an earthquake only when all three of these hold (Strasser 2017; GSA 2007):

- it lies inside the seismo-lineament;
- its orientation is within the nodal-plane uncertainty;
- its slip direction is within the slip-vector (rake) uncertainty.

## 2. Coordinate frame

All computation happens in a right-handed metric frame:

- x = grid east, y = grid north, z = up, in metres.
- Positions are UTM coordinates on the zone of the study area.

**Cronin's frame.** He puts the origin at the DEM's lower-left corner at sea level:

```
focusMean = { utmE - xllcorner,  utmN - yllcorner,  focalDepth }     focalDepth = -1000 * depthKm
```

**Ground nodes.** For DEM node (i, j), with row 1 at the north edge:

```
x' = cellsize*(j-1) - focusMean_x
y' = cellsize*(nrows-i) - focusMean_y
z' = elevation[i,j] - focalDepth
```

Here `P = (x', y', z')` is the position of a ground node **relative to the hypocentre**. GEOSLAM uses cell centres, `x_min + (j+½)·cellsize`. Cronin's node registration differs by half a cell; `registration="cronin"` reproduces it exactly.

**Depth datum.**

- Cronin and his students treat catalog depth as depth below **sea level**, so `z_h = -depth`.
- Network catalogs are not consistent on this. Some depths are relative to the geoid. Others are relative to a station-elevation datum or the model surface.
- GEOSLAM defaults to Cronin's convention. It offers `depth_datum="surface"`, which measures depth below the ground at the epicentre.
- In high-relief areas the two choices can shift the trace by *Δz / tan(dip)*. For example, 1 km of relief on a 60° fault moves the trace about 580 m.

**Horizontal datum.**

- Cronin's `convertToUTM` uses the Clarke 1866 ellipsoid (NAD27).
- Modern hypocentres (ComCat, NCEDC, SCEDC, TexNet) and DEMs (3DEP, SRTM, Copernicus) are WGS84/NAD83, so GEOSLAM projects with WGS84.
- Mixing NAD27 epicentres with a NAD83 DEM introduces errors of 10–100 m in CONUS.

**Projection accuracy.** UTM is implemented with the third-order Krüger series (Karney 2011). It agrees with the numerically integrated meridian arc to a few micrometres.

## 3. True north versus grid north (meridian convergence)

**Why it matters.** Catalog azimuths (strike, dip direction, error-ellipse azimuth) are measured from geographic north, but the DEM grid is aligned with grid north.

**Cronin's correction:**

```
gridConvergence = ArcTan[ Tan[(midLong - zoneMeridian) Degree] * Sin[midLat Degree] ] * 180/Pi
npDipTr = npDipTrend - gridConvergence          (also eh1Az = eh1Azimuth - gridConvergence)
```

**Definition.** γ is the angle from true north to grid north, positive clockwise. It is positive east of the central meridian in the northern hemisphere, and then **grid azimuth = true azimuth − γ**.

**What GEOSLAM 1.0 computes.**

- It uses the exact ellipsoidal convergence at the epicentre, not the spherical formula at the map centre.
- The two differ by less than 0.001° at Truckee.
- γ is +1.805° at Truckee (zone 10, CM −123°) and +0.43° at South Napa.

**Note on the 2020 port.** It computed `(CM − lon)·sin(lat)` and *added* it. That gives the same result to first order.

## 4. A nodal plane as vectors

**Conventions.** Dip direction α (azimuth, clockwise from grid north) and dip δ (degrees below horizontal).

```
dip vector     d = makeVector[δ, α] = ( cosδ sinα,  cosδ cosα,  −sinδ )
strike (RHR)   s = α − 90°  →  ŝ = ( sin s, cos s, 0 )
upward normal  n = d × ŝ = ( sinδ sinα,  sinδ cosα,  cosδ )
```

**Code transcription.** Cronin builds `n` with `Cross[npDipVect, npStrikeVect]`. GEOSLAM uses the closed form. A test checks that the two agree to 1e-12 over random orientations.

**Signed distance.** For a ground node with relative position P, the signed distance from the candidate fault plane is

    f = n · P.

The plane's ground trace is the contour `f = 0`.

**Fold-over for steep and shallow planes.**

- Cronin handles dips outside 0–90° explicitly, because δ ± σ can leave that range:
  - if δ − σ < 0, use |δ − σ| and add 180° to α;
  - if δ + σ > 90, use 180 − (δ + σ) and add 180° to α.
- The closed form of `n` is continuous for any real δ, and `n(δ, α) = −n(−δ, α+180°)`. So GEOSLAM can use δ ± σ directly. It describes the same plane, and distances are taken in absolute value.

**Auxiliary plane.** The second nodal plane follows by swapping the fault normal and the slip vector. The slip vector uses Aki & Richards' (strike, dip, rake) convention. It is checked against ComCat, for example 345/85/−170 ↔ 254/80/−5 for South Napa, and against TexNet.

## 5. Hypocentre uncertainty and swath half-width

**Uncertainty model.** The location uncertainty is an ellipsoid. Its semi-axes are aᵢ, along unit principal axes eᵢ.

**Tangent planes.** For a candidate fault plane with unit normal n, Cronin bounds the plane's position by the two planes parallel to it that are **tangent to the ellipsoid** (`findHalfWidth`). His module does three things:

1. It rotates n into the ellipsoid frame by `eh1Az − 90°`.
2. It finds the tangent point `t_i = λ n_i a_i²`, with λ chosen so that t lies on the ellipsoid.
3. It returns `n · t`.

**Closed form.** Algebraically, the result is the ellipsoid's *support function*:

    h(n) = sqrt( Σ (n·eᵢ)² aᵢ² ) = sqrt( nᵀ C n ),   C = R diag(a²) Rᵀ

**Slab test.** A ground node lies on the uncertainty slab of that plane when |n·P| ≤ h(n).

**Limits.** For horizontal error eh and vertical error ez (isotropic in the horizontal):

    h = sqrt( eh² sin²δ + ez² cos²δ )

- A vertical fault takes the horizontal error.
- A horizontal plane takes the vertical error.

**Bayliss's earlier box model.** Bayliss (2007, Eq. 7) used a box rather than an ellipsoid: `w = ez cosδ + eh sinδ`. This is always at least as large as the ellipsoid width. It is kept as `box_half_width_bayliss` for comparison.

**What GEOSLAM adds.** Cronin's sheet assumes a *horizontal* ellipse (eh1 along eh1Az, eh2 perpendicular) plus a vertical ez. GEOSLAM accepts the full 3-D QuakeML confidence ellipsoid (azimuth, plunge, rotation) as well. This matters in practice:

- The 2014 South Napa NCSN origin has semi-axes of 384/264/168 m.
- Its *major axis plunges 86°*, so it is nearly vertical.
- Treating 384 m as a horizontal semi-axis, as the 2020 port did, roughly doubles the horizontal half-width of a vertical fault.

**When nothing is known.**

- Cronin's thin trace uses `widthFactor · cellsize · sinδ`, with widthFactor = 1.5. That is a band about 1.5 DEM cells wide on the ground.
- Theses have used informal assumed values: ±2 km on every axis for the 1966 Truckee event (Strasser), and 10° on every angle for Oklahoma (Rasaka).
- GEOSLAM exposes `unknown_uncertainty` and `default_plane_uncertainty` for this purpose. It reports what was assumed.

**Confidence level.** These widths are as large as the published ellipse. Catalogs differ in what that means:

| Source | Confidence level |
|---|---|
| ComCat / NCSN confidence ellipsoids | 95 % |
| OGS errors | 90 % |
| Hypoinverse ERH/ERZ | Closer to 68 % |

`multiplier` rescales all widths. The theses used ×1 to ×1.5.

## 6. Orientation uncertainty: from planes to a swath

**Cronin's discrete planes.** The notebook evaluates seven slabs and overlays their ground intersections. Here σδ is the dip uncertainty and σα the dip-direction uncertainty, which equals the strike uncertainty.

| # | Name | Dip | Dip direction | Half-width |
|---|---|---|---|---|
| 1 | MThin | δ | α | widthFactor·cellsize·sinδ (the mean trace) |
| 2 | GS — gentle, smaller azimuth | δ−σδ | α−σα | h(n) |
| 3 | GL — gentle, larger azimuth | δ−σδ | α+σα | h(n) |
| 4 | SS — steep, smaller azimuth | δ+σδ | α−σα | h(n) |
| 5 | SL — steep, larger azimuth | δ+σδ | α+σα | h(n) |
| 6 | SM — steep, mean azimuth | δ+σδ | α | h(n) |
| 7 | Mean | δ | α | h(n) |

**How Cronin combines them.** He sums the binary grids and contours the union (`summary`). The map shows each swath's boundary lines over a hillshade. The interpreter reads the seismo-lineament as the region **bounded by the outermost lines**.

**Two practical consequences.**

- **Gaps between lines.** Location errors are now often only tens of metres (double-difference relocations). Then the slabs are thin, and the region *between* the gentle and steep planes is not covered by any slab. That region still holds plausible fault traces. The 2020 port tried to fill it with row-by-row "fill" heuristics (`fill_swaths`, `fill_corners`, `fill_middle`, plus a manual `corners_fill` argument). Those heuristics are not part of Cronin's method and fail for many orientations.
- **No "gentle, mean azimuth" plane.** It is absent from Cronin's set.

**GEOSLAM's continuous envelope** (`method="envelope"`, the default) computes the region Cronin draws by hand.

*Definition.* The seismo-lineament is the set of ground points P such that at least one plane with δ' ∈ [δ−σδ, δ+σδ] and α' ∈ [α−σα, α+σα] satisfies |n(δ',α')·P| ≤ h(n(δ',α')).

*Algorithm.* For fixed P, consider f(δ,α) = n·P = sinδ (Pₓ sinα + P_y cosα) + cosδ P_z. Its only interior critical points are two: where n ∥ ±P (f = ±|P|, the extrema), and saddles on horizontal planes. So |f| has **no interior minimum except its zeros**. A point is therefore inside when either condition holds:

1. f changes sign over the rectangle, meaning a plausible plane passes *through* it. The minimum and maximum of f are taken over a dense sampling of the rectangle's boundary plus the analytic critical points.
2. |f| ≤ h on some boundary sample.

The only approximation is the slow variation of h inside the rectangle.

*Verification.* A unit test compares the result with brute-force sampling of 97 × 121 interior planes. They disagree on fewer than 0.2 % of swath cells. The envelope always contains Cronin's 7-slab union.

`method="cronin"` remains available to reproduce published maps exactly.

## 7. Map size, DEM and hillshade

**Map size.** Cronin's minimum radius is d = (mean elevation + |depth|) / tan(gentlest plausible dip). `suggest_half_extent` adds the location error and a margin of 1.5.

**DEM resolution.** Theses used several resolutions:

- 10 m NED / 3DEP for SLAM;
- 1–2 m lidar for the geomorphology;
- thinning ("SLAM data diet") to keep Mathematica responsive.

Swath geometry is insensitive to the DEM cell size once cells are smaller than h. A 30 m cell is plenty for location errors of 100 m or more. Use 10 m or lidar for interpretation.

**Hillshade.** Cronin lights the surface from the direction **opposite the nodal plane's dip** at 15° sun elevation (`LightingAngle -> {lightDirection, Pi/12}`). That emphasises scarps that face the fault. The 2020 port passed a math-convention angle to matplotlib's compass-convention `azdeg`. GEOSLAM 1.0 uses the compass azimuth.

## 8. Data sources (verified October 2026)

| Need | Source | Notes |
|---|---|---|
| Hypocentres and mechanisms | USGS ComCat FDSN | GeoJSON `horizontal/vertical-error` are in **km**. `error-ellipse-*` axes are in **metres** and form a 3-D ellipsoid (95 %). Nodal-plane **uncertainties are only in each product's own `quakeml.xml`**: not in the GeoJSON, and not in the aggregated QuakeML. |
| N. California | NCEDC FDSN `includemechanisms=true` | Hypoinverse ellipse (min/max/azimuth) plus depth error, i.e. Cronin's eh1/eh2/eh1Az/ez. Relocations: Waldhauser NCAeqDD (LDEO). |
| S. California | SCEDC FDSN; Hauksson–Yang–Shearer relocations; YHS HASH mechanisms | HASH files carry fault- and auxiliary-plane uncertainty plus quality A–D. |
| Texas | TexNet ArcGIS MapServer layer 0 | Location errors and mechanisms; **no plane uncertainty** (assume one, e.g. 10°). |
| Oklahoma | OGS catalog (ComCat `catalog=ok`) | 90 % errors in km; no mechanisms. |
| DEM, US | USGS 3DEP ImageServer, exported directly in UTM | No key. Up to 8000 × 8000 px. |
| DEM, global | OpenTopography (`COP30`, `SRTMGL1`, `NASADEM`, …) | **Free API key required since 2022.** |
| DEM, global, no key | Copernicus GLO-30 COGs on AWS | A surface model (includes canopy and buildings). |
| Faults, US | USGS Qfaults ArcGIS REST, layer 21 | Bounding-box query. |
| Faults, global | GEM Global Active Faults | CC-BY-SA 4.0. |
| Basemap tiles | Stamen tiles are gone (Oct 2023); Stadia needs a key | GEOSLAM uses the DEM hillshade instead. |

## 9. Validation: 2014 M6.0 South Napa

**Inputs.**

- DEM: USGS 3DEP, 150 m, UTM 10N.
- Faults: USGS Quaternary faults.
- Aftershocks: 30 days from ComCat.
- Hypocentre and mechanism (Worrell 2016, Tables 1–2): double-difference hypocentre, eh = 28 m, ez = 29 m. NCSN FPFIT mechanism: dip direction 075, dip 85, rake −170, with σ of 5° on strike and 10° on dip.

**Result:**

| Case | Method | Swath area km² | Fraction of mapped West Napa fault inside |
|---|---|---|---|
| NP1 (strike 345), Worrell inputs | Cronin 7 slabs | 22.9 | 0.20 |
| NP1 (strike 345), Worrell inputs | envelope | 143.8 | **0.79** |
| NP2 (strike 254, auxiliary) | envelope | 111.0 | 0.10 |
| NP1, all inputs from ComCat | envelope | 178.8 | 0.69 |

**Interpretation.**

- The NP1 seismo-lineament contains the West Napa fault system and the 2014 ruptures, as Worrell found.
- The auxiliary plane does not.
- The 7-slab union, taken literally as a raster, covers only a fifth of the fault, because the 28 m slabs are thin lines. That is why the envelope is the default.

`tests/test_validation_napa.py` locks these results in.

## 10. Known limitations

- **Planarity and emergence.** SLAM assumes a planar fault that reaches the surface. Listric, segmented or blind faults violate this.
- **Gentle dips.** Results degrade quickly for gentle dips, because swath width grows as 1/tanδ.
- **Depth datum.** The datum ambiguity (§2) is often the largest unmodelled error in high-relief terrain.
- **Mechanism uncertainty.** Catalog uncertainties are often formal or missing, and FPFIT's ± values are not confidence intervals.
- **Interpretation.** A swath containing a fault is necessary but not sufficient. Orientation, slip and geomorphology must also agree.

## References

- Bayliss, B.G. (2007). MS thesis, Baylor University.
- Cronin, V.S. (2004–2022). SLAMcode20210103.nb. https://croninprojects.org/Vince/SLAM/
- Cronin, V.S., Millard, M.A., Seidman, L.E., Bayliss, B.G. (2008). Environmental & Engineering Geoscience 14(3):199–219.
- Cronin, V.S. (2014). INQUA PATA Days paper.
- Karney, C.F.F. (2011). J. Geodesy 85:475–485.
- Rasaka, B.M. (2016). MS thesis, Baylor University.
- Snyder, J.P. (1982, 1987). USGS Bulletin 1532; Professional Paper 1395.
- Strasser, M. (2017). MS thesis, Baylor University.
- Worrell, V.E. (2016). MS thesis, Baylor University.
- Yang, Hauksson & Shearer (2012). BSSA 102:1179.
- Hauksson, Yang & Shearer (2012). BSSA 102:2239.
- Waldhauser (2009). BSSA.
- Styron & Pagani (2020). Earthquake Spectra 36(1_suppl).
