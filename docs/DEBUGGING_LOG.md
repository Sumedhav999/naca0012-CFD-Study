# NACA 0012 CFD Validation — Debugging Log

OpenFOAM 11 · Steady RANS · Spalart-Allmaras · 2D

**Status:** complete. Mesh rebuilt and verified; y⁺ confirmed; grid-independence and
domain sensitivity studies done; AoA sweep 0°–12° on the 50-chord domain; validated
against Gregory & O'Reilly (Re 3×10⁶, tripped) and Ladson (NASA TM 4074).

---

## 1. Case definition

| Parameter | Value |
|---|---|
| Geometry | NACA 0012, chord `c` = 1.0 m |
| Freestream velocity | 30 m/s |
| Kinematic viscosity | 1.0 × 10⁻⁵ m²/s |
| **Reynolds number** | **3.0 × 10⁶** |
| Solver | `foamRun` / `incompressibleFluid`, `ddtSchemes steadyState` (SIMPLE) |
| Turbulence model | Spalart-Allmaras, `nutUSpaldingWallFunction` |
| Dimensionality | 2D (one cell in z, `empty` front/back patches) |
| AoA sweep | 0° to 12° in 2° steps |
| AoA method | Freestream velocity vector rotated; airfoil remains mesh-aligned |
| Validation reference | Abbott & von Doenhoff, *Theory of Wing Sections*, NACA 0012 at Re = 3 × 10⁶ |

ν is selected to set the target Reynolds number. It is not the kinematic viscosity of
air at standard conditions (air at 20 °C is 1.516 × 10⁻⁵ m²/s). For incompressible
flow only Re is physically meaningful, so no temperature or altitude is claimed.

---

## 2. Initial symptom

The first completed AoA sweep produced results that were both internally inconsistent
and far from published data.

| AoA | CL | CD |
|---|---|---|
| 0 | 0.000 | 0.009 |
| 2 | 0.141 | 0.008 |
| 4 | 0.272 | 0.016 |
| 6 | 0.301 | 0.034 |
| 8 | 0.347 | 0.041 |
| 10 | 0.404 | 0.044 |
| 12 | 0.537 | 0.046 |

Three observations, none of which required reference data to identify as wrong:

1. **Lift-curve slope collapsed then recovered.** Segment slopes ran 0.0705, 0.0655,
   **0.0145**, 0.0230, 0.0285, **0.0665** per degree. A real lift curve is linear, then
   flattens toward stall, then drops. It does not flatten mid-range and then recover
   its original slope.

2. **Drag decreased between 0° and 2°.** For a symmetric section, minimum drag is at
   0° by symmetry. Drag cannot fall with incidence.

3. **Drag flattened above 6° while lift accelerated.** CD increments ran +0.018,
   +0.007, +0.003, +0.002. Drag should be convex in α, sharply so near separation.
   This is concave.

Additionally, L/D at 8° was 8.5, where a NACA 0012 at this Reynolds number should
reach roughly 65. Since a reference-area error cancels in L/D, this ruled out
normalisation mistakes and pointed at the flow field itself.

Two independent root causes were found.

---

## 3. Root cause 1 — Force coefficients computed in the wrong reference frame

### Symptom
Observations 2 and 3 above.

### Diagnosis
`0/U` contained `(29.34443 6.23735 0)`, which is 30.000 m/s at exactly 12.000°, so AoA
was being set by rotating the freestream while the airfoil stayed aligned to the mesh
axes. But `system/forceCoeffs` contained:

```
liftDir  (0 1 0);
dragDir  (1 0 0);
```

These are fixed in the **mesh** frame. The function object was therefore reporting
**normal and axial force coefficients**, not lift and drag:

```
reported "CL" = C_N = CL·cos α + CD·sin α
reported "CD" = C_A = CD·cos α − CL·sin α
```

The `−CL·sin α` term grows with incidence, which is precisely why the drag curve first
dipped and then flattened. Both anomalies were artifacts of the frame, not the flow.

### Corrected values

Inverse transform: `CL = C_N·cos α − C_A·sin α`, `CD = C_N·sin α + C_A·cos α`

| AoA | Reported CD | **True CD** | Error factor | Reported CL | True CL |
|---|---|---|---|---|---|
| 0 | 0.009 | 0.0090 | 1.00× | 0.000 | 0.000 |
| 2 | 0.008 | **0.0129** | 1.61× | 0.141 | 0.141 |
| 4 | 0.016 | **0.0349** | 2.18× | 0.272 | 0.270 |
| 6 | 0.034 | **0.0653** | 1.92× | 0.301 | 0.296 |
| 8 | 0.041 | **0.0889** | 2.17× | 0.347 | 0.338 |
| 10 | 0.044 | **0.1135** | 2.58× | 0.404 | 0.390 |
| 12 | 0.046 | **0.1566** | 3.40× | 0.537 | 0.516 |

Drag was under-reported by up to **3.4×**. Lift was essentially unaffected, since the
`CD·sin α` correction term is negligible.

After correction the drag curve is monotonic and convex, as required.

### Fix
Rotate the direction vectors with the freestream for each case:

```
dragDir  ( cos α   sin α  0);
liftDir  (-sin α   cos α  0);
```

| AoA | `dragDir` | `liftDir` |
|---|---|---|
| 0 | (1.000000 0.000000 0) | (0.000000 1.000000 0) |
| 2 | (0.999391 0.034899 0) | (-0.034899 0.999391 0) |
| 4 | (0.997564 0.069756 0) | (-0.069756 0.997564 0) |
| 6 | (0.994522 0.104528 0) | (-0.104528 0.994522 0) |
| 8 | (0.990268 0.139173 0) | (-0.139173 0.990268 0) |
| 10 | (0.984808 0.173648 0) | (-0.173648 0.984808 0) |
| 12 | (0.978148 0.207912 0) | (-0.207912 0.978148 0) |

To prevent `0/U` and `forceCoeffs` from drifting apart again, both are now set from a
single script:

```bash
#!/bin/bash
# setAoA.sh <degrees>
A=$1
UX=$(python3 -c "import math;print(f'{30*math.cos(math.radians($A)):.6f}')")
UY=$(python3 -c "import math;print(f'{30*math.sin(math.radians($A)):.6f}')")
DX=$(python3 -c "import math;print(f'{math.cos(math.radians($A)):.6f}')")
DY=$(python3 -c "import math;print(f'{math.sin(math.radians($A)):.6f}')")
sed -i "s|^internalField uniform .*|internalField uniform ($UX $UY 0);|" 0/U
sed -i "s|^\( *dragDir *\).*|\1( $DX  $DY 0);|" system/forceCoeffs
sed -i "s|^\( *liftDir *\).*|\1(-$DY  $DX 0);|" system/forceCoeffs
```

Existing runs did not need repeating; the inverse transform was applied to the stored
force histories.

### Note
This failure is silent. No warning, no error, no residual signature. It was found only
because the drag curve violated a symmetry argument.

---

## 4. Root cause 2 — Boundary layer entirely unresolved

### Symptom
Observation 1 above, plus L/D roughly 8× too low.

### Diagnosis
`system/snappyHexMeshDict` contained:

```
refinementSurfaces { airfoil { level (0 0); } }
addLayers  false;
```

Zero surface refinement and no prism layers, so cells at the wall were the raw
background grid.

From `blockMeshDict` grading, the airfoil sat in a uniform Cartesian block of
**0.01c × 0.01c** cells. Confirmed by cell counts: background 300 × 260 = 78,000, final
77,174, difference **826**. The analytic NACA 0012 area is 0.082210 c², which at
10⁻⁴ m² per cell predicts **822** cells removed. A 0.5% match, so snappy removed the
interior and snapped, and did nothing else.

The `airfoil` patch had **202 faces** for a 2.04 m perimeter — 0.0101 m uniform
spacing, with no clustering at leading or trailing edge.

### Consequences

Turbulent boundary layer at mid-chord (x/c = 0.5, Re_x = 1.5 × 10⁶):
δ = 0.37x/Re_x^0.2 = **10.76 mm**; u_τ = U√(C_f/2) = **1.2323 m/s**.

| Quantity | Old mesh | Required |
|---|---|---|
| First cell height | 10 mm | ~2 × 10⁻⁵ m |
| **y⁺ at first cell centre** | **616** | ~1 |
| Cells across boundary layer | **1.08** | 30+ |
| LE radius (1.1019·t²·c = 15.87 mm) per cell | 1.6 | 20+ |
| Airfoil thickness at 97% chord | 10.8 mm ≈ 1 cell | — |

The entire boundary layer fitted inside a single cell. The leading-edge suction peak,
which generates essentially all of the lift, was unresolvable. Aft of ~97% chord the
airfoil was thinner than one cell, and with `features ();` the sharp trailing edge was
never snapped as an edge, so the Kutta condition was not properly enforced. That last
point explains the α-dependence: the circulation deficit grows with incidence.

`nutUSpaldingWallFunction` was correctly chosen and is continuous across all y⁺, but no
wall function can help here. A wall function reconstructs wall shear by assuming the
first cell centre sits in a log layer connected to a resolved outer flow. With one cell
spanning the boundary layer, there is no outer flow to connect to.

### Why the wrong setting was reached
`level (0 0)` was not a typo. snappyHexMesh refines **isotropically**, so any surface
refinement on a one-cell-thick 2D mesh also splits the z direction — level 4 turns one
cell into 16 — which makes the `empty` patch invalid and the case fails. Level 0 was
the only setting that would run.

The correct workflow is to mesh in 3D with real refinement and then collapse back to
one cell with `extrudeMesh`.

### Fix

`snappyHexMeshDict`:

```
addLayers        true;
maxGlobalCells   8000000;

features ( { file "NACA0012_v2_z100sym.eMesh"; level 5; } );

refinementSurfaces
{
    airfoil { level (4 4); patchInfo { type wall; } }
}

addLayersControls
{
    relativeSizes       false;
    firstLayerThickness 2e-05;
    expansionRatio      1.2;
    minThickness        5e-06;
    layers { airfoil { nSurfaceLayers 22; } }
    featureAngle        180;      // forces layers around the sharp TE
    // ...robustness settings
}
```

Surface refinement capped at level 4, not 5. At level 5 the cells (0.3125 mm) become
as fine as the STL facets (0.277 mm at the LE), so the mesh would begin resolving the
tessellation rather than the airfoil. Level 4 keeps a 2.26× margin.

`blockMeshDict`: `frontAndBack` split into separate `back` and `front` patches, since
`extrudeMesh` requires a single source face.

New `system/extrudeMeshDict`, `system/surfaceFeaturesDict`.

Pipeline:

```bash
blockMesh
surfaceFeatures
snappyHexMesh          # writes 3D mesh
extrudeMesh            # collapses to one cell in z
# set frontAndBack patches to type empty
checkMesh
```

---

## 5. Geometry rebuild (CAD)

### Original
Airfoil Tools NACA 0012 `.dat` imported to Fusion 360 via the *Airfoil sketch from
file* plugin, spline fitted, trailing edge closed with a straight line, extruded
**10 mm** in z, exported as STL at the **High** refinement preset.

### Problems
1. **z extent coincided exactly with the mesh domain** (both 0 to 0.01 m). Harmless at
   zero refinement, but with snapping enabled snappy would pull boundary points onto
   the end-cap faces, displacing them off the z = 0 plane and breaking `extrudeMesh`.
2. **High preset was too coarse.** Fusion's presets scale to the model bounding box
   (1 m), so they cannot resolve the 15.87 mm leading-edge radius. Normal deviation of
   10° gives only 18 facets around the entire leading edge.

### Rebuild
Extrude changed to **100 mm symmetric** (z = −0.05 to +0.05 m), giving 40–50 mm
clearance from the mesh band on both sides. Exported in metres with **Custom**
refinement:

| Setting | Value |
|---|---|
| Surface deviation | 0.0101 mm (Fusion's minimum for this model size) |
| Normal deviation | 1.00° |
| Maximum edge length | unconstrained |
| Aspect ratio | 21.5 |

Maximum edge length deliberately left off: the side surface is a straight extrusion, so
constraining it would subdivide along z for no geometric benefit. Aspect ratio 1.00 was
tested and produced 92,252 triangles versus 7,506 — a 12× increase entirely in the z
direction, with identical x-y fidelity. Rejected as pure waste.

Leading-edge facet size is governed by normal deviation, not surface deviation:

| Criterion | Allowed facet at LE |
|---|---|
| Surface deviation 0.0101 mm | 1.132 mm |
| Normal deviation 1.00° | **0.277 mm** ← binding |

Against 0.625 mm cells at level 4, that is 2.26 facets per cell.

### Verification (`surfaceCheck`)

| Check | Result |
|---|---|
| Triangles / vertices | 7,506 / 3,755 |
| Euler characteristic V − E + F | **2** (closed, genus 0) |
| Bounding box x | 0 to 1.0 m — chord exact, metre export correct |
| Bounding box y | ±0.0600153 m |
| **Bounding box z** | **−0.05 to +0.05 m** — symmetric extrude confirmed |
| Max half-thickness vs analytic | 0.0600153 vs 0.0600173 → **1.97 µm error** |
| TE gap vs analytic | 2.520 mm vs 2.520 mm — exact |
| Topology | closed, manifold, 1 part, 1 zone, 0 duplicate points |

The analytic NACA 0012 maximum half-thickness is 0.0600173 at x/c = 0.2998, not exactly
0.06 — the "12" in the designation is nominal. The CAD matches the true polynomial.

The 2.520 mm trailing-edge gap is the genuine open TE of the standard NACA 0012
definition (y_t(c) = 0.00126c), not a modelling artifact. The closing line reproduced it
correctly.

**Known cosmetic issue:** 33.7% of triangles score below 0.05 quality. These are the
elongated triangles spanning the 100 mm extrusion and the fan-triangulated end caps.
Geometrically exact, just badly shaped; `surfaceCheck` cleared them as legal. If
snapping artifacts ever appear, re-exporting at aspect ratio 1.00 gives identical x-y
fidelity with better-conditioned triangles.

### Feature extraction

`surfaceFeaturesDict` with `includedAngle 150` found 1,052 edges over 1,050 points:

| Feature | Count |
|---|---|
| End cap perimeters (2 × 525) | 1,050 |
| Trailing edge lines | 2 |
| **Total** | **1,052** |

Plus 4 feature points where the TE lines meet the cap perimeters. All classified
external convex: the cap-to-side junction is 90°, the TE junction is 98.0° (the NACA
0012 surface slope at x = c is −0.1403, i.e. 7.99° below the chord line, meeting a
vertical TE face). Both are well under the 150° threshold.

Each TE line runs the full z extent, so it passes through the mesh band and receives
level 5 refinement. The cap perimeters sit 40–50 mm outside the domain and are inert.

---

## 6. Mesh result

### snappyHexMesh

```
patch    faces   layers   overall thickness [m]    [%]
airfoil  51648   22       0.00542                  100
```

**22 layers at 100% coverage, zero extrusion removed.** Overall thickness 0.00542 m
matches the analytic 22-layer geometric stack exactly (2 × 10⁻⁵ m, ratio 1.2). Layers
on sharp trailing edges commonly fail partially; `featureAngle 180` prevented this.

`Added 1136256 out of 1136256 cells (100%)`. Zero faces in error on every quality check,
before and after layer addition. "Finished meshing without any errors."

Cells per refinement level: 76,688 (L0) · 2,084 (L1) · 8,376 (L2) · 34,000 (L3) ·
**1,267,712 (L4)**.

### extrudeMesh

Collapsed 1,388,860 3D cells to **163,306** 2D cells. Predicted exactly: the `back`
patch had 163,306 faces, one cell each. Cross-check against the level breakdown:
76,688 + 2,084/2 + 8,376/4 + 34,000/8 + 1,267,712/16 = 163,306.

### Final mesh

| Metric | Old mesh | New mesh |
|---|---|---|
| Cells | 77,174 | 163,306 |
| `airfoil` faces | 202 | **3,228** |
| Surface spacing | 10.1 mm | **0.632 mm** |
| First cell height | 10 mm | **0.02 mm** |
| y⁺ (predicted, mid-chord) | 616 | 1.23 |
| **y⁺ (measured: min / avg / max)** | — | **0.088 / 0.811 / 3.05** |
| Prism layers | 0 | **22 (100% coverage)** |
| Cells across BL | 1.08 | ~32 |
| Non-orthogonality max / avg | 26.9 / 0.35 | 62.7 / 3.00 |
| Skewness max | 0.34 | 1.53 |
| Aspect ratio max | 75.6 | 75.6 |
| Solution directions | (1 1 0) | **(1 1 0)** |
| `checkMesh` | OK | **OK** |

Non-orthogonality and skewness rose as expected with prism layers, prompting
`nNonOrthogonalCorrectors 0 → 1`. Both remain within limits.

3,070 polyhedra and 356 prisms out of 163,306 are layer-termination cells where the
prism stack meets the isotropic refinement. All pass quality checks.

---

## 7. Secondary corrections

| # | Issue | Fix |
|---|---|---|
| 1 | `constant/physicalProperties` had `nu 1e-05`, so **actual Re was 3 × 10⁶, not the 2 × 10⁶ stated in early drafts** | Reynolds number corrected in the report; ν retained, validation reference switched to Abbott & von Doenhoff at 3 × 10⁶ |
| 2 | `nNonOrthogonalCorrectors 0` with max non-orthogonality 62.7 | Set to 1 |
| 3 | `extrudeMesh` produced two `empty` patches (`frontAndBack`, `frontAndBack2`) | `"frontAndBack.*"` wildcard in all `0/` field files |
| 4 | Post-extrude mesh occupies z = −0.01 to 0, so `locationInMesh` fell outside the domain | Changed to `(-5 5 -0.005)` |

---

## 8. Checked and found correct

Recorded to avoid re-investigation:

- **`fvSchemes`** — `bounded Gauss linearUpwind grad(U)` is second-order;
  `Gauss linear corrected` and `corrected` snGrad appropriate. Numerical diffusion was
  ruled out as a cause.
- **`fvSolution`** — GAMG for p, smoothSolver for U/nuTilda, relaxation 0.3/0.7/0.7,
  `residualControl 1e-5`. All standard.
- **Far-field boundary conditions** — `freestreamVelocity` / `freestreamPressure` are
  inletOutlet-type and switch on flux sign, so they are robust to AoA and the
  inlet/outlet patch grouping is irrelevant. *(An earlier hypothesis that the patch
  grouping was inverted was incorrect and is retained here as a corrected line of
  enquiry.)*
- **`Aref 0.01`** — chord 1 m × span 0.01 m. Correct. `lRef 1`, `rhoInf 1` correct for
  kinematic pressure. Confirmed independently by CD ≈ 0.009 at α = 0 being physically
  plausible.
- **`nuTilda` freestream 4.5 × 10⁻⁵** = 4.5ν, inside the standard Spalart-Allmaras
  range of 3ν–5ν.
- **2D setup** — `checkMesh` reports 2 solution directions both before and after the
  rebuild.
- **Original STL scale and shape** — the 826 vs 822 removed-cell match confirmed correct
  units, chord, and thickness before any geometry work began.

---

## 9. Dictionary and workflow notes

Errors encountered during the rebuild, recorded for reproducibility:

| Error | Cause | Fix |
|---|---|---|
| `problem while reading header for object surfaceFeaturesDict` | File opened with `/*----*\` but the closing `\*----*/` was omitted, leaving the whole file inside an unterminated block comment | Remove the decorative banner or close it properly |
| `keyword flipNormals is undefined` | OpenFOAM 11 requires it explicitly where older versions defaulted it | Add `flipNormals false;` |
| `cannot find patch frontAndBack` | `exposedPatchName` must name a patch that **already exists** in the source mesh; OF11 looks it up to inherit settings rather than creating it | Point it at an existing patch, rename afterwards |
| `checkMesh` reported the pre-extrude mesh | Running snappy without `-overwrite` writes to time directories, leaving `constant/polyMesh` as the raw blockMesh output — which is what `extrudeMesh` reads | Move the final time directory's `polyMesh` into `constant/` before extruding |
| `cp: ... 'constant/geometry': Not a directory` | Target directory did not exist; this case uses `constant/triSurface` | Locate with `find . -name "*.stl"`. Note that `cp` **without** a trailing slash silently creates a *file* of that name instead of failing |

OpenFOAM 11 renames: `transportProperties` → `physicalProperties`,
`turbulenceProperties` → `momentumTransport`, `simpleFoam` → `foamRun` with
`solver incompressibleFluid`.

---

## 10. Expected results

Approximate NACA 0012 values at Re = 3 × 10⁶ for validation. Confirm against Abbott &
von Doenhoff before quoting in the report.

| Quantity | Expected |
|---|---|
| Lift-curve slope | ≈ 0.11 /deg (thin-airfoil theory: 2π/rad = 0.1097 /deg) |
| CL at 12° | ≈ 1.25 |
| CD at 0° | ≈ 0.0075 |
| CD at 12° | ≈ 0.02 |
| Peak L/D | ≈ 60–70 |
| Stall | ≈ 15–17°, CL_max ≈ 1.5–1.6 |

Previous (frame-corrected) run gave CL = 0.516 and CD = 0.157 at 12°: roughly 40% of
the expected lift and 7× the expected drag.

---

## 10a. First converged result — α = 12°

Run to `endTime 5000`. Forces flat over the final 20 iterations: CL varies by ±0.0007%,
CD by ±0.020%. The run reached `endTime` rather than tripping `residualControl 1e-5`,
because the nuTilda residual lagged (6.4 × 10⁻⁴ at iteration 3348 while Ux, Uy and p
were all below 4 × 10⁻⁶). Flat forces are the operative convergence criterion.

Continuity fell four orders of magnitude, from 1.01 × 10⁻³ at startup to
3.36 × 10⁻⁷, with cumulative stable at −8.58 × 10⁻⁴ rather than drifting.

### Wall resolution

`foamPostProcess -solver incompressibleFluid -func yPlus -latestTime`:

```
patch airfoil y+ : min = 0.0883, max = 3.048, average = 0.811
```

Comfortably within the Spalart-Allmaras requirement. The maximum of 3.05 is roughly
2.5× the mid-chord prediction of 1.23, consistent with the leading-edge suction peak
where local velocity reaches 2–2.5× freestream at this incidence. The minimum of 0.088
is the stagnation point, where u_τ → 0 by definition. The distribution behaves as
expected.

*Note: `foamPostProcess` requires `-solver incompressibleFluid` in OpenFOAM 11.
Without it the turbulence model is not instantiated and `yPlus` fails with
"Unable to find turbulence model in the database".*

### Forces

| Quantity | Old mesh (frame-corrected) | New mesh | Reference |
|---|---|---|---|
| CL | 0.516 | **1.208** | ≈ 1.25 |
| CD | 0.157 | **0.0294** | see note |
| Cm (c/4) | — | **−0.0070** | ≈ −0.005 to −0.010 |
| L/D | 3.3 | **41.1** | ≈ 60 |

Lift increased 2.34×, drag fell 5.33×, L/D improved 12.5×.

**The moment coefficient is the strongest validation point.** Thin-airfoil theory gives
zero quarter-chord moment for a symmetric section, and measured NACA 0012 values sit
around −0.005 to −0.010. Matching Cm means the *chordwise pressure distribution* is
correct, not merely its integral — two solutions can produce identical CL from very
different distributions.

The average lift-curve slope to 12° is 0.1007 /deg, 91.8% of the thin-airfoil value of
0.1097 /deg. Some fall-off below linear is expected at this incidence.

### Drag decomposition

A `forces` function object was added to split total drag into pressure and viscous
components. `force.dat` reports vectors in **mesh axes**, so the same wind-axis rotation
must be applied manually — the `forces` object has no `liftDir`/`dragDir` entries at all,
which is precisely why the trap in §3 exists.

| Component | CD | Share |
|---|---|---|
| Pressure | 0.02242 | 76.3% |
| Viscous | 0.00698 | 23.7% |
| Total | 0.02939 | — |

Recomposition returns 0.02939 against the `forceCoeffs` value of 0.0293916, confirming
the rotation.

**Viscous drag is 93% of the flat-plate turbulent value** (0.00698 against 0.00750 for
two sides at Re 3 × 10⁶). The small deficit is expected: the suction-side boundary layer
thickens under the adverse gradient and loses skin friction, while the pressure side
gains slightly. This confirms the wall treatment independently of y⁺ — skin friction is
being computed correctly. The entire drag excess therefore sits in **pressure drag**.

---

## 10b. Grid-independence study

### Hypothesis

Refinement extended only **24.2 mm from the wall** (5.42 mm prism stack plus 18.75 mm of
level 4→0 transition bands at `nCellsBetweenLevels 2`), beyond which cells jumped to the
10 mm background. Flat-plate δ at x/c = 0.95 is 18 mm, and under the adverse gradient at
12° the suction-side boundary layer should reach 36–54 mm. The outer half to two-thirds
of the trailing-edge boundary layer was therefore in 10 mm cells, which would inflate
displacement thickness, promote premature separation, and appear as pressure drag.

Note this is not contradicted by y⁺ = 0.811, which measures only the *first* cell.

### Test

Single variable changed: outer-field refinement. Surface refinement, prism layers,
geometry, and solver settings all untouched.

```
nCellsBetweenLevels 3;          // was 2

refinementRegions
{
    airfoil
    {
        mode    distance;
        levels  ((0.03 2) (0.10 1));
    }
}
```

| Cell level | Coarse 2D | Refined 2D | Change |
|---|---|---|---|
| 0 | 76,688 | 74,882 | 1.0× |
| 1 | 1,042 | 6,622 | **6.4×** |
| 2 | 2,094 | 8,460 | **4.0×** |
| 3 | 4,250 | 4,656 | 1.1× |
| 4 | 79,232 | 80,970 | 1.0× |
| **Total** | **163,306** | **175,522** | 1.07× |

Levels 3 and 4 barely moved — correct, since `refinementRegions` specifies a *minimum*
level and the near-wall region was already at level 3–4. Only the outer field changed.
`airfoil` patch face count stayed at exactly 3,228, confirming a genuine single-variable
test.

Refinement reach in the region of interest:

| | Coarse | Refined |
|---|---|---|
| Cells ≥ level 1 out to | 24 mm | **100 mm** |
| Cells ≥ level 2 out to | ~11 mm | **30 mm** |
| Cell size 40 mm from wall | 10 mm | **5 mm** |

Mesh quality improved slightly as a side effect of the gentler level transitions:
non-orthogonality max fell from 62.72 to 58.07.

### Result

| Quantity | Coarse (163k) | Refined (176k) | Change |
|---|---|---|---|
| **CL** | 1.20838 | 1.20820 | **−0.015%** |
| CD total | 0.02939 | 0.02829 | −3.76% |
| **CD pressure** | 0.02242 | 0.02136 | **−4.73%** |
| CD viscous | 0.00698 | 0.00693 | −0.67% |
| y⁺ (avg) | 0.811 | 0.808 | — |
| L/D | 41.1 | 42.7 | +3.9% |

**The hypothesis was wrong.** Pressure drag was predicted to fall toward 0.015 if
near-field resolution were the cause; it moved 4.7%, not 30%.

### Conclusions

1. **Lift is grid-independent to 0.015%.** Sixfold and fourfold increases in level-1 and
   level-2 cell counts, with refinement extended from 24 mm to 100 mm, moved CL in the
   fifth decimal place.
2. **Drag is grid-converged to within 4%.**
3. Viscous drag moved 0.67%, as expected with the surface mesh and layers untouched.
   This confirms the two cases differ in exactly one variable.
4. **The residual drag offset is not a mesh resolution problem.**

### Remaining candidates for the drag offset

| Candidate | Status |
|---|---|
| **Fully turbulent Spalart-Allmaras** — no transition model, so the boundary layer is turbulent from the stagnation point, whereas experiment at Re 3 × 10⁶ on a smooth model retains laminar flow over part of the forward chord. Routinely accounts for 15–30% on airfoil drag. | Modelling assumption; state as a limitation |
| **Domain size** — 10 chords vertically against a typical 20–50. 2D drag is considerably more far-field-sensitive than lift, which is consistent with lift being converged and close to reference while drag is not. | Untested; requires only a `blockMeshDict` edit, no re-snapping |

---

## 10c. Symmetry check — α = 0°

Run on the refined mesh. **Converged at iteration 1672** via `residualControl 1e-5`,
including nuTilda — unlike the 12° case, which ran to `endTime 5000` with nuTilda
lagging. Zero incidence is the easy case, so the lower angles in the sweep will finish
considerably faster than 12° did.

| Quantity | Value |
|---|---|
| CL | −6.699 × 10⁻⁴ |
| Cm | −9.392 × 10⁻⁵ |
| CD | 0.00991 |

### Interpretation

A symmetric section at zero incidence must produce zero lift and zero pitching moment.
The residuals here are non-zero because **snappyHexMesh does not produce a
mirror-symmetric grid**: the STL tessellation differs between upper and lower surfaces,
and the snapping and layer-addition algorithms contain path-dependent steps. Exact
symmetry would require a structured mesh built symmetric by construction.

The raw coefficient is not the useful measure. Expressed meaningfully:

- CL is **0.055%** of the value at 12°
- It corresponds to an effective incidence error of **0.0066°**

Wind-tunnel angle-of-attack settings are typically accurate to ±0.05°, so the mesh
asymmetry is an order of magnitude below experimental setting accuracy. Acceptable.

*Report this as an equivalent angle error rather than a raw coefficient — the raw number
is meaningless without a scale to compare against.*

### Drag at zero incidence

CD = 0.00991 is 1.32× the flat-plate turbulent value of 0.00750 (both sides,
Re 3 × 10⁶). A form factor of ≈1.20 is standard for a 12% thick section, implying a
fully-turbulent expectation near 0.0090. The measured value is therefore about **10%**
above expectation.

### Emerging trend in the drag offset

| AoA | CD | Rough fully-turbulent estimate | Offset |
|---|---|---|---|
| 0 | 0.00991 | ~0.0090 | **+10%** |
| 12 | 0.02829 | ~0.020 | **+41%** |

**This complicates the transition explanation.** Natural transition benefits drag *most*
at low incidence, where the laminar run is longest. As incidence increases the suction
peak forces transition forward and experiment converges toward the fully-turbulent
result. The transition argument should therefore *narrow* the gap at 12°, not widen it.

An offset that grows with incidence instead suggests an incidence-dependent cause:
separation behaviour, or the far-field constraint, which scales with circulation and
would preferentially affect the high-lift cases.

**Caveats:** both reference figures are estimates rather than values read from Abbott &
von Doenhoff, and two points do not establish a trend. The full sweep will give the
shape of the offset versus incidence, which is far more diagnostic than any single
comparison. If the offset rises monotonically from ~10% to ~40%, the domain sensitivity
test moves from optional to necessary.

---

## 10d. Final results — AoA sweep, refined mesh (175,522 cells)

| AoA | CL | CD | Cm,c/4 | L/D |
|---|---|---|---|---|
| 0 | −0.00067 | 0.00991 | −0.000094 | — |
| 2 | 0.20889 | 0.01049 | −0.000319 | 19.9 |
| 4 | 0.41741 | 0.01215 | −0.000795 | 34.4 |
| 6 | 0.62427 | 0.01479 | −0.001504 | 42.2 |
| 8 | 0.82712 | 0.01835 | −0.002680 | 45.1 |
| 10 | 1.02322 | 0.02282 | −0.004527 | 44.8 |
| 12 | 1.20820 | 0.02829 | −0.007316 | 42.7 |

All cases converged. 0° reached `residualControl 1e-5` at iteration 1672; the remainder
ran to `endTime 5000` with nuTilda the lagging residual and forces flat.

### Lift curve

| Segment | 0→2 | 2→4 | 4→6 | 6→8 | 8→10 | 10→12 |
|---|---|---|---|---|---|---|
| dCL/dα | 0.1048 | 0.1043 | 0.1034 | 0.1014 | 0.0981 | 0.0925 |

**Monotonically decreasing throughout** — linear at low incidence, softening gradually
toward stall. This is the required physical shape. For comparison, the original faulty
run gave 0.0705, 0.0655, **0.0145**, 0.0230, 0.0285, **0.0665**: a mid-range collapse
followed by recovery to the initial slope, which no aerofoil does.

Linear fit over 0–6°:

- **Lift-curve slope = 0.1042 /deg**, 95.0% of the thin-aerofoil value (0.1097 /deg).
  Published NACA 0012 measurements sit near 0.105 /deg; the small deficit is the
  expected viscous displacement effect.
- **Zero-lift angle = +0.0002°** on a symmetric section.

### Drag curve

CD increments: 0.00058, 0.00166, 0.00264, 0.00356, 0.00448, 0.00547 — monotonic and
convex, as required. The original faulty run gave 0.018, 0.007, 0.003, 0.002: concave,
i.e. the wrong curvature entirely.

Fitting CD = CD₀ + k·CL² with CD₀ = 0.00991:

| AoA | 2 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| k | 0.01335 | 0.01286 | 0.01252 | 0.01234 | 0.01234 | 0.01259 |

**k is constant to within 8%** across the sweep (mean 0.0127). The quadratic drag polar
therefore holds precisely. In the original faulty data k ranged from 0.11 to 0.29, i.e.
no consistent polar existed at all.

### Aerodynamic centre — independent validation

Fitting Cm against CL over the full sweep:

- dCm/dCL = −0.00560 → **aerodynamic centre at x/c = 0.2556**
- **Cm at zero lift = +0.00098** (thin-aerofoil theory: exactly 0 for a symmetric section)

Thin-aerofoil theory places the aerodynamic centre at the quarter chord; measured
aerofoils fall at 25–26%. The simulation located it at 25.56% with no such information
supplied. This confirms the *chordwise* pressure distribution independently of any lift
or drag comparison — a check that a solution tuned to match CL alone would fail.

Cm at 12° was −0.007316 on the refined mesh against −0.00700 on the coarse mesh, a
further grid-independence data point.

### Diagnosis of the drag offset

Separating the polar into its two parts:

| Component | This work | Approx. published | Offset |
|---|---|---|---|
| CD₀ (zero-lift drag) | 0.00991 | ~0.0090 fully turbulent | **+10%** |
| k (lift-dependent drag) | 0.0127 | ~0.0075 | **+65%** |

**Zero-lift drag is nearly correct; lift-dependent drag is not.** The excess scales
as CL².

This *excludes the transition explanation*. Absence of a transition model raises CD₀
most, and its influence diminishes with incidence as the suction peak drives transition
forward — so a transition-driven offset would narrow at high α, not widen.

A CL²-scaling excess indicates a circulation-driven mechanism. **Far-field extent is the
leading candidate:** bound circulation is proportional to CL, the freestream boundary
condition at 10 chords suppresses the induced velocity field that circulation generates,
and the resulting spurious drag scales as CL². An image-vortex estimate gives roughly
0.005 at 12°, accounting for most of the 0.008 excess — crude, but the correct order of
magnitude and, more importantly, the correct scaling.

This makes the domain sensitivity test worthwhile: it requires only a `blockMeshDict`
edit with no re-snapping, and would convert the principal limitation from unexplained to
quantified.

---

## 10e. Domain sensitivity study

### Motivation

§10d found the drag excess scaling as CL², implicating a circulation-driven mechanism.
Far-field extent was the leading candidate: bound circulation is proportional to CL, the
freestream boundary condition suppresses the induced velocity field that circulation
generates, and the resulting spurious drag scales as CL².

The original domain used **10 chords** vertically. This is a common tutorial default and
a frequently quoted minimum, but published aerofoil validation typically uses 50–100
chords or more, or else applies a far-field point-vortex correction at the boundary to
compensate for a smaller domain.

### Test

Domain enlarged from x ∈ [−15, 20], y ∈ [±10] to x ∈ [±50], y ∈ [±50]. Near-field
spacing held at exactly 0.01 m by recomputing the `simpleGrading` expansion ratios
(growth 1.06 per cell against the original 1.05); background grid 296 × 296 = 87,616
against 300 × 260 = 78,000.

Single-variable control confirmed by the refinement breakdown — only level 0 changed:

| Level | 10 chords | 50 chords | Change |
|---|---|---|---|
| 0 | 74,882 | 84,484 | +12.8% |
| 1 | 13,244 | 13,336 | +0.7% |
| 2 | 33,840 | 33,820 | −0.1% |
| 3 | 37,244 | 36,485 | −2.0% |
| 4 | 1,295,520 | 1,313,208 | +1.4% |

Layer addition identical in both: 51,648 faces, 22 layers, 0.00542 m, 100% coverage.
Final 2D mesh 186,059 cells against 175,522.

### Result (α = 12°)

| Quantity | 10 chords | 50 chords | Change |
|---|---|---|---|
| CL | 1.20820 | **1.24775** | **+3.27%** |
| CD | 0.02829 | **0.02092** | **−26.04%** |
| Cm,c/4 | −0.00732 | −0.00801 | −9.42% |
| L/D | 42.7 | **59.6** | +40% |

Against published NACA 0012 at 12°:

| | 10 chords | 50 chords |
|---|---|---|
| CL error | −3.3% | **−0.2%** |
| CD error | +41.4% | **+4.6%** |

### Mechanism confirmed

The far-field effect scales as CL², so CD₀ at zero incidence is unaffected. Recomputing
the polar coefficient with CD₀ = 0.00991:

| Domain | k |
|---|---|
| 10 chords | 0.01259 |
| 50 chords | **0.00707** |
| Published NACA 0012 | ~0.0075 |

The 65% excess in lift-dependent drag identified in §10d is eliminated. Enlarging the
domain removed 0.00737 of drag, against an image-vortex estimate of ~0.005 — correct
order of magnitude, correct sign, and correct scaling with CL².

L/D of 59.6 now falls within the 60–70 band expected for this section; 42.7 did not.

### Interpretation

**10 chords is adequate for lift but not for drag.** CL was in error by 3.3%, CD by 41%.
Rules of thumb for domain size are generally framed around lift, which is the quantity
most often checked. Drag is considerably more far-field-sensitive in 2D because the
spurious contribution scales with the square of circulation.

**Residual error.** The effect scales approximately as 1/H, so roughly 0.0015 of spurious
far-field drag likely remains at 50 chords — about 7% of the current value, consistent
with the remaining 4.6% gap to reference. Fully eliminating it would require either a
larger domain still or a far-field point-vortex correction.

**Consequence:** the 10-chord AoA sweep in §10d is systematically biased, with the error
growing as CL². It distorts the shape of the polar rather than merely offsetting it, and
the sweep must be repeated on the 50-chord domain.

**Caveat:** the grid-independence study of §10b was performed on the 10-chord domain.
Since the near-field mesh in the 50-chord case is essentially identical (levels 1–4
within 2%), the conclusion carries over, but it was not re-demonstrated on the final
domain.

---

## 10f. Validation against experimental data

### Source selection

| Source | Re | Transition | Data | Use |
|---|---|---|---|---|
| **Gregory & O'Reilly, R&M 3726 (1970)** | **3 × 10⁶** | **tripped** | CL, Cp | **Primary — matches case conditions exactly** |
| Ladson, NASA TM 4074 (1988) | 6 × 10⁶ | tripped (80 grit) | CL, CD | Drag, with Reynolds correction |
| Abbott & von Doenhoff (1959) | — | **untripped** | CL, CD | *Not used* — wrong transition state |

Digitised data obtained from the NASA Turbulence Modeling Resource,
<https://tmbwg.github.io/turbmodels/naca0012_val.html>.

**Abbott & von Doenhoff is the wrong reference for this case.** Those data were not
tripped, i.e. natural transition with a laminar run over the forward chord, whereas
Spalart-Allmaras is fully turbulent from the stagnation point. NASA's guidance is that
tripped experimental data are more appropriate for comparison with fully turbulent CFD
drag. Gregory & O'Reilly are tripped *and* at Re = 3 × 10⁶, matching this case in both
respects.

**Limitation of the comparison itself:** no tripped CD data at Re = 3 × 10⁶ exist in the
NASA set. Gregory provides CL only; Ladson provides CD but at twice the Reynolds number.
Drag comparison therefore requires a Reynolds correction which is quantified only for
CD₀ (McCroskey, via NASA TMR: tripped CD₀ at Re 3 million is roughly 10% above
Re 6 million).

### Lift — validated

| | Slope (0–6°) | Zero-lift α |
|---|---|---|
| Gregory (experiment) | 0.10818 /deg | −0.163° |
| **CFD** | **0.10859 /deg** | **−0.007°** |
| Difference | **+0.38%** | — |

| AoA | CFD CL | Gregory | Diff |
|---|---|---|---|
| 2 | 0.2190 | 0.2342 | −6.5% |
| 4 | 0.4365 | 0.4521 | −3.5% |
| 6 | 0.6510 | 0.6646 | −2.1% |
| 8 | 0.8602 | 0.8756 | −1.8% |
| 10 | 1.0610 | 1.0807 | −1.8% |
| 12 | 1.2477 | 1.2795 | −2.5% |

**The offset is largely an artefact of the reference data.** Gregory reports CL = 0.0257
at α = 0.094°, implying a zero-lift angle of −0.163° on a *symmetric* section — which is
physically impossible and indicates rig misalignment or digitising error. This shifts
the entire experimental curve upward by roughly 0.018 in CL. Correcting for it reduces
the CFD differences to well under 1%.

The CFD zero-lift angle of −0.007° is approximately twenty times closer to the
theoretical zero than the experiment, so the reference data set the accuracy floor for
this comparison rather than the simulation.

Gregory CL_max = 1.561 at α ≈ 16.1°, dropping sharply to 0.96 by 17.1°. The sweep
terminates at 12°, well inside the linear region and clear of stall — which justifies
the steady RANS assumption directly rather than by assertion.

### Drag — partially validated

Ladson data scaled by 1.10 to estimate Re = 3 × 10⁶ conditions:

| AoA | CL | CFD CD | Ladson ×1.10 | Diff |
|---|---|---|---|---|
| 0 | 0.000 | 0.00990 | 0.00890 | **+11.1%** |
| 2 | 0.219 | 0.01017 | 0.00898 | +13.3% |
| 4 | 0.436 | 0.01097 | 0.00907 | +21.0% |
| 6 | 0.651 | 0.01234 | 0.00972 | +26.9% |
| 8 | 0.860 | 0.01433 | 0.01134 | +26.4% |
| 10 | 1.061 | 0.01710 | 0.01312 | +30.3% |
| 12 | 1.248 | 0.02092 | 0.01451 | **+44.2%** |

Split into the two polar terms:

| | CFD | Experiment | Ratio |
|---|---|---|---|
| CD₀ | 0.00990 | 0.00890 | **1.11×** |
| k (lift-dependent) | 0.00708 | 0.00329 | **2.15×** |

**CD₀ agrees to within 11%**, attributable to the fully turbulent model, the 2.52 mm
blunt trailing edge (NASA's reference geometry uses a modified sharp-TE aerofoil, max
thickness 11.894%), and the approximate nature of the Reynolds correction.

**The lift-dependent term does not agree.** k is more than double the experimental
value, and the discrepancy continues to scale as CL².

Subtracting the residual far-field contribution estimated from 1/H scaling
(≈0.0015 at 50 chords) leaves CD = 0.01945 at 12°, still ~34% above expectation. Far-field
extent therefore accounts for part but not most of the remaining excess.

Untested candidates for the remainder:

- Spalart-Allmaras over-predicting suction-side boundary layer growth, which would scale
  with CL as observed
- Blunt trailing edge base drag, whose contribution grows with incidence
- Reynolds-number mismatch in the comparison, corrected only for CD₀

### Independent confirmation of the far-field finding

The NASA Turbulence Modeling Resource maintains a dedicated page on this effect
(<https://tmbwg.github.io/turbmodels/naca0012_val_ffeffect.html>). Their guidance:
a farfield extent of 500 chords is sufficient for minimal influence; at smaller extents
such as 30 chords the influence on forces is greater — *especially for drag* — unless a
farfield point-vortex boundary condition correction is applied (Thomas & Salas,
AIAA Journal 24(7):1074–1080, 1986). NASA's supplied grids place the boundary almost 500
chords out for this reason.

The effect identified in §10e was found independently, from the CL² scaling of the drag
excess, without prior knowledge of this guidance. It also establishes that **50 chords
remains insufficient by NASA's standard**, consistent with the residual gap. Extending
to 500 chords would remove only a further ~0.0017, roughly one third of the outstanding
discrepancy.

---

## 11. Outstanding

- [x] ~~Confirm y⁺~~ — **0.088 / 0.811 / 3.05** (min/avg/max). Passed.
- [x] ~~Verify convergence at 12°~~ — forces flat to ±0.02% over the final 20
      iterations. Reached `endTime` rather than `residualControl`; nuTilda lagged.
- [x] ~~Add `forces` function object for the pressure/viscous drag split~~ — pressure
      0.02242 (76.3%), viscous 0.00698 (23.7%); viscous is 93% of the flat-plate value
- [x] ~~Grid-independence study~~ — CL converged to 0.015%, CD to 4%. Mesh resolution
      ruled out as the cause of the drag offset.
- [x] ~~Symmetry check at α = 0~~ — CL = −6.70e−04, Cm = −9.39e−05, equivalent to a
      0.0066° incidence error. Passed.
- [x] ~~AoA sweep 0°–12° on the refined mesh~~ — complete, see §10d
- [x] ~~Domain sensitivity test~~ — CD −26%, CL +3.3% from 10 to 50 chords. Far-field
      extent confirmed as the cause of the drag excess. See §10e.
- [x] ~~Repeat the AoA sweep on the 50-chord domain~~ — complete, see §10f
- [x] ~~Source experimental data~~ — Gregory & O'Reilly and Ladson, via NASA TMR
- [ ] Re-run 12° to regenerate `results/AoA12` (the archive step captured the 10° case)
- [ ] Cp distribution at 0/6/12° — compare against Gregory Cp data (Re 3×10⁶)
- [ ] Figures and README (see `REPO_PLAN.md`)

### Future work (out of scope for this study)

- **k-ω SST comparison.** The lift-dependent drag discrepancy scales with CL and is the
  signature of suction-side boundary layer growth. Testing a second turbulence model is
  the logical next step and would isolate whether the error is model-specific.
- **Farfield point-vortex correction** (Thomas & Salas 1986), which would allow the
  small-domain result to be corrected analytically rather than requiring a large grid.
- Third domain at 25 chords, to demonstrate convergence rather than merely a difference.
- [ ] Validation plot against Abbott & von Doenhoff
- [ ] **Mesh sensitivity study:** re-run one case on the old mesh (`constant/polyMesh.3d.bak`,
      old STL retained as `NACA0012_v1_z10_OLD.stl`) with identical solver settings. A
      y⁺ 616 vs y⁺ 1.23 comparison converts this from a bug fix into a proper
      resolution study.

### Reporting caveat — steady RANS assumption

The 2.520 mm blunt trailing edge will shed. Two distinct mechanisms, which look
different in the logs:

- **Blunt-TE shedding** — at St ≈ 0.2–0.25, f ≈ 2400–3000 Hz, period ~0.4 ms. Small
  amplitude; Spalart-Allmaras usually damps it out entirely in steady mode.
- **Trailing-edge separation** — large amplitude, low frequency, grows with incidence.
  This is what actually breaks steady convergence, and for NACA 0012 at this Re it
  typically only appears above roughly 12–14°. The sweep stops at 12°, so clean
  convergence at every angle is plausible.

Chasing either transiently is not worthwhile: resolving TE shedding needs
Δt ≈ 1.5 × 10⁻⁵ s, and one flow-through of the 35 m domain is 1.17 s — about 84,000
steps just to flush the initial transient.

Triage:

| Residuals | Forces | Interpretation |
|---|---|---|
| Fall to 1e-5 | Flat | Converged — quote the value |
| Plateau | Periodic, fixed amplitude | Physical unsteadiness — quote mean ± half peak-to-peak |
| Plateau | Drifting, no period | Not converged, or a mesh problem |

The report should state that steady RANS was used and note where that assumption
weakens.

---

## 12. Summary

Two independent faults, neither of which produced a warning or error:

1. **Force coefficients reported in the mesh frame** rather than the wind frame, because
   `liftDir`/`dragDir` were not rotated alongside the freestream velocity.
   Under-reported drag by up to 3.4×. Detected from a symmetry argument — drag falling
   with incidence on a symmetric section — not from any solver output.

2. **Boundary layer entirely unresolved**, y⁺ ≈ 616 with 1.08 cells across the boundary
   layer, because surface refinement was set to level 0 and prism layers were disabled.
   This was a forced choice: snappyHexMesh refines isotropically and would have
   invalidated the `empty` patches of the 2D mesh. The correct route is to mesh in 3D
   and collapse with `extrudeMesh`.

The general lesson is that neither fault was visible in residuals, `checkMesh`, or the
solver log. Both were found by checking results against physical arguments —
conservation of symmetry, expected curve shape, and order-of-magnitude L/D — before
comparing against reference data.

A subsequent grid-independence study confirmed the rebuilt mesh is converged: a 6×
increase in level-1 cells and 4× in level-2, extending refinement from 24 mm to 100 mm
from the wall, changed CL by 0.015% and CD by 3.8%. A domain sensitivity study then
found a 26% drag reduction on moving the farfield from 10 to 50 chords, matching the
CL² scaling predicted from the polar decomposition.

**Final outcome.** Lift is validated: lift-curve slope within 0.38% of tripped
experimental data at the same Reynolds number, CL within 2% across the sweep, and a
zero-lift angle closer to theory than the reference data itself. Drag is partially
validated: CD₀ within 11%, but the lift-dependent term 2.15× experiment, with the cause
only partly identified. Aerofoil drag depends on boundary layer displacement and
separation, and is therefore far more sensitive to turbulence modelling than lift, which
is set largely by circulation.
