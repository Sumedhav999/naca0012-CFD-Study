# NACA 0012 CFD Validation in OpenFOAM

Steady RANS simulation of a NACA 0012 aerofoil at Re = 3×10⁶ in OpenFOAM 11,
built from CAD through post-processing and validated against tripped wind-tunnel
data.

![CL vs alpha](docs/figures/cl_alpha.png)

## Key Results

| Metric | CFD | Reference | Error |
|---|---|---|---|
| CL @ 12° | 1.248 | 1.250 (Gregory) | −0.2% |
| Lift-curve slope | 0.1086 /deg | 0.1082 /deg (Gregory) | +0.38% |
| CL @ 0° | −0.0005 | 0 (symmetric) | ≈0 |
| Maximum L/D | 62.05 (at 10°) | — | — |

## Validation

Compared against Gregory & O'Reilly (lift) and Ladson, NASA TM 4074
(drag, Reynolds-corrected) — both tripped datasets, matched to the fully-turbulent
Spalart-Allmaras setup used here.

- Symmetry check: CL and Cm vanish at α = 0° on a symmetric section
- y⁺: 0.81 average, 3.05 max — measured, not assumed
- Grid independence: doubling refinement changed CL by 0.015%
- Domain sensitivity: 10 vs 50 chords changed CD by 26% — identified and corrected

## Method

| | |
|---|---|
| Aerofoil | NACA 0012, 1 m chord |
| Reynolds number | 3×10⁶ |
| Mesh | blockMesh + snappyHexMesh, ~186k cells, 22 prism layers |
| Solver | OpenFOAM 11, incompressibleFluid, SIMPLE |
| Turbulence | Spalart-Allmaras, fully turbulent |

## Mesh

![mesh overview](docs/figures/mesh_overview.png)
![mesh leading edge](docs/figures/mesh_le.png)

## Flow Visualisation

### Pressure Contours (12° AoA)

![pressure](docs/figures/pressure_12deg.png)

## Limitations

- Lift is well-validated; the lift-dependent drag term is not — still roughly
  2× experiment, candidates narrowed but not resolved
- Steady RANS only; the blunt trailing edge sheds, untested near/at stall
- Fully turbulent — no transition model
- 50-chord domain, short of the ~500-chord guidance for fully domain-independent drag

## Debugging Highlights

Two implementation faults were found through physical-consistency checks rather
than solver diagnostics: force coefficients computed in the wrong reference frame
(caught via a symmetry argument — drag cannot fall with incidence on a symmetric
section) and a boundary layer left completely unresolved by the mesh (caught via
an order-of-magnitude L/D check). [Full log →](docs/DEBUGGING_LOG.md)

## Reproduce

```bash
cd case
./Allrun
```

## References

1. Gregory, N. & O'Reilly, C.L. (1970) *Low-Speed Aerodynamic Characteristics of NACA 0012*, NASA R&M 3726.
2. Ladson, C.L. (1988) *Effects of Independent Variation of Mach and Reynolds Numbers on the Low-Speed Aerodynamic Characteristics of the NACA 0012*, NASA TM 4074.
3. OpenFOAM User Guide, OpenFOAM 11.
