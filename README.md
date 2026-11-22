# Serial D3Q19 LBM: 3D flow around a generic ogive body

This is a **working single-threaded C++17 3D CFD simulator**, not a 3D mesh generator. The solid body is computed **inside the fluid solver** as a voxelized bounce-back obstacle. There is no STL input or external surface-meshing step.

The geometry is a **generic streamlined axisymmetric test body** with tangent-ogive nose, cylindrical section, and linear boattail ending in a flat rear cap. No Mk 82 or other weapon's actual dimensions or special components are modeled.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build --output-on-failure
./build/lbm_3d --steps 500 --every 250 --out output_3d
```

No dependencies beyond a C++17 toolchain and CMake. Serial execution only: no CUDA, MPI, OpenMP, GPU, or Kokkos.

## Model and flow axes

The **longitudinal body axis is Z** and the inlet velocity points in **+Z**. `nx,ny` are the crossflow lattice extents, `nz` is streamwise. For each lattice cell `(x,y,z)`:

```cpp
const double r = profile_.radius_at(double(z));
const double dx = x - (nx - 1.0)/2.0;
const double dy = y - (ny - 1.0)/2.0;
solid[cell] = (r > 0.0 && dx*dx + dy*dy <= r*r);
```

The continuous underlying surface is `x = xc + R(z)*cos(phi)`, `y = yc + R(z)*sin(phi)`. A D3Q19 Cartesian LBM simulation does not need triangles or sampling in `phi` to implement the solid boundary; it needs a fluid/solid mask. Each link meeting a solid voxel uses bounce-back to approximate a stationary no-slip boundary.

`radius_at()` uses a true tangent-ogive formula at the nose, a constant-radius middle, and a linearly narrowing rear section, and stops abruptly at the rear face to make a flat base. All dimensions are **generic lattice-unit proportions**, not measurements of a particular object:

- Maximum radius: `min(nx,ny)/8.5`
- Nose length: 2.8 radii
- Midsection length: 2.0 radii
- Taper length: 2.7 radii
- Rear radius: 0.60 of the maximum radius
- Front tip begins at 16% of the computational Z length

For better boundary resolution, increase all three lattice extents while preserving the relevant ratios.

## Parameters

| Argument | Default | Description |
|---|---:|---|
| `--nx` | 48 | Lattice cells across X |
| `--ny` | 48 | Lattice cells across Y |
| `--nz` | 112 | Lattice cells along Z (flow direction) |
| `--steps` | 500 | Number of time steps |
| `--every` | 250 | Snapshot interval |
| `--u` | 0.045 | Uniform +Z inlet velocity in lattice units |
| `--re` | 35 | Nominal Reynolds number based on maximum diameter |
| `--smagorinsky` | 0 | Experimental, unvalidated sub-grid coefficient |
| `--out` | output_3d | VTK volume output directory |

The model uses D3Q19 BGK collision and pull streaming, no-slip halfway bounce-back at solid links, first-order non-equilibrium extrapolation for Z inlet/outlet, and periodic transverse X/Y boundaries. Viscosity in lattice units is `nu = U*D/Re` and `tau = 0.5+3*nu`.

## Visualize the *simulation* (not a standalone surface mesh)

Open `output_3d/flow_000500.vtk` in ParaView. `solid` is a voxelized obstacle field. Use **Contour** at `solid=0.5` to visualize it, and use **Slice** or **Stream Tracer** and color by `speed` or `vorticity_magnitude` to analyze the surrounding flow. Volume arrays: `velocity`, `speed`, `density`, `vorticity_magnitude`, and `solid`.

## Verified

- Compiles successfully with C++17 and passes CTest smoke test at `28 x 28 x 64` for 20 steps.
- Default `48 x 48 x 112` run completed 500 steps. Last snapshot: density 0.99483–1.00376, maximum speed 0.05036 lattice units.
- Checked the stored `solid` field: 2952 solid voxels across Z planes 19–60. The occupied cross-sectional area grows through the ogive, stays constant in the cylindrical portion, and then decreases through the boattail.

**Limitations:** This is an educational low-Reynolds-number solver. It is not validated for actual turbulent airflow, drag prediction, weapons aerodynamics, or flight performance. The voxelized surface has staircasing, the boundary conditions are approximate, and periodic crossflow boundaries may interact with the wake. No fins, moving parts, deployment geometry, or weapon-specific features are implemented.
