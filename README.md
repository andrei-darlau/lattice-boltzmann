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

The solver requires only a C++17 toolchain and CMake. Serial execution only: no CUDA, MPI, OpenMP, GPU, or Kokkos. The optional reader regression test requires Python 3 with the VTK Python modules. CMake enables it when those modules are importable by the Python interpreter it finds; set `-DPython3_EXECUTABLE=/path/to/python` when configuring to select an interpreter. CMake reports when the reader test is disabled. Run `ctest --test-dir build -R vti_reader --output-on-failure` to run it separately.

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

## Visualize the simulation in ParaView

The solver writes `flow_000000.vti` at initialization and subsequent numbered `.vti` snapshots at the requested interval and final step. Open the `.vti` file series in ParaView and click **Apply**; the animation controls select snapshots. The six-digit numbers are simulation step indices, not physical times. Use a fresh output directory for each run to avoid mixing old snapshots with new ones.

The XML ImageData volume has unit spacing and origin `(0,0,0)` in lattice units, with dimensions `nx × ny × nz`. All arrays are **point data**, sampled at integer lattice coordinates, with X varying fastest. `velocity` is a three-component vector; `speed`, `density`, `vorticity_magnitude`, and `solid` are scalars. `solid` is 1 inside the body and 0 in the fluid. Flow is in +Z, so the wake is downstream toward larger Z.

1. Select the volume, add **Contour**, choose the point array `solid`, and set the isovalue to `0.5`. Apply and show the resulting surface. It interpolates between lattice samples and approximates the voxelized boundary.
2. Select the original volume and add **Slice** through the center of the body, with a normal along X or Y. Color it by `speed`, `density`, or `vorticity_magnitude`; rescale the color range if necessary.
3. Select the original volume and add **Stream Tracer**. Choose `velocity` as the vectors and place a seed line or plane upstream of the nose, across X/Y at a small Z. Color the resulting lines by `speed`. Adjust the seed and integration settings to inspect the wake.
4. For a vortex view, color a slice by `vorticity_magnitude` or apply a contour to that array and choose an isovalue based on its displayed range. Use the animation controls to see how the wake evolves.

Solid points have zero velocity, speed, and vorticity, and density is set to 1 as a placeholder; they are not fluid measurements. The curl uses periodic X/Y differences and one-sided differences at the Z ends. Its finite differences use zero velocity at solid neighbors, so near-wall values reflect the voxelized boundary. Streamlines and contours near the surface are interpolated from these point samples.

## Verification

`ctest` runs a small solver smoke test. When VTK Python modules are available it also runs `lbm_3d_vti_reader`, which generates fresh initial and final snapshots and uses VTK's XML reader to check dimensions, coordinates, arrays, values, and a nonempty obstacle contour. The original legacy files were readable by VTK in the checked sample; no specific corruption was established. XML ImageData makes the array types, association, byte order, and binary block lengths explicit.

**Limitations:** This is an educational low-Reynolds-number solver. It is not validated for actual turbulent airflow, drag prediction, weapons aerodynamics, or flight performance. The voxelized surface has staircasing, the boundary conditions are approximate, and periodic crossflow boundaries may interact with the wake. No fins, moving parts, deployment geometry, or weapon-specific features are implemented.
