# Serial D3Q19 LBM: 3D flow around a generic ogive body

This is a **working single-threaded C++17 3D CFD simulator**, not a 3D mesh generator. The solid body is computed **inside the fluid solver** as a voxelized bounce-back obstacle. There is no STL input or external surface-meshing step.

The central body is a **generic streamlined axisymmetric test body** with tangent-ogive nose, cylindrical section, and linear boattail ending in a flat rear cap.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build --output-on-failure
./build/lbm_3d --steps 500 --every 250 --out output_3d
```

The solver requires only a C++17 toolchain and CMake. Serial execution only: no CUDA, MPI, OpenMP, GPU, or Kokkos. The optional reader regression test requires Python 3 with the VTK Python modules. CMake enables it when those modules are importable by the Python interpreter it finds; set `-DPython3_EXECUTABLE=/path/to/python` when configuring to select an interpreter. CMake reports when the reader test is disabled. Run `ctest --test-dir build -R vti_reader --output-on-failure` to run it separately.

## Benchmark the serial optimizations

```sh
python3 benchmarks/compare.py
```

The script builds four single-threaded binaries with the same compiler and `-O0`. The original comes from Git revision `aaad885`. Applying `benchmarks/pre-simd.patch` to that revision reconstructs the earlier scalar optimizations. The working tree builds twice to measure the new storage layout with scalar collision (`-DLBM_DISABLE_SIMD`) and with runtime SIMD dispatch. The script requires Python 3, Git with the original revision available locally, and a GCC- or Clang-compatible compiler.

Each row stores the X cells for a given D3Q19 direction consecutively. On AVX-capable x86 systems with GCC or Clang, collision loads four neighboring fluid cells into each vector. Both BGK and Smagorinsky collision use this path. Short runs of fluid cells use scalar collision, including any remainder at the end of a row. Other supported architectures and compilers use the scalar implementation. Streaming follows the contiguous X spans and sums reflected momentum in the original cell and direction order.

The default workload is 48 × 48 × 112 for 24 steps. After one warm-up per binary, four timed rounds place each binary in each execution position once. The script prints each wall time, the medians, and speedup ratios. Timing includes initialization, stepping, diagnostics, and writing two snapshots. Compilation and byte comparisons of both snapshots and the force history happen outside the timed interval.

Use `--grid NX NY NZ`, `--steps N`, `--smagorinsky CS`, `--repeats N`, and `--cxx COMPILER` to change the run. On Linux, prefix the command with `taskset -c CPU` to keep all four binaries on the same logical CPU. The [SIMD benchmark report](benchmarks/simd-report.txt) includes the commands and complete output from this machine.

## Model and flow axes

The **longitudinal body axis is Z** and the inlet velocity points in **+Z**. `nx,ny` are the crossflow lattice extents, `nz` is streamwise. For the central body alone, at each lattice cell `(x,y,z)`:

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

Four rectangular plates extend from the rear taper along the ±X and ±Y directions. Each span tilts 30° upstream from radial (60° to the nose direction, −Z), with the wide face exposed to the flow. Using the maximum body radius R, the nominal exposed length is 1.5 × 1.25R, the lateral width is 2.0R, and the thickness is 0.2R with a minimum of 1.5 lattice units for voxelization. The roots overlap the taper half a radius ahead of the rear face. The solver includes the plates in the same solid mask as the body.

For better boundary resolution, increase all three lattice extents while preserving the relevant ratios.

## Parameters

| Argument | Default | Description |
|---|---:|---|
| `--nx` | 288 | Lattice cells across X |
| `--ny` | 288 | Lattice cells across Y |
| `--nz` | 672 | Lattice cells along Z (flow direction) |
| `--steps` | 500 | Number of time steps |
| `--every` | 250 | Snapshot interval |
| `--u` | 0.045 | Uniform +Z inlet velocity in lattice units |
| `--re` | 35 | Nominal Reynolds number based on maximum diameter |
| `--smagorinsky` | 0 | Experimental, unvalidated sub-grid coefficient |
| `--out` | output_3d | VTK volume output directory |

The model uses D3Q19 BGK collision and pull streaming, no-slip halfway bounce-back at solid links, first-order non-equilibrium extrapolation for Z inlet/outlet, and periodic transverse X/Y boundaries. Viscosity in lattice units is `nu = U*D/Re` and `tau = 0.5+3*nu`.

## Visualize the simulation in ParaView

The solver writes `flow_000000.vti` at initialization and subsequent numbered `.vti` snapshots at the requested interval and final step. Open the `.vti` file series in ParaView and click **Apply**; the animation controls select snapshots. The six-digit numbers are simulation step indices, not physical times. Use a fresh output directory for each run to avoid mixing old snapshots with new ones.

With 500 steps and an interval of 250, the three snapshots show initialization, step 250 and step 500. The solid geometry is fixed throughout. The terminal prints the instantaneous `Cd` at each snapshot interval and its mean over the last quarter of the run at completion. The accompanying `forces.csv` records step, force components and drag coefficient every `--force-every` steps and at the last step. Legacy `.vtk` snapshots in older output directories belong to earlier runs.

The default grid is 288 × 288 × 672. Each snapshot occupies about 1.3 GiB. The solver advances one Z plane at a time and keeps two extra planes for streaming, with about 9 GiB of array storage at snapshot time. At the same lattice inlet speed, doubling each grid dimension requires twice as many steps to cover the same distance in body lengths.

The XML ImageData volume has unit spacing and origin `(0,0,0)` in lattice units, with dimensions `nx × ny × nz`. All arrays are **point data**, sampled at integer lattice coordinates, with X varying fastest. `velocity` is a three-component vector; `speed`, `density`, `vorticity_magnitude`, and `solid` are scalars. `solid` is 1 inside the body and 0 in the fluid. Flow is in +Z, so the wake is downstream toward larger Z.

1. Select the volume, add **Contour**, choose the point array `solid`, and set the isovalue to `0.5`. Apply and show the resulting surface. It interpolates between lattice samples and approximates the voxelized boundary.
2. Select the original volume and add **Slice** through the center of the body, with a normal along X or Y. Color it by `speed`, `density`, or `vorticity_magnitude`; rescale the color range if necessary.
3. Select the original volume and add **Stream Tracer**. Choose `velocity` as the vectors and place a seed line or plane upstream of the nose, across X/Y at a small Z. Color the resulting lines by `speed`. Adjust the seed and integration settings to inspect the wake.
4. For a vortex view, color a slice by `vorticity_magnitude` or apply a contour to that array and choose an isovalue based on its displayed range. Use the animation controls to see how the wake evolves.

Solid points have zero velocity, speed, and vorticity, and density is set to 1 as a placeholder; they are not fluid measurements. The curl uses periodic X/Y differences and one-sided differences at the Z ends. Its finite differences use zero velocity at solid neighbors, so near-wall values reflect the voxelized boundary. Streamlines and contours near the surface are interpolated from these point samples.

## Verification

`ctest` runs a small solver smoke test. When VTK Python modules are available it also runs `lbm_3d_vti_reader`, which generates fresh initial and final snapshots and uses VTK's XML reader to check dimensions, coordinates, arrays, values, and a nonempty obstacle contour. The original legacy files were readable by VTK in the checked sample; no specific corruption was established. XML ImageData makes the array types, association, byte order, and binary block lengths explicit.

**Limitations:** This is an educational low-Reynolds-number solver. It is not validated for actual turbulent airflow, drag prediction, or flight performance. The voxelized surface has staircasing, the boundary conditions are approximate, and periodic crossflow boundaries may interact with the wake. The four rear plates are fixed throughout the simulation.
