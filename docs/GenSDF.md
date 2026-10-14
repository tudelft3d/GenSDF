# GenSDF user guide

GenSDF is an MPI-parallel Fortran program for generating signed-distance fields
(SDFs) from triangulated surface geometry on a Cartesian grid. It supports OBJ
and binary STL input, multiple solids, uniform or stretched grids, optional
face-centred output, and CPU or NVIDIA GPU execution.

This guide is the canonical documentation for building, preparing inputs,
running examples, and interpreting outputs.

## Contents

- [Requirements and build](#requirements-and-build)
- [Running GenSDF](#running-gensdf)
- [Input file](#input-file)
- [Grid and geometry files](#grid-and-geometry-files)
- [Single and multiple solids](#single-and-multiple-solids)
- [Algorithms and sign convention](#algorithms-and-sign-convention)
- [Output files](#output-files)
- [GPU build and execution](#gpu-build-and-execution)
- [Examples](#examples)
- [Visualization](#visualization)
- [MPI, performance, and troubleshooting](#mpi-performance-and-troubleshooting)

## Requirements and build

### CPU build

The CPU build requires a Fortran compiler, `make`, and an MPI implementation
with its Fortran wrapper (for example, `mpifort`). Build from `src/`:

```sh
cd src
make
```

The executable is `src/gensdf`. To select a particular MPI Fortran compiler:

```sh
make clean
make FC=mpifort
```

Use the same MPI implementation to build and run GenSDF. The root-level
`src/parameters.in` is a sample input; when running the program, it reads
`parameters.in` from the **current working directory**, not automatically from
`src/`.

### NVIDIA GPU build

GPU support uses NVIDIA HPC SDK Fortran and OpenACC. Build with:

```sh
cd src
make clean
make GPU=1
```

If necessary, select the NVHPC MPI wrapper explicitly:

```sh
make GPU=1 FC=/path/to/nvhpc/comm_libs/openmpi/bin/mpifort
```

The Makefile's default GPU target is `cc86` (compute capability 8.6). Change
the target in the Makefile for a different GPU architecture supported by the
installed NVHPC release. The MPI wrapper, compiler, and linked OpenACC runtime
must come from a compatible NVHPC installation.

## Running GenSDF

Run the executable from a case directory containing `parameters.in` and its
`data/` directory. For example:

```sh
cd examples/armadillo
mpirun -np 4 ../../src/gensdf
```

Relative geometry and grid paths are resolved from the run directory. Grid
metadata and coordinates are read from its `data/` directory, and generated
fields are written to that `data/` directory. Ensure the case's
`data/geometry.out` and any referenced geometry/grid files are present before
launching.

For GPU execution, use the GPU-built executable and a compatible MPI runtime.
Start with one MPI rank for a small validation case, then scale out as needed;
multiple local MPI ranks can contend for the same GPU.

## Input file

`parameters.in` has nine data records, in the order below. Blank lines and
lines whose first non-whitespace character is `#` are ignored. Keep the
records in order; do not add explanatory text to a data record.

| Record | Contents | Meaning |
| --- | --- | --- |
| 1 | Geometry path | One OBJ/STL file, or a `.list`/`.lst` manifest of geometry paths |
| 2 | `scalarvalue buffer_points progressbarwidth` | Positive far-field sentinel, geometry-AABB padding in grid points, and progress display width |
| 3 | `nx ny nz` | Grid dimensions; must match `data/geometry.out` |
| 4 | `x0 y0 z0` | Physical origin |
| 5 | `.true.` or `.false.` | Enable nonuniform wall-normal grid coordinates |
| 6 | `.true.` or `.false.` | Enable Fast Sweeping distance propagation |
| 7 | `narrow_band_width` | Retained input parameter; it is not currently used as a distance cutoff |
| 8 | `vertical_axis` | `2` for solver-y wall-normal, `3` for solver-z wall-normal |
| 9 | `.true.` or `.false.` | Write optional face-centred SDF files |

Fortran logical values such as `.true.` and `.false.` are case-insensitive.
The reference files are [the root sample](https://github.com/tudelft3d/GenSDF/blob/main/src/parameters.in),
[the Armadillo case](https://github.com/tudelft3d/GenSDF/blob/main/examples/armadillo/parameters.in),
and the two
[sphere cases](https://github.com/tudelft3d/GenSDF/tree/main/examples).

`data/geometry.out` is required in every case, including uniform-grid cases.
It supplies grid dimensions and physical domain lengths. Its dimensions must
agree with record 3. `scalarvalue` initializes points without a seeded
distance, so choose a positive value larger than the greatest distance
expected in the domain. For use with dopamine IBM, choose `buffer_points` at
least `n_image_layers + 2`. See [Grid and geometry files](#grid-and-geometry-files)
for the file layout.

## Grid and geometry files

### Grid metadata

`data/geometry.out` contains the grid metadata used by the solver: the grid
dimensions and physical domain lengths. Use the layout already present in
the reference case files; the dimensions must match `nx ny nz` in
`parameters.in`.

For uniform grids, the origin comes from record 4 and the spacing is derived
from the dimensions and physical lengths. No `grid.out` file is needed.

For a stretched wall-normal grid, set record 5 to `.true.` and provide
`data/grid.out`. It must contain one five-column row per wall-normal grid
point. The solver uses columns 2 and 3 as the point-centre and face
coordinates, respectively; the other columns are read for format
compatibility:

```text
column_1  point_centre  face_coordinate  column_4  column_5
```

`data/geometry.out` has two data lines after any comments: first `nx ny nz`,
then the physical domain lengths `Lx Ly Lz`. For example:

```text
255 127 127
2.0 0.8 0.8
```

The dimensions must match record 3. The wall-normal coordinate is interpreted
according to `vertical_axis`. When `vertical_axis=2`, an optional
`data/grid_z.out` can provide stretched solver-z/spanwise coordinates using
the same five-column format and row count for the solver-z dimension;
otherwise the spanwise coordinates follow the uniform-grid setup. Use the
included case data as format references rather than assuming that a grid
file from another case has compatible dimensions or axis ordering.

### Surface geometry

Supported input is triangulated OBJ or binary STL:

- OBJ faces must be triangles. The parser accepts face references in
  `v//vn` and `v/vt/vn` form; provide vertex normals as expected by the
  parser. Consistent, outward-facing normals are recommended for distance
  initialization, even though final inside/outside sign is independent of
  triangle winding.
- Binary STL is recognized by the `.stl` or `.STL` extension.
- Geometry coordinates must use the same physical coordinate system and
  units as the grid.

The mesh must enclose a volume for meaningful inside/outside classification.
Ray-parity voting reduces sensitivity to an individual ray or triangle hit,
but it does not make an open or self-intersecting surface a valid solid.
For other CAD formats, a Python `trimesh` workflow can export a triangulated
mesh with normals, for example:

```python
import trimesh

mesh = trimesh.load("mygeometry.step")
mesh.fix_normals()
mesh.export("mygeometry.obj", include_normals=True)
```

## Single and multiple solids

Record 1 can name one surface file or a manifest with a lowercase `.list` or
`.lst` extension. A manifest contains one geometry path per nonblank,
noncomment line; paths are relative to the run directory:

```text
data/solid_a.obj
data/solid_b.stl
```

Solids receive integer IDs in manifest order, starting at 1. The supported
maximum is 15 solids. Background points have ID 0. If solids overlap, the
first manifest entry wins the inside/object-ID classification at an
overlapping point.

The object-ID field is written to `sdfp_objid.bin` on every run.

## Algorithms and sign convention

GenSDF computes unsigned distances from the triangulated surface and
propagates them over the grid with Fast Sweeping when enabled. Initial
surface-distance work is restricted to the geometry bounding box padded by
`buffer_points`; Fast Sweeping then propagates distance over the full domain.
The `narrow_band_width` record is retained for input compatibility but is not
currently applied as a truncation distance.

Inside/outside sign and solid IDs are determined on the root rank by
three-axis ray-crossing majority voting. This classification is independent
of triangle winding. The resulting signed distance is negative inside a
solid, positive outside, and zero on the surface (subject to grid and floating
point resolution).

## Output files

The program writes big-endian, double-precision binary data in the run
directory:

| File | Description |
| --- | --- |
| `data/sdfp.bin` | Point-centred signed distance; always written |
| `data/sdfp_objid.bin` | Point-centred solid ID (`0` for background); always written |
| `data/sdfu.bin` | Optional face-centred SDF in the streamwise direction |
| `data/sdfv.bin` | Optional face-centred SDF in the solver-y direction |
| `data/sdfw.bin` | Optional face-centred SDF in the solver-z direction |

The face-centred files are written only when record 9 is `.true.`. All outputs
are stored under `data/`. Point data are padded along the streamwise and
wall-normal dimensions. For interpretation and visualization, account for the
selected `vertical_axis`, case dimensions, padding, Fortran array ordering,
and big-endian byte order. Do not assume the flat binary files are
interchangeable with a native-endian or unpadded array.

## GPU build and execution

`GPU=1` enables OpenACC code compiled by NVIDIA HPC SDK. The accelerated paths
include narrow-band distance computation (including spatial-hash seeding)
and Fast Sweeping. Ray-parity sign/object-ID classification remains a
root-rank CPU operation.

Each local MPI rank selects a GPU by rank modulo the number of visible GPUs.
This assignment is convenient for simple runs, but does not account for
multiple MPI jobs sharing a device. Set GPU visibility and MPI rank placement
appropriately for the target system.

The spatial-hash path checks its per-cell candidate capacity. If a case is too
coarse or has too many triangles in a hash cell, the GPU path reports that the
candidate limit was exceeded instead of silently truncating candidates.
Improve the grid/geometry resolution match or reduce mesh density and retry;
do not treat raising `narrow_band_width` as a fix.

The GPU build and a compact runtime smoke test have been exercised with
NVHPC 23.3, OpenMPI 4.0.5, and an RTX 3060 (compute capability 8.6). CUDA
launch notifications confirmed execution of the distance and GPU sweep
kernels. This is a smoke test, not an exhaustive GPU/CPU equivalence or
performance certification; validate representative production cases on the
target hardware.

## Examples

### Armadillo workflow

The current case lives in `examples/armadillo/` (not the older singular
`example/` path). From the repository root:

```sh
cd examples/armadillo/data
bash fetchArmadillo.sh
python3 rescale_geometry.py
cd ..
mpirun -np 4 ../../src/gensdf
```

The script downloads the source mesh and rescaling uses `trimesh`; install
that Python package in the environment used to prepare the geometry if it is
not already available. The solver itself does not require Python. Inspect
`examples/armadillo/parameters.in` and `examples/armadillo/data/` before
running to confirm the case files and mesh path are set up for your checkout.

### Sphere cases

The two self-contained reference workflows are:

- [DNS sphere example](https://github.com/tudelft3d/GenSDF/tree/main/examples/case05_sphere_DNS_Re200):
  sphere case with `vertical_axis=2`.
- [LES sphere example](https://github.com/tudelft3d/GenSDF/tree/main/examples/case06_sphere_LES_Re3700):
  sphere case with `vertical_axis=3`.

Each directory provides `generate_sphere.py`, `parameters.in`, and `run.sh`.
From either case directory, run `bash run.sh 4` (replace `4` with the desired
rank count). The script generates the case geometry and launches GenSDF with
the required working directory. Python and NumPy are required to generate the
sphere mesh.

## Visualization

The repository includes helper scripts and XMF examples under `python/` and
`examples/armadillo/data/`. These are starting points for reading the raw
binary output in tools such as ParaView; they are not general metadata files
generated for every case. Check the dimensions, byte order, axis mapping,
padding, and output filename in an XMF descriptor against the case being
visualized.

`sdfp_objid.bin` can be used to distinguish background and solid regions.
`python/plot_sdf.py` provides a simple SDF inspection path; review its
expected inputs and array shape before using it with a different grid layout.

## MPI, performance, and troubleshooting

- **Input or geometry file not found:** Run from the case directory. Paths in
  `parameters.in` and manifests are run-directory-relative.
- **Dimension mismatch:** Make `nx ny nz` agree with `data/geometry.out`, and
  verify that the grid coordinate files correspond to those dimensions and
  the selected wall-normal axis.
- **Small geometry features are missed or SDF is noisy:** Ensure each solid
  spans several grid cells in its thinnest direction. Features narrower than
  about two cells may be missed by grid-based sign classification. A surface
  that is much coarser than the computational grid can also introduce
  distance artifacts; refine or remesh it as appropriate.
- **GPU hash candidate limit:** The geometry is too dense for the selected
  hash-cell population. Use a better-matched or finer grid, or a suitably
  reduced mesh, then rerun and validate.
- **Unsupported GPU/compiler target:** Confirm that `FC`, NVHPC, OpenMPI, the
  OpenACC runtime, and the Makefile's GPU architecture target are compatible.
- **Unexpected sign or solid IDs:** Check that each mesh is closed, has
  consistent physical units and placement, and that manifest order is
  intentional. Interior overlap is assigned to the first manifest entry.
- **Bad visualization:** Confirm big-endian doubles, Fortran storage order,
  padding, and whether `vertical_axis` is 2 or 3 before interpreting values.

For MPI, start with a small rank count and increase it only after checking
that the case is correct. Use no more ranks than streamwise x-cells so each
rank receives a non-empty slab. A GPU run should also account for
rank-to-device placement and GPU sharing on the host.

## Publication

Associated manuscript: [SoftwareX article](https://www.sciencedirect.com/science/article/pii/S2352711025000846).

Contributors: Akshay and Udhaya Chandiran.

## Manually publishing this guide to the GitHub Wiki

Run `./scripts/sync-wiki.sh` from a checkout of the branch whose documentation
you want to publish. The script clones the wiki using the same HTTPS or SSH
transport as the checkout's `origin`, copies this guide to the wiki's
`Home.md`, updates repository links from `main` to the checked-out branch,
then commits and pushes the change. It only updates `Home.md`; other wiki
pages are left unchanged.

The wiki must be enabled and initialized with its first page. Configure Git
credentials with wiki write access before running the script. For SSH, the
SSH key used for the source repository must also have wiki write access. For
HTTPS, use a credential helper or another secure Git credential mechanism;
do not put a personal access token in the script or its URL. Configure
`git user.name` and `git user.email` so the script can create the wiki commit.
If the checkout is in detached-HEAD state, check out the desired branch
before running.
