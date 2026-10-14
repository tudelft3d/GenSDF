# Input parameters

GenSDF reads `parameters.in` from the current working directory. Blank lines and
lines whose first non-space character is `#` are ignored; all other lines are
the data records below, in order.

```text
# Geometry: one OBJ/STL file or a manifest of multiple files
'data/armadillo_withnormals.obj'
# scalarvalue  buffer_points  progressbarwidth
1d10  8  40
# nx ny nz (must match data/geometry.out)
512 128 128
# Grid origin x y z
0.0 0.0 0.0
# Read stretched wall-normal grid from data/grid.out
.false.
# Use spatial-hash distance seeds and Fast Sweeping
.true.
# narrow_band_width (retained in the input format; currently not a distance cutoff)
8
# Solver wall-normal axis: 2=y, 3=z
3
# Also write face-staggered sdfu/sdfv/sdfw
.false.
```

## Records

1. **Geometry path.** A single triangulated OBJ or binary STL is one solid with
   ID 1. A filename ending in `.list` or `.lst` is a manifest: each non-blank,
   non-comment line names one OBJ or STL file. Paths are interpreted relative to
   the run directory, and manifest order assigns solid IDs starting at 1. OBJ
   geometry must contain vertex normals (`vn`); both `v//vn` and `v/vt/vn`
   triangular face tokens are supported. Binary STL is detected by `.stl` or
   `.STL`.

   The multi-solid path supports up to 15 solids to match the dopamine solver's
   IBM boundary-condition arrays. Solids should not overlap; where they do, the
   first manifest entry claims interior grid points in the overlap.

2. **`scalarvalue buffer_points progressbarwidth`.**
   - `scalarvalue` is a positive, large sentinel for points without a seeded
     distance. Choose a value larger than distances expected in the domain.
   - `buffer_points` expands the geometry bounding box in grid cells for local
     distance initialization and sign classification. For dopamine IBM, retain
     at least `n_image_layers + 2` cells.
   - `progressbarwidth` controls the terminal progress bar (typically 20–40).

3. **`nx ny nz`.** Grid dimensions in solver coordinates. These must agree with
   the dimensions in `data/geometry.out`; the latter also supplies the physical
   domain lengths. GenSDF permutes the dimensions internally when
   `vertical_axis = 2`.

4. **Grid origin.** The `x y z` origin, consistent with the solver's
   `geometry.out` and grid coordinate files.

5. **`non_uniform_grid`.** `.true.` reads the stretched wall-normal coordinates
   from `data/grid.out` (columns 2 and 3 are cell-centre and face coordinates).
   With `vertical_axis = 3`, this is the solver's z grid; with `vertical_axis =
   2`, it is the solver's y grid. `.false.` uses uniform spacing in that
   direction. When `vertical_axis = 2`, an available `data/grid_z.out` is also
   read automatically for a stretched solver-z/spanwise grid.

6. **`use_fast_sweep`.** `.false.` uses the original brute-force point-to-
   triangle distance initialization. `.true.` uses spatial-hash candidates for
   local distance seeds, ray-crossing parity to determine signs, then the
   Godunov Fast Sweeping Method to propagate distances across the domain. The
   hash-seeded initialization is not clipped by `narrow_band_width`.

7. **`narrow_band_width`.** Kept in the input format for compatibility; the
   current spatial-hash initializer does not use it as a distance cutoff.

8. **`vertical_axis`.** The solver wall-normal axis: `3` means z (default), and
   `2` means y. For axis 2, GenSDF permutes coordinates internally so its
   internal z axis is the wall-normal direction. When this axis is selected,
   stretched solver-z coordinates can additionally be supplied in
   `data/grid_z.out`.

9. **`compute_face_sdf`.** `.false.` writes only the cell-centre SDF.
   `.true.` additionally writes the u-, v-, and w-face SDFs.

The sign is classified independently of triangle winding by majority voting
over ray crossings along all three coordinate axes. This improves robustness
for multi-object scenes and meshes with a missing face, but geometry should
still represent a meaningful solid.

## Generated files

Outputs are big-endian, double-precision binary arrays in `data/`, padded in x
and the solver wall-normal direction for the IBM reader:

- `sdfp.bin`: cell-centre signed-distance field.
- `sdfp_objid.bin`: cell-centre solid ID (zero outside solids).
- `sdfu.bin`, `sdfv.bin`, `sdfw.bin`: face-staggered signed-distance fields,
  written only when `compute_face_sdf = .true.`.
