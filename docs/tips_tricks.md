## Tips & Tricks

With any software there are some tips and tricks that allow the user to efficiently and effectively get the results as desired. In this document, some details are provided to avoid the common pitfalls in the workflow.

1. **Supported geometry formats**: GenSDF reads triangulated OBJ files (with vertex normals) and binary STL files. The format is detected from the file extension (`.stl` / `.STL` → binary STL; anything else → OBJ). For other CAD formats, load the geometry using trimesh and export as OBJ or binary STL:
   ```python
   import trimesh
   mesh = trimesh.load("mygeometry.step")
   mesh.export("mygeometry.obj")          # OBJ with normals
   # or
   mesh.export("mygeometry.stl")          # binary STL
   ```
2. OBJ geometry must contain vertex normals. The final inside/outside sign is determined by three-axis ray-crossing parity, rather than by the winding of the nearest triangle, but consistent outward-facing normals are still recommended for well-behaved distance initialization. If normals are missing, pass the OBJ through trimesh and re-export:
   `mymesh.export("mygeometry.obj", include_normals=True)`
3. **Multiple solids**: put one geometry path per line in a `.list` manifest. Blank lines and lines beginning with `#` are ignored, and the order determines IDs (1–15). Paths are relative to the run directory. Avoid overlapping solids because the first manifest entry wins for inside cells in overlaps.
4. **Grid resolution**: each solid should span at least a few cells in its narrowest direction. GenSDF warns when an object is under-resolved; a feature narrower than two cells can be missed by any grid-based sign classifier.
5. **Stretched grids**: `data/grid.out` provides the stretched wall-normal grid when `non_uniform_grid = .true.`. For `vertical_axis = 2`, an available `data/grid_z.out` is read for a stretched solver-z/spanwise grid.
6. **Surface resolution**: if the computational grid is much finer than a coarse surface triangulation, the SDF can have artifacts. Consider using the `wrapwrap` [https://github.com/ipadjen/wrapwrap] tool to shrink-wrap the geometry with a finer triangulation, especially for cities, square blocks, and other flat or box-like geometry.
