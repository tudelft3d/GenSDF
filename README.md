# GenSDF - Signed-Distance-Field Generator
MPI + Fortran software to compute the signed-distance-field for surface geometry (OBJ or binary STL) on a Cartesian grid.

Associated Manuscript: [Software-X - Click Here](https://www.sciencedirect.com/science/article/pii/S2352711025000846)

The code assumes a regular grid along the streamwise and spanwise directions, while a non-uniform grid can be used in the wall-normal direction. The wall-normal direction can be set to $y$ or $z$ via the `vertical_axis` parameter.

## How to compile the code

```
cd src
make
```

A debug build (bounds checking, no optimisation) can be triggered with:
```
cd src
make FFLAGS="-O0 -g -fbounds-check -ffree-line-length-none -cpp"
```

### GPU build (OpenACC / NVHPC)
```
cd src
make GPU=1 FC=/path/to/nvhpc/comm_libs/openmpi4/openmpi-4.0.5/bin/mpif90
```
Requires the NVHPC SDK (`nvfortran` on `PATH`) and its MPI wrapper. Targets `cc86` by default (`-acc -gpu=cc86,rdc`); edit `GPU=1`'s `FFLAGS`/`LDFLAGS` in `src/Makefile` for a different compute capability. `GPU=0` (default) is unaffected. See `docs/gpu_build.md` for what the port covers and how it was validated.

## How to change parameters
Edit `parameters.in` following the format described in `docs/input_parameters.md`. Comment lines begin with `#`.

## How to run the code
Edit `parameters.in`, then:
```
mpirun -np 8 ./gensdf
```

Example output (8 MPI ranks, cell-centre SDF only). The elapsed time shown on the progress bar is indicative; computational load may differ across ranks.

```
  ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓███████▓▒░▒▓███████▓▒░░▒▓████████▓▒░ 
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░   
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
  ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓█▓▒░        
  *** Starting with            8 MPI ranks ***
 *** Input file successfully read ***
 *** Successfully read the grid ***
 *** Successfully set up grid spacing ***
 Successfully read OBJ: data/lucy.obj
 Vertices:        49987   Normals:        49987   Faces:        99970
 Geometry bbox min:  0.325  0.400  0.200
 Geometry bbox max:  0.675  0.600  0.800
 AABB x:  156  0.305  |  355  0.693
 AABB y:  195  0.380  |  317  0.618
 AABB z:   92  0.179  |  420  0.819
 -- Pre-processing done in   0.27 s
 -- Estimated minimum memory:  0.14  GiB(s)
 *** Calculating SDF | cell-centres ***
||||||||||||||||||||||100.00%    99970/   99970 Elapsed:     3.57s Remaining:     0.00s
 *** Writing | cell-centres ***
 -- Write done in   0.23 s | cell-centres
 *** SDF complete in   23.87 s ***
```
<center><img src="armadillo.png" height=400></center>

<center> 
Figure: Comparison of the calculated SDF around the Stanford Armadillo geometry scaled to a smaller size
</center>

## Output files
All SDF fields are written to the `data/` directory as big-endian, double-precision binary files, ghost-padded in the streamwise ($x$) and wall-normal directions to match the IBM reader convention of the target solver:
- `sdfp.bin` — cell-centre SDF (always written)
- `sdfu.bin`, `sdfv.bin`, `sdfw.bin` — face-staggered SDFs (written only when `compute_face_sdf = .true.`)

## How to visualise the results
1. Use `python/read_xmf_constant_grid.xmf` (uniform grid spacing) or `python/read_binary_output_in_paraview.xmf` (non-uniform grid) to load the binary output files directly in Paraview via the XDMF reader. Use the default XDMF reader — other readers may crash the session.
2. For non-uniform grids, generate a matching XMF descriptor with `python/write_xmf_grid.py`.
3. Use `python/readmask.py` to read the binary SDF output in Python and plot slices.
4. To visualise the geometry surface, plot the contour at a value of 0.

# Contributors
Akshay & Udhaya [https://github.com/udhaya-chandiran]
