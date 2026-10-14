# GenSDF examples

Build the MPI executable once:

```sh
cd src
make
```

Each example must be run from its own directory because GenSDF reads
`parameters.in` and the `data/` files relative to the current working
directory.

## Armadillo

1. Fetch the Stanford Armadillo OBJ:

   ```sh
   cd example/data
   chmod +x fetchArmadillo.sh
   ./fetchArmadillo.sh
   ```

2. From `example/`, scale the mesh and write the required OBJ vertex normals:

   ```sh
   cd ..
   python3 data/rescale_geometry.py
   ```

3. Run GenSDF using the supplied `example/parameters.in`:

   ```sh
   mpirun -np 4 ../src/gensdf
   ```

The output fields are written under `example/data/`. The `parameters.in`
geometry path is already set to `data/armadillo_withnormals.obj`; there is no
need to copy the executable or edit the input format.

## Sphere solver cases

[`../examples/case05_sphere_DNS_Re200`](../examples/case05_sphere_DNS_Re200)
and [`../examples/case06_sphere_LES_Re3700`](../examples/case06_sphere_LES_Re3700)
provide binary-STL generators, matching `geometry.out` grid references,
parameters, and run scripts. From either case directory, run:

```sh
bash run.sh 4
```

The scripts generate `data/sphere.stl`, run `../../src/gensdf`, and write
cell-centre and face-staggered fields. NumPy and an MPI runtime are required.
The cell-centre SDF and its object-ID field can be copied to a dopamine solver
run directory as required by that solver's IBM configuration.

## MPI notes

GenSDF is MPI-only; running with a single rank is supported. Keep the rank
count commensurate with the x-grid size. Fast Sweeping exchanges data between
neighboring x-slabs, so use a non-empty slab per rank (no more ranks than
x-cells).

Always ensure the MPI library used to build GenSDF is available in the runtime
environment. For all input records and supported grid files, see
[`input_parameters.md`](input_parameters.md).
