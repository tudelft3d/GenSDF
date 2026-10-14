# GenSDF

GenSDF is an MPI-parallel Fortran program that generates signed-distance fields
from triangulated OBJ or binary STL geometry on a Cartesian grid.

The complete build, input, geometry, grid, output, GPU, visualization, and
example guide is maintained in one place:
[`docs/GenSDF.md`](docs/GenSDF.md).

## Quick start

```sh
cd src
make
cd ../examples/armadillo
mpirun -np 4 ../../src/gensdf
```

Prepare the Armadillo OBJ first as described in the guide. The executable reads
`parameters.in` and the `data/` files relative to the directory from which it
is launched.
