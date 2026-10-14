## GPU build (OpenACC) -- This document was created with the help of AI

`GPU=1` builds an OpenACC path (`-DGPU_SDF`) alongside the default CPU/MPI build. Same
toolchain as NVHPC 23.3 (`nvfortran`, `-acc -gpu=cc86`). The CPU path (`GPU=0`, default) is
untouched — all GPU code is guarded by `#ifdef GPU_SDF` or lives in plain `!$acc` comment
lines that gfortran ignores.

### Validation

Built and run on this machine (NVHPC 23.3, RTX 3060, cc86) against a synthetic UV-sphere
test case:

- `vertical_axis=2` vs `vertical_axis=3` on the same geometry: interior cell counts and
  minimum SDF value agree closely (discretization-level differences only), confirming the
  y-z winding fix in `utils_io.f90` is correct.
- CPU vs GPU output: identical away from a handful of cells at the test sphere's pole
  (a UV-sphere singularity where many near-degenerate triangles converge to a point);
  those disagree at the FP-tie level, not a systematic offset.
- CPU `nprocs=1` vs `nprocs=2`: byte-identical for the brute-force distance path;
  Gauss-Seidel FSM sweep order is decomposition-dependent at near-tie cells, a pre-existing
  property of the algorithm rather than something the GPU port introduced.

### Known limitation

The spatial hash's cell size is tuned for grid resolution, not triangle density. A geometry
whose bounding box is small relative to the grid can concentrate many triangles into one hash
cell; `max_cands_gpu=32768` bounds the per-thread candidate buffer, and `build_hash` will
`error stop` if a denser mesh or tighter bounding box would exceed it, rather than silently
truncating candidates.
