#!/usr/bin/env bash
# =============================================================
#  Run GenSDF for solver Example 06 — sphere LES Re_D ≈ 3700
#
#  Sphere D = 0.1 m, centre at (0.4, 0.4, 0.4), domain 1.6×0.8×0.8 m
#
#  Usage:
#    cd examples/case06_sphere_LES_Re3700
#    bash run.sh [N_PROCS]     # N_PROCS default = 4
# =============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

NPROCS="${1:-4}"
GENSDF="../../src/gensdf"

if [[ ! -x "${GENSDF}" ]]; then
    echo "ERROR: ${GENSDF} not found or not executable."
    echo "       Compile first:  cd ../../src && make"
    exit 1
fi
if [[ ! -f data/geometry.out ]]; then
    echo "ERROR: data/geometry.out is required and must match parameters.in."
    exit 1
fi

echo "Step 1: generating sphere geometry (data/sphere.stl)..."
python3 generate_sphere.py

echo ""
echo "Step 2: running GenSDF with ${NPROCS} MPI ranks..."
mpirun -np "${NPROCS}" "${GENSDF}"

echo ""
echo "Done.  SDF files written to data/:"
ls -lh data/sdf*.bin 2>/dev/null || true
echo ""
echo "To use with solver Example 06:"
echo "  cp data/sdfp.bin  <run_dir>/SDF_in"
echo "  Ensure ibm_sdf_input = 1 in solver input_parameters."
