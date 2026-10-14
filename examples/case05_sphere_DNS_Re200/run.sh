#!/usr/bin/env bash
# =============================================================
#  Run GenSDF for solver Example 05 — sphere DNS Re_D = 200
#
#  Steps performed:
#    1. Generate data/sphere.stl  (UV sphere, D=0.1 m at (0.5,0.4,0.4))
#    2. Run gensdf to compute SDF on the solver grid
#    3. Print instructions for copying output to solver run directory
#
#  Prerequisites:
#    - NumPy installed (for generate_sphere.py):   pip install numpy
#    - GenSDF compiled:  cd ../../src && make
#    - MPI runtime available (mpirun / mpiexec)
#
#  Usage:
#    cd examples/case05_sphere_DNS_Re200
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

# --- Step 1: generate sphere geometry ---
echo "Step 1: generating sphere geometry (data/sphere.stl)..."
python3 generate_sphere.py

# --- Step 2: run GenSDF ---
echo ""
echo "Step 2: running GenSDF with ${NPROCS} MPI ranks..."
mpirun -np "${NPROCS}" "${GENSDF}"

# --- Step 3: usage hint ---
echo ""
echo "Done.  SDF files written to data/:"
ls -lh data/sdf*.bin 2>/dev/null || true
echo ""
echo "To use with solver Example 05:"
echo "  cp data/sdfp.bin  <run_dir>/SDF_in"
echo "  cp data/sdfu.bin  <run_dir>/sdfu.bin   # optional: future staggered IBM"
echo "  Ensure ibm_sdf_input = 1 in solver input_parameters."
