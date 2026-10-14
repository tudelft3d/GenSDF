#!/usr/bin/env python3
"""
generate_sphere.py — create a UV sphere as binary STL.

Case 06: Flow around a sphere — LES, Re_D ≈ 3700
  Sphere diameter D = 0.1 m, centre at (0.4, 0.4, 0.4) m
  Domain: Lx=1.6 × Ly=0.8 × Lz=0.8 m  (sphere at 4D from inlet)

Output: data/sphere.stl

Usage:
    python3 generate_sphere.py

The generated STL is consumed by GenSDF (gensdf) to compute the
signed-distance field on the solver's Cartesian grid.  The normals point
outward (away from the sphere centre); GenSDF determines the final sign by
three-axis ray-crossing parity.
"""

import os
import struct
import numpy as np

# ---------------------------------------------------------------------------
# Sphere parameters — edit to match your case
# ---------------------------------------------------------------------------
CX, CY, CZ = 0.4, 0.4, 0.4   # sphere centre (m)  — 4D from inlet
R          = 0.05              # sphere radius (m)  — D = 0.1 m
N_PHI      = 32                # polar angle divisions
N_THETA    = 64                # azimuthal angle divisions
OUT_FILE   = 'data/sphere.stl'
# ---------------------------------------------------------------------------


def sphere_pt(phi: float, theta: float) -> np.ndarray:
    """
    Return the Cartesian point on the sphere surface.
    y = vertical axis (phi=0 → north pole at y = CY + R).
    """
    return np.array([
        CX + R * np.sin(phi) * np.cos(theta),
        CY + R * np.cos(phi),
        CZ + R * np.sin(phi) * np.sin(theta),
    ], dtype=np.float64)


def outward_normal(v0: np.ndarray, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Unit normal pointing away from the sphere centre."""
    centre   = np.array([CX, CY, CZ], dtype=np.float64)
    centroid = (v0 + v1 + v2) / 3.0
    radial   = centroid - centre
    n        = np.cross(v1 - v0, v2 - v0)
    if np.dot(n, radial) < 0.0:
        n = -n
    return n / max(np.linalg.norm(n), 1.0e-30)


def build_triangles():
    tris = []
    phi_vals = np.linspace(0.0, np.pi, N_PHI + 1)

    for i in range(N_PHI):
        phi1, phi2 = phi_vals[i], phi_vals[i + 1]
        for j in range(N_THETA):
            t1 = 2.0 * np.pi * j       / N_THETA
            t2 = 2.0 * np.pi * (j + 1) / N_THETA

            A = sphere_pt(phi1, t1)
            B = sphere_pt(phi2, t1)
            C = sphere_pt(phi2, t2)
            D = sphere_pt(phi1, t2)

            if i == 0:
                n = outward_normal(A, B, C)
                tris.append((n, A, B, C))
            elif i == N_PHI - 1:
                n = outward_normal(A, B, D)
                tris.append((n, A, B, D))
            else:
                tris.append((outward_normal(A, B, C), A, B, C))
                tris.append((outward_normal(A, C, D), A, C, D))

    return tris


def write_stl_binary(tris: list, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    hdr = b'Binary STL sphere - generate_sphere.py (dopamine IBM example)'
    hdr = hdr[:80].ljust(80, b'\x00')
    with open(path, 'wb') as f:
        f.write(hdr)
        f.write(struct.pack('<I', len(tris)))
        for n, v0, v1, v2 in tris:
            data = np.array([*n, *v0, *v1, *v2], dtype=np.float32)
            f.write(data.tobytes())
            f.write(struct.pack('<H', 0))


if __name__ == '__main__':
    print(f"Generating UV sphere STL:")
    print(f"  Centre  = ({CX}, {CY}, {CZ}) m")
    print(f"  Radius  = {R} m   Diameter = {2*R} m")
    print(f"  Mesh    = {N_PHI} polar × {N_THETA} azimuthal divisions")

    tris = build_triangles()
    write_stl_binary(tris, OUT_FILE)

    expected = 2 * N_THETA + (N_PHI - 2) * N_THETA * 2
    print(f"  Triangles written: {len(tris)}  (expected {expected})")
    print(f"  Output: {OUT_FILE}")
    print()
    print("Next step: bash run.sh [N]   (run from this directory)")
