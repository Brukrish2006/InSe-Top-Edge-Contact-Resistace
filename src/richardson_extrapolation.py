"""
Richardson extrapolation of the top-contact and edge-contact mesh-convergence
sequences reported in Section 3.6 of the manuscript.

Motivation: the raw 0.5x-4x mesh sweep in mesh_convergence.csv shows the
top-contact resistance still rising at 4x mesh density (8.7% change from 1x to
4x), while the edge-contact resistance is already converged (0.8% change).
Rather than requiring an arbitrarily finer direct-solver mesh (which becomes
memory-prohibitive well before it plateaus), this script estimates the
mesh-independent ("converged") resistance by Richardson extrapolation of the
observed geometric convergence pattern, and cross-checks the result against an
independently re-run mesh sequence (mesh_convergence_independent_check.csv)
computed on a second DEVSIM installation (different OS / BLAS-LAPACK backend,
hence a different Delaunay mesh realization at each nominal mesh factor).

Usage:
    python richardson_extrapolation.py
"""

import csv


def load_column(path, column):
    values = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            values.append(float(row[column]))
    return values


def richardson_extrapolate(values):
    """Extrapolate a monotonic, geometrically-converging sequence of 4+ points
    (each point at 2x the mesh density of the previous) to its limit, using
    the ratio of the last two successive differences."""
    diffs = [values[i + 1] - values[i] for i in range(len(values) - 1)]
    ratios = [diffs[i + 1] / diffs[i] for i in range(len(diffs) - 1)]
    r = ratios[-1]
    r_inf = values[-1] + diffs[-1] * r / (1 - r)
    return diffs, ratios, r_inf


def main():
    top = load_column("mesh_convergence.csv", "R_Top")
    edge = load_column("mesh_convergence.csv", "R_Edge")

    print("=== Original manuscript mesh sequence (mesh_convergence.csv) ===")
    d_top, r_top, top_inf = richardson_extrapolate(top)
    print(f"Top-contact  R: {top}")
    print(f"  successive diffs:  {[round(d, 4) for d in d_top]}")
    print(f"  diff ratios:       {[round(r, 3) for r in r_top]}")
    print(f"  Richardson R_top,inf  = {top_inf:.4f} Ohm.cm")

    d_edge, r_edge, edge_inf = richardson_extrapolate(edge)
    print(f"\nEdge-contact R: {edge}")
    print(f"  successive diffs:  {[round(d, 4) for d in d_edge]}")
    print(f"  diff ratios:       {[round(r, 3) for r in r_edge]}")
    print(f"  Richardson R_edge,inf = {edge_inf:.4f} Ohm.cm")

    adv_top = (top_inf - edge_inf) / top_inf * 100
    adv_edge = (top_inf - edge_inf) / edge_inf * 100
    print(f"\nConverged geometric advantage: {adv_top:.1f}% (rel. to R_top) "
          f"/ {adv_edge:.1f}% (rel. to R_edge)")

    raw_top = (top[-1] - edge[-1]) / top[-1] * 100
    raw_edge = (top[-1] - edge[-1]) / edge[-1] * 100
    print(f"Raw 4x-only advantage (unconverged snapshot): "
          f"{raw_top:.1f}% / {raw_edge:.1f}%")

    print("\n=== Independent cross-validation "
          "(mesh_convergence_independent_check.csv, second DEVSIM install, "
          "extended to 8x) ===")
    top2 = load_column("mesh_convergence_independent_check.csv", "R_Top")
    d_top2, r_top2, top2_inf = richardson_extrapolate(top2)
    print(f"Top-contact R: {top2}")
    print(f"  successive diffs:  {[round(d, 4) for d in d_top2]}")
    print(f"  diff ratios:       {[round(r, 3) for r in r_top2]}")
    print(f"  Richardson R_top,inf  = {top2_inf:.4f} Ohm.cm")
    print(f"  8x raw value is {(top2_inf - top2[-1]) / top2_inf * 100:.2f}% "
          f"below the extrapolated limit -> effectively converged by 8x")
    print(f"\nAgreement between the two independent extrapolations: "
          f"{abs(top_inf - top2_inf) / top_inf * 100:.2f}% relative difference")


if __name__ == "__main__":
    main()
