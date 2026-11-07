"""Run the measured-pinning studies (barrier sweep + real metals) in parallel.

Usage:  python run_parallel.py [n_workers]
Writes data_v3/barrier_sweep.csv and data_v3/metals_pinning.csv.
Each task runs in a fresh process, so DEVSIM state never leaks between tasks.
"""
import multiprocessing as mp
import os
import sys
import time

# one BLAS thread per worker, otherwise N workers x M threads oversubscribe the cores
for _v in ("MKL_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def tasks():
    import studies as s
    NM = s.NM
    T = []
    Ev = s.lit_Ev()
    geoms = [("top", 50 * NM), ("edge", 50 * NM)]
    for phiB in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        for kind, Lc in geoms:
            Phi = Ev - phiB
            T.append(("barrier", "point", kind, Phi, Lc, dict(s.LIT), dict(study="ungated", phiB_set=phiB), None))
            T.append(("barrier", "transfer", kind, Phi, Lc, dict(s.LIT), dict(study="gated", phiB_set=phiB),
                      [0.0, -1.0, -2.0, -3.0, -4.0]))
            T.append(("barrier", "transfer", kind, Phi, Lc, dict(s.LIT, tunneling="off"),
                      dict(study="gated", phiB_set=phiB, tunneling="off"), [-3.0]))
    for case, c in s.pin_cases().items():
        for metal, Phi in s.METALS.items():
            for kind, Lc in geoms:
                P = dict(s.LIT, S_top=c["S"], S_edge=c["S"], Phi_CNL=c["CNL"])
                tag = dict(case=case, metal=metal, S=c["S"], CNL=c["CNL"])
                T.append(("metals", "point", kind, Phi, Lc, P, dict(tag, study="ungated"), None))
                T.append(("metals", "transfer", kind, Phi, Lc, P, dict(tag, study="gated"), [-3.0]))
    return T


def run(task):
    import studies as s
    out, how, kind, Phi, Lc, P, tag, VGs = task
    if how == "point":
        rows = [s.point(kind, Phi, Lc=Lc, P=P, extra=tag)[0]]
    else:
        rows = s.transfer(kind, Phi, 0.05, P=P, Lc=Lc, VGs=VGs, tag=tag)
    return out, rows


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, (os.cpu_count() or 2) - 2)
    import studies as s
    T = tasks()
    print(f"{len(T)} tasks on {n} workers", flush=True)
    res = {"barrier": [], "metals": []}
    t0 = time.time()
    with mp.get_context("spawn").Pool(n, maxtasksperchild=1) as pool:
        for i, (out, rows) in enumerate(pool.imap_unordered(run, T), 1):
            res[out] += rows
            print(f"[{i}/{len(T)}] {out} done, {time.time() - t0:.0f} s", flush=True)
    s.write("barrier_sweep", res["barrier"])
    s.write("metals_pinning", res["metals"])
    print("finished in", round(time.time() - t0), "s", flush=True)
