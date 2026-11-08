"""Run the revision studies in parallel (one DEVSIM device per fresh process).

Usage:  python run_revision.py <n_workers> <study> [<study> ...]
        studies: uq_gated thermal_lit mesh_convergence thickness_lit inse_electron   (or 'all')
Smoke test first:  UQ_N=2 python run_revision.py 4 uq_gated
Writes data_v3/<study>.csv (partial results are re-written every 20 tasks).
"""
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def work(spec):
    import studies_rev as r
    return r.run_dev(spec)


if __name__ == "__main__":
    import studies as s
    import studies_rev as r
    n = int(sys.argv[1])
    names = sys.argv[2:]
    if names == ["all"]:
        names = list(r.STUDIES)
    for name in names:
        T = r.STUDIES[name]()
        print(f"== {name}: {len(T)} tasks on {n} workers", flush=True)
        rows, t0 = [], time.time()
        with mp.get_context("spawn").Pool(n, maxtasksperchild=1) as pool:
            for i, (_, rr) in enumerate(pool.imap_unordered(work, T), 1):
                rows += rr
                print(f"[{name} {i}/{len(T)}] {time.time() - t0:.0f} s", flush=True)
                if i % 20 == 0:
                    s.write(name, rows)
        s.write(name, rows)
        print(f"== {name} finished in {time.time() - t0:.0f} s", flush=True)
