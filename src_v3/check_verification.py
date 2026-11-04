"""Fail (exit 1) if any v3 verification metric exceeds its tolerance."""
import csv, os, sys
TOL = {"V1": 1e-6, "V2": 2e-3, "V4": 1.5e-2, "V5": 2e-2, "V6": 2e-2, "V7": 1e-3, "V8": 1e-2, "V9": 1e-5}
SKIP = ("J/J_MG at 5 V", "J/J_MG at 20 V", "J/J_MG at 80 V", "depletion", "share", "|I(V=0)|", "max|phi_p|")
rows = list(csv.DictReader(open(os.path.join(os.path.dirname(__file__), "..", "data_v3", "verification_summary.csv"))))
bad = []
for r in rows:
    t = r["test"]; q = r["quantity"]
    if any(s in q for s in SKIP) or t not in TOL:
        continue
    if abs(float(r["rel_error"])) > TOL[t]:
        bad.append((t, q, r["rel_error"]))
mg = [r for r in rows if r["test"] == "V3" and "160 V" in r["quantity"]]
bad += [(r["test"], r["quantity"], r["rel_error"]) for r in mg if abs(float(r["rel_error"])) > 1e-2]
print("failed:", bad if bad else "none")
sys.exit(1 if bad else 0)
