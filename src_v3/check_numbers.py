"""Regenerate the headline numbers quoted in the paper from data_v3/*.csv.
Run:  python src_v3/check_numbers.py   (needs only numpy and pandas)."""
import json, os, sys
import numpy as np, pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_v3")
rd = lambda f: pd.read_csv(os.path.join(D, f))
ok = True
def show(label, val, quoted=None, tol=0.06):
    global ok
    flag = ""
    if quoted is not None:
        good = abs(val - quoted) <= tol * abs(quoted)
        ok &= good
        flag = "  [quoted %g]  %s" % (quoted, "OK" if good else "MISMATCH")
    print("%-58s %.4g%s" % (label, val, flag))

# --- Sec. 3.5 : measured pinning, gated on state (VG=-3 V, VD=0.05 V)
m = rd("metals_pinning.csv"); g = m[(m.study == "gated") & (m.VG == -3)]
print("cases:", sorted(m.case.unique()))
def Ron(case, metal, geom):
    r = g[(g.case == case) & (g.metal == metal) & (g.geom == geom)]
    return r.R_tot.iloc[0] if len(r) else np.nan
direct = [c for c in m.case.unique() if c.startswith("direct")][0]
for c in m.case.unique():
    if "interlayer" in c or "oxide" in c: oxide = c
    if "calc" in c or "ideal" in c: calc = c
show("edge on-state R, direct Pt (Ohm cm)", Ron(direct, "Pt", "edge"), 5.2e3)
show("edge on-state R, direct Ti (Ohm cm)", Ron(direct, "Ti", "edge"), 1.2e4)
show("Ti/Pt ratio, direct", Ron(direct, "Ti", "edge") / Ron(direct, "Pt", "edge"), 2.4)
show("edge R, Pt, oxide interlayer (Ohm cm)", Ron(oxide, "Pt", "edge"))
show("improvement oxide vs direct Pt", Ron(direct, "Pt", "edge") / Ron(oxide, "Pt", "edge"), 700, 0.1)
show("edge R, Pt, calculated ideal (Ohm cm)", Ron(calc, "Pt", "edge"), 7.1)
show("top/edge ratio, direct Pt (phiBp=1.05 eV)", Ron(direct, "Pt", "top50") / Ron(direct, "Pt", "edge"))
print("non-converged / missing gated points (top):")
for c in m.case.unique():
    for me in m.metal.unique():
        if np.isnan(Ron(c, me, "top50")) and not np.isnan(Ron(c, me, "edge")):
            print("   top missing:", c, me)

# --- edge/top ratio versus hole barrier
b = rd("barrier_sweep.csv"); b = b[(b.study == "gated") & (b.VG == -3) & (b.tunneling != "off")]
print("top/edge on-state ratio vs barrier (paper: 4.1,14,35,51,74,150):")
for p in sorted(b.phiB_set.unique()):
    t = b[(b.geom == "top50") & (b.phiB_set == p)].R_tot; e = b[(b.geom == "edge") & (b.phiB_set == p)].R_tot
    if len(t) and len(e): print("   phiBp=%.1f eV  ratio=%.3g" % (p, t.iloc[0] / e.iloc[0]))

# --- Table 3 and energy balance
t = rd("thermal_gated.csv"); t = t[t.case == "full"]
print("Table 3 (full model):")
for _, r in t.iterrows():
    print("   %-6s VD=%.1f  I=%.3f A/cm  IV=%.3f  dTmax=%.1f K  dTsrc=%.1f  dTdrn=%.1f  Pout/IV-1=%+.2f%%" % (
        r.geom, r.V, r.I, r.P_el, r.Tmax - 300, r.T_source - 300, r.T_drain - 300, 100 * (r.P_out / r.P_el - 1)))

# --- Au/MoS2 benchmark, UQ
s = json.load(open(os.path.join(D, "derived_summary.json")))
f = s["mos2"]["fits"]["300.0"]
show("MoS2 Rc 300 K (Ohm um)", f["Rc_long_ohm_um"], 731, 0.01)
show("MoS2 LT 300 K (nm)", f["L_T_nm"], 34, 0.03)
show("MoS2 rho_c 300 K (1e-7 Ohm cm2)", f["rho_c"] * 1e7, 2.7, 0.03)
for k, v in s["uq"].items(): print("UQ %-10s p5 %.3g  median %.3g  p95 %.3g" % (k, v["p5"], v["median_R"], v["p95"]))
print("\nALL QUOTED NUMBERS REPRODUCED" if ok else "\nSOME QUOTED NUMBERS DO NOT MATCH")
sys.exit(0 if ok else 1)
