"""Summaries of the revision studies from data_v3/*.csv (numpy + pandas only).
Run after run_revision.py:  python analyze_revision.py"""
import os
import numpy as np
import pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_v3")
def rd(n):
    p = os.path.join(D, n + ".csv")
    return pd.read_csv(p) if os.path.exists(p) else None
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
PAR = ["dEV", "S", "S_edge_indep", "cnl_off", "d_vdW", "Phi_vdW", "mu_p", "rmu", "NA", "eps_r", "mh"]

def first_order(x, y, nb=8):
    """Binned first-order index Var(E[y|x])/Var(y) minus the same quantity for a shuffled x (noise floor)."""
    def idx(xx):
        q = pd.qcut(pd.Series(xx).rank(method="first"), nb, labels=False)
        g = pd.Series(y).groupby(q.values)
        return (g.mean().var(ddof=0) * 1.0) / np.var(y)
    rng = np.random.default_rng(0)
    floor = np.mean([idx(rng.permutation(x)) for _ in range(50)])
    return idx(x), floor

u = rd("uq_gated")
if u is not None:
    u = u[u.get("failed", pd.Series(0, index=u.index)).fillna(0) == 0]
    print("\n=== uq_gated: gated on state, Pt, VG=-3 V, VD=0.05 V (n samples per geometry/variant shown)")
    last = u[u.VG == -3.0]
    for (geom, var), g in last.groupby(["geom", "variant"]):
        y = np.log10(g.R_tot.clip(lower=1e-12))
        print(f"\n{geom}/{var}: n={len(g)}  log10 R_tot: p5={y.quantile(.05):.2f} med={y.median():.2f} p95={y.quantile(.95):.2f}"
              f"   tunnelling share median={g.tun_share.median():.3f}")
        rows = []
        for p in PAR:
            if p in g and g[p].nunique() > 3:
                x = g[p].values; si, fl = first_order(x, y.values)
                rows.append((p, pd.Series(x).rank().corr(y.reset_index(drop=True).rank()), si, fl))
        print(pd.DataFrame(rows, columns=["param", "spearman", "S1_binned", "noise_floor"]).round(3).to_string(index=False))
    sh = last[last.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot")
    if {"top50", "edge"} <= set(sh.columns):
        r = (sh["edge"] / sh["top50"]).dropna()
        print(f"\nedge/top ratio (shared pinning): median={r.median():.3g} p5={r.quantile(.05):.3g} p95={r.quantile(.95):.3g} "
              f"edge lower in {(r < 1).mean() * 100:.0f}% of {len(r)} samples")
    ie = last[last.variant == "indep_edge"].set_index("sample").R_tot
    if len(ie) and "top50" in sh:
        r2 = (ie / sh["top50"]).dropna()
        print(f"edge with independent S_edge vs top: edge lower in {(r2 < 1).mean() * 100:.0f}% of {len(r2)} samples")

t = rd("thermal_lit")
if t is not None:
    t = t[t.get("failed", 0).fillna(0) == 0] if "failed" in t else t
    t = t.assign(dT=t.Tmax - 300, balance_pct=100 * (t.P_out / t.P_el - 1))
    print("\n=== thermal_lit (Pt, VG=-3 V)")
    print(t[["case", "variant", "geom", "V", "I", "P_el", "dT", "balance_pct"]].round(3).to_string(index=False))

mc = rd("mesh_convergence")
if mc is not None:
    print("\n=== mesh_convergence: change relative to mf=1")
    key = [c for c in ["case", "geom", "V", "VG"] if c in mc]
    base = mc[mc.mf == 1.0].set_index(key)
    for mf in (2.0, 3.0):
        cur = mc[mc.mf == mf].set_index(key)
        d = pd.DataFrame({"dR_tot_%": 100 * (cur.R_tot / base.R_tot - 1), "dI_%": 100 * (cur.I / base.I - 1),
                          "dTmax_K": cur.Tmax - base.Tmax}).dropna(how="all")
        print(f"\nmf={mf}"); print(d.round(3).to_string())

th = rd("thickness_lit")
if th is not None:
    print("\n=== thickness_lit: gated on-state R_tot (Ohm cm) vs layers")
    th = th[th.VG == -3.0] if "VG" in th else th
    cols = [c for c in ["case", "dEv_frac"] if c in th]
    for k, g in th.groupby(["study"] + cols, dropna=False):
        pv = g.pivot_table(index="layers", columns="geom", values="R_tot")
        best = {c: int(pv[c].idxmin()) for c in pv}
        print(f"\n{k}  best layer count: {best}"); print(pv.round(4).to_string())

e = rd("inse_electron")
if e is not None:
    print("\n=== inse_electron: Arrhenius slope of ln(I/T^2) vs 1/T -> effective barrier (eV) per VG")
    e = e[e.get("failed", 0).fillna(0) == 0] if "failed" in e else e
    kB = 8.617e-5
    out = []
    for (metal, geom, vg), g in e.groupby(["metal", "geom", "VG"]):
        if len(g) >= 3:
            sl = np.polyfit(1 / (kB * g.T0.values), np.log(g.I.values / g.T0.values ** 2), 1)[0]
            out.append((metal, geom, vg, -sl))
    print(pd.DataFrame(out, columns=["metal", "geom", "VG", "phi_eff_eV"]).round(3).to_string(index=False))
    print("Compare the flat-band barrier (VG where phi_eff stops decreasing) with the S = 0.05 / 0.3 eV values of Ref. [11].")

c = rd("chen2019")
if c is not None:
    c = c[c["failed"].fillna(0) == 0] if "failed" in c else c
    kB = 8.617e-5
    a = c[c.part == "A_arrhenius"]
    print("\n=== chen2019 part A: Arrhenius barrier (300-400 K) per V_G, top contacts, 15 nm n-InSe, 300 nm SiO2")
    rows = []
    for (case, metal, vg), g in a.groupby(["case", "metal", "VG"]):
        if len(g) >= 3:
            rows.append((case, metal, g.Phi_M.iloc[0], vg, -np.polyfit(1 / (kB * g.T0.values), np.log(g.I.values / g.T0.values ** 2), 1)[0]))
    A = pd.DataFrame(rows, columns=["case", "metal", "Phi_M", "VG", "phi_B"])
    print(A.pivot_table(index=["case", "metal"], columns="VG", values="phi_B").round(3).to_string())
    print("apparent pinning factor (slope of phi_B vs Phi_M across metals), reported: control 0.05, OML 0.5")
    for (case, vg), g in A.groupby(["case", "VG"]):
        if len(g) >= 3:
            print(f"   {case:20s} VG={vg:+5.0f} V  S_app={np.polyfit(g.Phi_M, g.phi_B, 1)[0]:+.3f}")
    print("reported: In with oxidized monolayer, 65 meV")
    b = c[c.part == "B_oncurrent"]
    if len(b):
        print("\n=== chen2019 part B: In contact, V_D=1 V; reported 103 uA/um (1.03 A/cm) for OML at V_G=80 V, 3-6x control")
        pv = b.pivot_table(index=["Lch_nm", "VG"], columns="case", values="I")
        pv["ratio_OML/control"] = pv["oxidized monolayer"] / pv["control"]
        print(pv.to_string(float_format=lambda x: f"{x:.4g}"))
