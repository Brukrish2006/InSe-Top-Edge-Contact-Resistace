"""Check every number quoted in the revised sections of the paper against data_v3/*.csv.

Run:  python src_v3/check_paper.py          (numpy + pandas only)
Each line prints the value recomputed from the CSVs, the value quoted in the
paper and OK / MISMATCH.  Sections whose CSV is missing are skipped.
check_numbers.py covers the baseline sections (MoS2 benchmark, measured
pinning, Table of self-heating); this script covers the revision studies.
"""
import os
import sys
import numpy as np
import pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_v3")
N_OK = N_BAD = 0
CASE_PHI = {"phiB=0.2": 0.2, "unpinned": 0.25, "interlayer (measured)": 0.52, "phiB=0.6": 0.6,
            "phiB=0.8": 0.8, "phiB=1.0": 1.0, "direct (measured)": 1.05}


def rd(name):
    p = os.path.join(D, name + ".csv")
    if not os.path.exists(p):
        print(f"\n[skip] {name}.csv not found")
        return None
    d = pd.read_csv(p)
    if "failed" in d:
        d = d[d["failed"].fillna(0) == 0]
    return d


def head(t):
    print("\n=== " + t)


def q(label, val, quoted=None, tol=0.06, where=""):
    """Compare val with the quoted value (relative tolerance tol; tol<0 means absolute |tol|)."""
    global N_OK, N_BAD
    s = f"  {label:62s} {val:11.4g}"
    if quoted is not None:
        good = (abs(val - quoted) <= -tol) if tol < 0 else (abs(val - quoted) <= tol * abs(quoted))
        N_OK += good
        N_BAD += (not good)
        s += f"   quoted {quoted:<9g} {'OK' if good else 'MISMATCH'}  {where}"
    print(s)


def ratio(df, idx):
    pv = df.pivot_table(index=idx, columns="geom", values="R_tot")
    return pv["top50"] / pv["edge"]


def crossing(x, y, logx=False):
    """x where log(y) crosses 0 (y = top/edge ratio); log(y) is interpolated linearly in x
    (in log x if logx); nan if there is no crossing.  The paper uses linear x."""
    x = np.asarray(x, float)
    ly = np.log(np.asarray(y, float))
    xx = np.log(x) if logx else x
    for i in range(len(x) - 1):
        if ly[i] == 0:
            return x[i]
        if ly[i] * ly[i + 1] < 0:
            xc = xx[i] - ly[i] * (xx[i + 1] - xx[i]) / (ly[i + 1] - ly[i])
            return float(np.exp(xc)) if logx else float(xc)
    return float("nan")


# ---------------------------------------------------------------- Sec. 3.2: oblique paths, isotropic mass
ob = rd("oblique")
if ob is not None:
    head("Sec. 3.2  most favourable straight path, light isotropic mass (oblique.csv)")
    ob = ob.assign(oblique=ob.oblique.astype(int))
    r = ratio(ob, ["case", "oblique"]).unstack("oblique")
    pv = ob.pivot_table(index=["case", "oblique"], columns="geom", values="R_tot").unstack("oblique")
    gain_top = pv["top50"][0] / pv["top50"][1]
    gain_edge = 100 * (1 - pv["edge"][1] / pv["edge"][0])
    q("max fall of top R with best path (x)", gain_top.max(), 2.4, where="Sec 3.2")
    q("min fall of edge R (%)", gain_edge.min(), 6, tol=-1.5, where="Sec 3.2")
    q("max fall of edge R (%)", gain_edge.max(), 31, tol=-1.5, where="Sec 3.2")
    q("ratio phiB=1.0, normal", r.loc["phiB=1.0", 0], 154, 0.02)
    q("ratio phiB=1.0, best path", r.loc["phiB=1.0", 1], 97, 0.02)
    q("ratio direct Pt, normal", r.loc["direct (measured)", 0], 170, 0.02)
    q("ratio direct Pt, best path", r.loc["direct (measured)", 1], 90, 0.03)
    q("ratio phiB=0.6, normal", r.loc["phiB=0.6", 0], 51, 0.03)
    q("ratio phiB=0.6, best path", r.loc["phiB=0.6", 1], 70, 0.03)

# ---------------------------------------------------------------- Sec. 3.3: heavier isotropic mass
mh = rd("mh_sweep")
if mh is not None:
    head("Sec. 3.3  heavier isotropic mass 0.5-2 m0, parabolic (mh_sweep.csv, tunnelling mass alone or with DOS mass)")
    m = mh[mh.variant.isin(["both", "tun"]) & (mh.mh >= 0.5)]
    r = ratio(m, ["case", "variant", "mh"])
    for case, lo, hi in (("direct (measured)", 670, 1800), ("interlayer (measured)", 310, 750), ("unpinned", 57, 98)):
        q(f"min ratio {case}", r[case].min(), lo, 0.04)
        q(f"max ratio {case}", r[case].max(), hi, 0.04)

# ---------------------------------------------------------------- Sec. 3.3: anisotropic mass, Fig. 17
km = rd("kane_map")
if km is not None:
    head("Sec. 3.3  mass tensor m_c=0.17 (kane_map.csv)")
    r = ratio(km, ["variant", "case", "m_par"])
    ko = r["kane_oblique"]
    q("edge better, phiB=0.2 at m_par=0.9 (top/edge)", ko["phiB=0.2"][0.9], 11.5, 0.03)
    q("edge better, unpinned at 0.9", ko["unpinned"][0.9], 12.5, 0.03)
    q("edge better, interlayer at 0.9", ko["interlayer (measured)"][0.9], 2.1, 0.04)
    q("ratio phiB=0.6 at 0.9", ko["phiB=0.6"][0.9], 0.9, 0.06)
    q("top better, phiB=1.0 at 0.9", 1 / ko["phiB=1.0"][0.9], 20, 0.04)
    q("top better, direct at 0.9", 1 / ko["direct (measured)"][0.9], 27, 0.04)
    q("edge wins at every m_par for barriers <= 0.52 (min edge advantage)",
      min(ko[c].min() for c in ("phiB=0.2", "unpinned", "interlayer (measured)")), None)
    quoted = {"parabolic": {"interlayer (measured)": 0.77, "phiB=0.6": 0.63, "phiB=1.0": 0.43, "direct (measured)": 0.42},
              "kane": {"interlayer (measured)": 0.67, "phiB=0.6": 0.54, "phiB=1.0": 0.41, "direct (measured)": 0.40},
              "kane_oblique": {"phiB=0.6": 0.86, "phiB=1.0": 0.44, "direct (measured)": 0.42}}
    for var, d in quoted.items():
        for case, qv in d.items():
            s = r[var][case].sort_index()
            q(f"crossover m_par, {var}, {case}", crossing(s.index, s.values), qv, tol=-0.011)
    s = r["kane_oblique"]["interlayer (measured)"]
    q("kane_oblique interlayer: no crossover up to 1.2 (ratio at 1.2 > 1)", s[1.2], None)
    # path geometry: kane vs kane_oblique at heaviest m_par
    f = (r["kane"] / r["kane_oblique"]).xs(1.2, level="m_par")
    q("max change of ratio from best path at m_par=1.2 (x)", max(f.max(), 1 / f.min()), 40, 0.08, "Limitations")
    pv = km[km.variant == "kane_oblique"].pivot_table(index=["case", "m_par"], columns="geom", values="R_tot")
    q("depinning gain, top, m_par=0.9 (direct/interlayer)",
      pv.loc[("direct (measured)", 0.9), "top50"] / pv.loc[("interlayer (measured)", 0.9), "top50"], 290, 0.05)
    q("depinning gain, edge, m_par=0.9",
      pv.loc[("direct (measured)", 0.9), "edge"] / pv.loc[("interlayer (measured)", 0.9), "edge"], 1.7e4, 0.05)
    q("direct-pinned min R at 0.9 (Ohm cm), > 3e4",
      pv.loc[("direct (measured)", 0.9)].min(), 3e4, tol=-1e9)

mc = rd("mass_controls")
if mc is not None:
    head("Sec. 3.3  controls at 0.9 m0, two-band (mass_controls.csv)")
    r = ratio(mc, ["variant", "case"])
    q("iso: min top/edge", r["iso"].min(), 58, 0.04)
    q("iso: max top/edge", r["iso"].max(), 850, 0.04)
    q("swapped: min top/edge", r["swapped"].min(), 180, 0.05)
    q("swapped: max top/edge", r["swapped"].max(), 9e6, 0.06)

et = rd("edge_T")
if et is not None and km is not None:
    head("Sec. 3.3  edge-interface transmission at which top = edge (edge_T.csv, kane normal paths)")
    top = km[(km.variant == "kane") & (km.m_par == 0.9) & (km.geom == "top50")].set_index("case").R_tot
    for case, qv in (("phiB=0.2", 0.17), ("unpinned", 0.22)):
        g = et[et.case == case].sort_values("T_edge")
        x = np.log(g.T_edge.values); y = np.log(g.R_tot.values / top[case])
        q(f"T_edge threshold {case}", float(np.exp(np.interp(0, y[::-1], x[::-1]))), qv, 0.1)

# ---------------------------------------------------------------- Sec. 3.4: phase diagram
def phase_section(name, quoted_thr):
    mp = rd(name)
    if mp is None:
        return
    head(f"Sec. 3.4  geometry phase diagram ({name}.csv)  n={len(mp)} converged runs")
    r = ratio(mp, ["case", "m_c", "m_par"])
    for case in r.index.get_level_values(0).unique():
        thr = []
        for mc_ in sorted(r[case].index.get_level_values(0).unique()):
            s = r[case][mc_].sort_index()
            x = crossing(s.index, s.values)
            if np.isfinite(x):
                thr.append(x / mc_)
            print(f"    {case:22s} m_c={mc_:<6g} crossover m_par={x:.3g}  ratio m_par/m_c={x / mc_:.3g}")
        if thr and case in quoted_thr:
            q(f"min threshold ratio, {case}", min(thr), quoted_thr[case][0], 0.08)
            q(f"max threshold ratio, {case}", max(thr), quoted_thr[case][1], 0.08)
    if 2.5 in r.index.get_level_values(1):
        hv = r.xs(2.5, level="m_c")
        q("m_c=2.5: min edge advantage", hv.min(), 100, 0.15)
        q("m_c=2.5: max edge advantage", hv.max(), 5000, 0.15)
    for mc_, vals in ((0.17, (27, 12.5, 2.1)), (0.078, (13, 5.2, 1.0))):
        if (mc_ in r.index.get_level_values(1)) and (0.9 in r.index.get_level_values(2)):
            q(f"m_c={mc_}: top better, direct", 1 / r["direct (measured)"][mc_][0.9], vals[0], 0.05)
            q(f"m_c={mc_}: edge better, unpinned", r["unpinned"][mc_][0.9], vals[1], 0.05)
            q(f"m_c={mc_}: edge better, interlayer", r["interlayer (measured)"][mc_][0.9], vals[2], 0.06)
    return r


phase_section("mass_phase", {"direct (measured)": (2.1, 2.6), "phiB=0.8": (2.9, 3.3)})

_mf = rd("mass_phase_fine")
if _mf is not None:
    head("Sec. 3.4  WKB form ln(Rtop/Redge) = A + B sqrt(m_c) - C sqrt(m_par) fitted to mass_phase_fine")
    rr = np.log(ratio(_mf, ["case", "m_c", "m_par"]))
    for case, (qr2lo, qr2hi), qthr, qI in (("direct (measured)", (0.90, 0.97), 2.3, 1.5), ("phiB=0.8", (0.90, 0.97), 4.5, 2.1),
                                             ("interlayer (measured)", (0.90, 0.97), 9.0, 3.0), ("unpinned", (0.90, 0.97), 32, 5.7)):
        d = rr[case].reset_index(); d.columns = ["mc", "mp", "y"]
        X = np.c_[np.ones(len(d)), np.sqrt(d.mc), np.sqrt(d.mp)]
        b, *_ = np.linalg.lstsq(X, d.y, rcond=None); pred = X @ b
        R2 = 1 - ((d.y - pred) ** 2).sum() / ((d.y - d.y.mean()) ** 2).sum()
        q(f"R2 of WKB fit, {case} (within {qr2lo}-{qr2hi})", R2, (qr2lo + qr2hi) / 2, tol=-0.0375)
        q(f"(B/C)^2 threshold, {case}", (b[1] / -b[2]) ** 2, qthr, 0.04)
        q(f"I_top/I_edge = B/C, {case}", b[1] / -b[2], qI, 0.04)
phase_section("mass_phase_fine", {"direct (measured)": (2.1, 2.8), "phiB=0.8": (2.7, 3.0), "interlayer (measured)": (11.3, 11.8)})

# ---------------------------------------------------------------- Sec. 3.6: Monte Carlo with anisotropic mass
def mc_section(name, quoted):
    uk = rd(name)
    if uk is None:
        return
    head(f"Sec. 3.6  Monte Carlo, anisotropic mass ({name}.csv)")
    raw = pd.read_csv(os.path.join(D, name + ".csv"))
    if "failed" in raw:
        q("failed runs", raw["failed"].fillna(0).sum(), quoted.get("failed"), tol=-0.5)
    sh = uk[uk.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
    r = sh["edge"] / sh["top50"]
    q("converged paired samples (shared)", len(r), quoted.get("n"), tol=-0.5)
    q("edge lower, shared (%)", 100 * (r < 1).mean(), quoted.get("shared"), tol=-0.6)
    for var in ("indep_edge", "indep_cnl"):
        e = uk[uk.variant == var].set_index("sample").R_tot
        j = e.index.intersection(sh.index)
        q(f"edge lower, {var} (%)", 100 * (e[j] / sh.loc[j, "top50"] < 1).mean(), quoted.get(var), tol=-0.6)
    q("median edge/top", r.median(), quoted.get("median"), 0.05)
    q("p5 edge/top", r.quantile(0.05), quoted.get("p5"), 0.06)
    q("p95 edge/top", r.quantile(0.95), quoted.get("p95"), 0.08)
    mi = uk[uk.variant == "shared"].groupby("sample").m_in.first().reindex(r.index)
    for (lo, hi), qv in zip(((0.3, 0.5), (0.5, 0.8), (0.8, 1.5)), quoted.get("bins", (None,) * 3)):
        sel = (mi >= lo) & (mi < hi)
        q(f"edge lower, m_par {lo}-{hi} (%) [n={int(sel.sum())}]", 100 * (r[sel] < 1).mean(), qv, tol=-0.6)


mc_section("uq_kane", dict(failed=4, n=239, shared=20, indep_edge=19, indep_cnl=22, median=13, p5=0.25, p95=2700,
                           bins=(57, 11, 0)))
mc_section("uq_obl", dict(failed=0, n=240, shared=29, indep_edge=28, indep_cnl=33, median=3.8, p5=0.11, p95=82,
                          bins=(66, 24, 7)))

# ---------------------------------------------------------------- Sec. 3.7: thickness with anisotropic mass
tk = rd("thickness_kane")
if tk is not None:
    head("Sec. 3.7  thickness, m_par=0.9, two-band, normal paths (thickness_kane.csv)")
    pv = tk.pivot_table(index=["case", "layers"], columns="geom", values="R_tot")
    il = pv.loc["interlayer (measured)"]; di = pv.loc["direct (measured)"]
    q("interlayer top R, 4 layers", il.loc[4, "top50"], 126, 0.03)
    q("interlayer top R, 12 layers", il.loc[12, "top50"], 1.9e4, 0.04)
    q("interlayer edge R, 4 layers", il.loc[4, "edge"], 1130, 0.03)
    q("interlayer edge R, 12 layers", il.loc[12, "edge"], 101, 0.03)
    q("interlayer top better at 4 (edge/top)", il.loc[4, "edge"] / il.loc[4, "top50"], 9, 0.06)
    q("interlayer top better at 6", il.loc[6, "edge"] / il.loc[6, "top50"], 2.5, 0.06)
    q("interlayer edge better at 8", il.loc[8, "top50"] / il.loc[8, "edge"], 4.6, 0.05)
    q("interlayer edge better at 12", il.loc[12, "top50"] / il.loc[12, "edge"], 190, 0.05)
    q("direct top better at 6", di.loc[6, "edge"] / di.loc[6, "top50"], 250, 0.1)
    q("direct top better at 8", di.loc[8, "edge"] / di.loc[8, "top50"], 250, 0.1)
    q("direct ratio at 12 (equal)", di.loc[12, "edge"] / di.loc[12, "top50"], 1.0, 0.25)

to = rd("thickness_obl")
if to is not None:
    head("thickness with best straight paths (thickness_obl.csv)")
    pv = to.pivot_table(index=["case", "layers"], columns="geom", values="R_tot")
    il = pv.loc["interlayer (measured)"]; di = pv.loc["direct (measured)"]
    q("interlayer top R, 4 layers", il.loc[4, "top50"], 120, 0.02)
    q("interlayer top R, 12 layers", il.loc[12, "top50"], 1.8e4, 0.03)
    q("interlayer edge R, 3 layers", il.loc[3, "edge"], 1060, 0.01)
    q("interlayer edge R, 12 layers", il.loc[12, "edge"], 21, 0.03)
    q("interlayer top better at 4 (edge/top)", il.loc[4, "edge"] / il.loc[4, "top50"], 1.8, 0.03)
    q("interlayer edge better at 3", il.loc[3, "top50"] / il.loc[3, "edge"], 1.5, 0.03)
    q("interlayer edge better at 6", il.loc[6, "top50"] / il.loc[6, "edge"], 2.0, 0.03)
    q("interlayer edge better at 8", il.loc[8, "top50"] / il.loc[8, "edge"], 22, 0.03)
    q("interlayer edge better at 12", il.loc[12, "top50"] / il.loc[12, "edge"], 860, 0.01)
    q("direct min top advantage 4-8 layers", min(di.loc[n, "edge"] / di.loc[n, "top50"] for n in (4, 6, 8)), 12, 0.03)
    q("direct max top advantage 4-8 layers", max(di.loc[n, "edge"] / di.loc[n, "top50"] for n in (4, 6, 8)), 27, 0.03)
    q("direct edge better at 12", di.loc[12, "top50"] / di.loc[12, "edge"], 12, 0.03)

# ---------------------------------------------------------------- self-heating transients and Peltier constant
hs = rd("heat_step")
if hs is not None:
    head("Sec. 3.8  thermal response time (heat_step.csv)")
    t63 = []
    for _, g in hs.groupby(["case", "tunneling", "geom"]):
        g = g.sort_values("t"); dT = g.Tmax.values - 300
        if dT[-1] > 0:
            t63.append(np.interp(0.632 * dT[-1], dT, g.t.values) * 1e9)
    q("min t63 (ns)", min(t63), 0.19, 0.05)
    q("max t63 (ns)", max(t63), 0.25, 0.05)
tc = rd("thermo_const")
if tc is not None:
    head("Sec. 3.8  contact Peltier constant 2 vs 1.5 (thermo_const.csv)")
    g = tc.pivot_table(index=["geom", "V"], columns="pelt_const", values="Tmax")
    q("max |dTmax| between constants (K)", (g[2.0] - g[1.5]).abs().max(), 3, tol=-0.6)

# ---------------------------------------------------------------- Sec. 3.3: metals with the anisotropic model
am = rd("aniso_metals")
if am is not None:
    head("Sec. 3.3  metals x pinning, anisotropic model (aniso_metals.csv)")
    raw = pd.read_csv(os.path.join(D, "aniso_metals.csv"))
    q("failed runs", raw["failed"].fillna(0).sum() if "failed" in raw else 0, 4, tol=-0.5)
    b = am[am.study == "barrier"]; rb = ratio(b, ["phiB_set"])
    for pb, qv in ((0.0, 4.2), (0.2, 11.5), (0.4, 6.7), (0.6, 0.90), (0.8, 0.12), (1.0, 0.051)):
        q(f"barrier sweep ratio at {pb} eV", rb[pb], qv, 0.04)
    mm = am[am.study == "metals"]; pv = mm.pivot_table(index=["case", "metal"], columns="geom", values="R_tot")
    r = pv["top50"] / pv["edge"]
    d = pv.loc["direct (measured)"]
    q("direct: edge R Ti/Pt", d.loc["Ti", "edge"] / d.loc["Pt", "edge"], 2.9, 0.04)
    q("direct: top R Ti/Pt", d.loc["Ti", "top50"] / d.loc["Pt", "top50"], 1.8, 0.04)
    q("direct: min top advantage", 1 / r["direct (measured)"].max(), 27, 0.04)
    q("direct: max top advantage", 1 / r["direct (measured)"].min(), 44, 0.04)
    il = r["interlayer (measured)"]
    q("interlayer Pt edge advantage", il["Pt"], 2.1, 0.04)
    q("interlayer Pd/Au min top advantage", 1 / max(il["Pd"], il["Au"]), 7.5, 0.04)
    q("interlayer Pd/Au max top advantage", 1 / min(il["Pd"], il["Au"]), 8.0, 0.05)
    q("interlayer Cr top advantage", 1 / il["Cr"], 39, 0.03)
    q("interlayer Ti top advantage", 1 / il["Ti"], 67, 0.03)
    dft = r["ideal MIGS (DFT)"]
    q("DFT Pt edge advantage", dft["Pt"], 2.1, 0.04)
    q("DFT Pd/Au min top advantage", 1 / max(dft["Pd"], dft["Au"]), 1.3, 0.04)
    q("DFT Pd/Au max top advantage", 1 / min(dft["Pd"], dft["Au"]), 1.4, 0.04)
    q("DFT Cr top advantage", 1 / dft["Cr"], 4.3, 0.05)
    q("DFT Ti top advantage", 1 / dft["Ti"], 5.8, 0.05)
    q("unpinned Pt edge advantage", r["unpinned"]["Pt"], 12.5, 0.03)

ps_ = rd("physics_sens")
if ps_ is not None:
    head("Sec. 3.3  two-band only along c (physics_sens.csv, variant kane_c)")
    r = ratio(ps_[ps_.variant == "kane_c"], ["case"])
    q("kane_c interlayer", r["interlayer (measured)"], 1.3, 0.04)
    q("kane_c phiB=0.6", r["phiB=0.6"], 0.47, 0.04)
    q("kane_c min (1.0-1.05 eV)", min(r["phiB=1.0"], r["direct (measured)"]), 0.0039, 0.05)
    q("kane_c max (1.0-1.05 eV)", max(r["phiB=1.0"], r["direct (measured)"]), 0.0068, 0.05)
    q("kane_c phiB=0.2", r["phiB=0.2"], 11.1, 0.03)
    q("kane_c unpinned", r["unpinned"], 11.6, 0.03)

ath = rd("aniso_thermal")
if ath is not None:
    head("Sec. 3.8  self-heating, anisotropic model (aniso_thermal.csv)")
    a = ath.assign(dT=ath.Tmax - 300).set_index(["case", "geom", "V"])
    q("reference edge I at 1 V (A/cm)", a.loc[("reference", "edge", 1.0), "I"], 4.1, 0.02)
    q("reference top I at 1 V (A/cm)", a.loc[("reference", "top50", 1.0), "I"], 2.0, 0.03)
    q("reference edge dT at 1 V (K)", a.loc[("reference", "edge", 1.0), "dT"], 172, 0.01)
    q("reference top dT at 1 V (K)", a.loc[("reference", "top50", 1.0), "dT"], 34, 0.02)
    q("interlayer edge dT (K)", a.loc[("interlayer (measured)", "edge", 1.0), "dT"], 14, 0.05)
    q("interlayer top dT (K)", a.loc[("interlayer (measured)", "top50", 1.0), "dT"], 1.0, 0.06)
    q("unpinned edge dT (K)", a.loc[("unpinned", "edge", 1.0), "dT"], 126, 0.01)
    q("unpinned top dT (K)", a.loc[("unpinned", "top50", 1.0), "dT"], 10, 0.05)
    q("energy balance edge 1 V (%)", 100 * (a.loc[("reference", "edge", 1.0), "P_out"] / a.loc[("reference", "edge", 1.0), "P_el"] - 1), 1.4, 0.05)

imf = rd("image_force"); hyb = rd("hybrid")
if imf is not None and hyb is not None:
    head("Sec. 3.3  image-force lowering (image_force.csv vs hybrid.csv aniso references)")
    r = ratio(imf, ["case"])
    base = hyb[(hyb.variant == "aniso") & (hyb.geom != "hybrid50")].pivot_table(index="case", columns="geom", values="R_tot")
    pi = imf.pivot_table(index="case", columns="geom", values="R_tot")
    q("image interlayer ratio", r["interlayer (measured)"], 4.8, 0.03)
    q("image phiB=0.6 ratio", r["phiB=0.6"], 2.0, 0.03)
    q("image min ratio 1.0-1.05", min(r["phiB=1.0"], r["direct (measured)"]), 0.055, 0.03)
    q("image max ratio 1.0-1.05", max(r["phiB=1.0"], r["direct (measured)"]), 0.070, 0.03)
    q("image min ratio 0.2-0.25", min(r["phiB=0.2"], r["unpinned"]), 26, 0.03)
    q("image max ratio 0.2-0.25", max(r["phiB=0.2"], r["unpinned"]), 38, 0.03)
    hi = ["interlayer (measured)", "phiB=0.6", "phiB=1.0", "direct (measured)"]
    fe = base.loc[hi, "edge"] / pi.loc[hi, "edge"]; ft = base.loc[hi, "top50"] / pi.loc[hi, "top50"]
    q("edge R fall min (>=0.52 eV)", fe.min(), 1.9, 0.04); q("edge R fall max", fe.max(), 4.6, 0.03)
    q("top R fall min", ft.min(), 1.3, 0.04); q("top R fall max", ft.max(), 2.0, 0.03)

ho = rd("hybrid_overlap")
if ho is not None and hyb is not None:
    head("Sec. 3.4  hybrid contacts (hybrid.csv, hybrid_overlap.csv)")
    a = hyb[hyb.variant == "aniso"].pivot_table(index="case", columns="geom", values="R_tot")
    o = ho[ho.variant == "aniso"].pivot_table(index="case", columns="Lc_set", values="R_tot")
    hi = ["phiB=0.6", "phiB=1.0", "direct (measured)"]; lo = ["phiB=0.2", "unpinned", "interlayer (measured)"]
    q("hybrid50/top min (0.6-1.05)", (a.loc[hi, "hybrid50"] / a.loc[hi, "top50"]).min(), 1.02, 0.01)
    q("hybrid50/top max (0.6-1.05)", (a.loc[hi, "hybrid50"] / a.loc[hi, "top50"]).max(), 1.04, 0.01)
    q("hybrid50/edge min (0.2-0.52)", (a.loc[lo, "hybrid50"] / a.loc[lo, "edge"]).min(), 1.6, 0.03)
    q("hybrid50/edge max (0.2-0.52)", (a.loc[lo, "hybrid50"] / a.loc[lo, "edge"]).max(), 1.9, 0.03)
    q("hybrid2/edge min (0.2-0.25)", (o.loc[["phiB=0.2", "unpinned"], 2] / a.loc[["phiB=0.2", "unpinned"], "edge"]).min(), 1.13, 0.01)
    q("hybrid2/edge max (0.2-0.25)", (o.loc[["phiB=0.2", "unpinned"], 2] / a.loc[["phiB=0.2", "unpinned"], "edge"]).max(), 1.16, 0.01)
    hh = ["phiB=1.0", "direct (measured)"]
    q("hybrid2/edge min (1.0-1.05)", (o.loc[hh, 2] / a.loc[hh, "edge"]).min(), 2.7, 0.03)
    q("hybrid2/edge max (1.0-1.05)", (o.loc[hh, 2] / a.loc[hh, "edge"]).max(), 2.8, 0.03)
    q("hybrid2/top min (1.0-1.05)", (o.loc[hh, 2] / a.loc[hh, "top50"]).min(), 54, 0.02)
    q("hybrid2/top max (1.0-1.05)", (o.loc[hh, 2] / a.loc[hh, "top50"]).max(), 78, 0.02)
    best = pd.concat([a["edge"], a["top50"]], axis=1).min(axis=1)
    allh = pd.concat([o, a["hybrid50"].rename(50)], axis=1)
    worst = (allh.div(best, axis=0)).max(axis=0)
    q("best overlap: worst-case hybrid/best-geometry (must be > 1.5)", worst.min(), 1.5, tol=-1e9)
    q("any hybrid better than best geometry? (min ratio, must be > 1)", allh.div(best, axis=0).min().min(), 1.0, tol=-1e9)

# ---------------------------------------------------------------- revision batch rev15 (printed for the text)
for name, idx in (("physics_sens", ["case", "variant"]), ("hybrid", ["case", "variant"])):
    d = rd(name)
    if d is not None:
        head(f"{name}.csv")
        print(d.pivot_table(index=idx, columns="geom", values="R_tot").to_string(float_format=lambda x: f"{x:.4g}"))

print(f"\n{N_OK} quoted numbers reproduced, {N_BAD} mismatches")
sys.exit(1 if N_BAD else 0)
