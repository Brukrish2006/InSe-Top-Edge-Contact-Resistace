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
c2 = rd("chen2019_300K")
if c is not None and c2 is not None:
    c = pd.concat([c, c2], ignore_index=True)
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


cm = rd("chen2019_mu")
if cm is not None:
    cm = cm[cm["failed"].fillna(0) == 0] if "failed" in cm else cm
    if c is not None:
        b0 = c[(c.part == "B_oncurrent")].copy(); b0["mu_n"] = 1000.0
        cm = pd.concat([cm, b0], ignore_index=True)
    print("\n=== chen2019_mu: In contact, V_D=1 V, on-current (uA/um) vs mu_n and L_ch; reported OML 103 uA/um, ratio 3-6x")
    cm = cm.assign(I_uA_um=cm.I * 100.0)
    pv = cm[cm.VG == 80.0].pivot_table(index=["Lch_nm", "mu_n"], columns="case", values="I_uA_um")
    pv["ratio"] = pv["oxidized monolayer"] / pv["control"]
    print(pv.to_string(float_format=lambda x: f"{x:.4g}"))


mu_ = rd("mesh_ungated")
if mu_ is not None:
    mu_ = mu_[mu_["failed"].fillna(0) == 0] if "failed" in mu_ else mu_
    print("\n=== mesh_ungated: Rtot (Ohm cm) at 0.05 V, ungated, reference alignment, vs mesh factor")
    pv = mu_.pivot_table(index=["geom", "Phi_M"], columns="mf", values="R_tot")
    print(pv.to_string(float_format=lambda x: f"{x:.5g}"))
    rel = pv.div(pv[1.0], axis=0) - 1
    print("max |change| relative to factor 1 (%):", round(100 * rel.abs().max().max(), 2))
    print("max |change| between 0.5x and 3x (%):", round(100 * ((pv[3.0] / pv[0.5] - 1).abs().max()), 2))


tr = rd("top_retry")
if tr is not None:
    print("\n=== top_retry: previously non-converged gated on-state cases (V_G=-3 V, V_D=0.05 V)")
    if "failed" in tr:
        print("failed runs:", int(tr["failed"].fillna(0).sum()), "of", len(tr))
    ok_ = tr[tr["failed"].fillna(0) == 0] if "failed" in tr else tr
    if len(ok_):
        print(ok_[["geom", "case", "metal", "variant", "phiBp", "R_tot", "tun_share"]].round(4).to_string(index=False))

mh = rd("mh_sweep")
if mh is not None:
    print("\n=== mh_sweep: gated on state (Pt, V_G=-3 V, V_D=0.05 V) versus hole mass")
    if "failed" in mh:
        print("failed runs:", int(mh["failed"].fillna(0).sum()), "of", len(mh))
        mh = mh[mh["failed"].fillna(0) == 0]
    if len(mh):
        pv = mh.pivot_table(index=["case", "variant", "mh"], columns="geom", values="R_tot")
        if "top50" in pv and "edge" in pv:
            pv["top/edge"] = pv["top50"] / pv["edge"]
        print(pv.to_string(float_format=lambda x: f"{x:.4g}"))

ob = rd("oblique")
if ob is not None:
    print("\n=== oblique: gated on state (V_G=-3 V), normal-path vs minimum-action oblique tunnelling")
    if "failed" in ob:
        print("failed runs:", int(ob["failed"].fillna(0).sum()), "of", len(ob))
        ob = ob[ob["failed"].fillna(0) == 0]
    if len(ob):
        pv = ob.pivot_table(index=["case", "oblique"], columns="geom", values="R_tot")
        if "top50" in pv and "edge" in pv:
            pv["top/edge"] = pv["top50"] / pv["edge"]
        print(pv.to_string(float_format=lambda x: f"{x:.4g}"))

uc = rd("uq_cnl")
ug = rd("uq_gated")
if uc is not None and ug is not None:
    print("\n=== uq_cnl: edge with independent pinning level vs top (shared) from uq_gated")
    okc = uc[uc["failed"].fillna(0) == 0] if "failed" in uc else uc
    okg = ug[ug["failed"].fillna(0) == 0] if "failed" in ug else ug
    okg = okg[(okg.variant == "shared") & (okg.geom == "top50")]
    if "VG" in okg:
        okg = okg[np.isclose(okg.VG, -3.0)]
    if "VG" in okc:
        okc = okc[np.isclose(okc.VG, -3.0)]
    top = okg.set_index("sample").R_tot
    for var, g in okc.groupby("variant"):
        e = g.set_index("sample").R_tot
        j = e.index.intersection(top.index)
        r = (e[j] / top[j])
        print(f"{var}: n={len(j)}  edge<top in {100 * (r < 1).mean():.1f}%  median edge/top={r.median():.3g}  "
              f"5-95%: {r.quantile(0.05):.3g}-{r.quantile(0.95):.3g}")
    print("failed edge runs:", int(uc["failed"].fillna(0).sum()) if "failed" in uc else 0)

tc = rd("thermo_const")
if tc is not None:
    print("\n=== thermo_const: contact Peltier constant 2 vs 1.5 (Table-3 point)")
    cols = [c_ for c_ in ("geom", "pelt_const", "V", "I", "Tmax", "P_el", "P_out", "T_source", "T_drain") if c_ in tc]
    t2 = tc[cols].copy()
    if "Tmax" in t2:
        t2["dT"] = t2["Tmax"] - 300
    if "P_out" in t2 and "P_el" in t2:
        t2["balance_%"] = 100 * (t2["P_out"] / t2["P_el"] - 1)
    print(t2.round(4).to_string(index=False))

hs = rd("heat_step")
if hs is not None:
    print("\n=== heat_step: thermal response, time to 63% of final rise")
    if "failed" in hs:
        print("failed runs:", int(hs["failed"].fillna(0).sum()))
        hs = hs[hs["failed"].fillna(0) == 0]
    for (case, tun, geom), g in hs.groupby(["case", "tunneling", "geom"]):
        g = g.sort_values("t")
        dT = g.Tmax.values - 300
        fin = dT[-1]
        t63 = float(np.interp(0.632 * fin, dT, g.t.values)) if fin > 0 else float("nan")
        print(f"{case:12s} {tun:8s} {geom:6s} final dT={fin:.2f} K  t63={t63:.3g} s")

md = rd("mh_dir_map")
if md is not None:
    print("\n=== mh_dir_map: direction-dependent hole mass (top: m_c=0.17 tunnelling; edge: m_in), V_G=-3 V")
    if "failed" in md:
        print("failed runs:", int(md["failed"].fillna(0).sum()), "of", len(md))
        md = md[md["failed"].fillna(0) == 0]
    if len(md):
        pv = md.pivot_table(index=["case", "m_in"], columns="geom", values="R_tot")
        if "top50" in pv and "edge" in pv:
            pv["top/edge"] = pv["top50"] / pv["edge"]
        print(pv.to_string(float_format=lambda x: f"{x:.4g}"))

ud = rd("uq_dir")
if ud is not None:
    print("\n=== uq_dir: gated Monte Carlo with direction-dependent hole mass")
    if "failed" in ud:
        print("failed runs:", int(ud["failed"].fillna(0).sum()), "of", len(ud))
        ud = ud[ud["failed"].fillna(0) == 0]
    sh = ud[ud.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
    r = sh["edge"] / sh["top50"]
    print(f"shared S: n={len(r)} edge<top in {100 * (r < 1).mean():.1f}%  median edge/top={r.median():.3g}  "
          f"5-95%: {r.quantile(0.05):.3g}-{r.quantile(0.95):.3g}")
    ie = ud[ud.variant == "indep_edge"].set_index("sample").R_tot
    j = ie.index.intersection(sh.index)
    r2 = ie[j] / sh.loc[j, "top50"]
    print(f"indep S_edge: n={len(r2)} edge<top in {100 * (r2 < 1).mean():.1f}%")
    ic = ud[ud.variant == "indep_cnl"].set_index("sample").R_tot
    j3 = ic.index.intersection(sh.index)
    r3 = ic[j3] / sh.loc[j3, "top50"]
    print(f"indep pinning level (edge): n={len(r3)} edge<top in {100 * (r3 < 1).mean():.1f}%")
    mi = ud[ud.variant == "shared"].groupby("sample").m_in.first().reindex(r.index)
    for lo, hi in ((0.3, 0.5), (0.5, 0.8), (0.8, 1.5)):
        sel = (mi >= lo) & (mi < hi)
        if sel.any():
            print(f"  m_in {lo}-{hi}: n={int(sel.sum())} edge<top in {100 * (r[sel] < 1).mean():.1f}%  "
                  f"median edge/top={r[sel].median():.3g}")

td = rd("thickness_dir")
if td is not None:
    print("\n=== thickness_dir: Pt, V_G=-3 V, direction-dependent mass (m_par=0.9)")
    if "failed" in td:
        print("failed runs:", int(td["failed"].fillna(0).sum()), "of", len(td))
        td = td[td["failed"].fillna(0) == 0]
    pv = td.pivot_table(index=["case", "layers"], columns="geom", values="R_tot")
    if "top50" in pv and "edge" in pv:
        pv["top/edge"] = pv["top50"] / pv["edge"]
    print(pv.to_string(float_format=lambda x: f"{x:.4g}"))


def _ratio_table(df, idx, title):
    print(title)
    if df is None:
        return
    if "failed" in df:
        print("failed runs:", int(df["failed"].fillna(0).sum()), "of", len(df))
        df = df[df["failed"].fillna(0) == 0]
    pv = df.pivot_table(index=idx, columns="geom", values="R_tot")
    if "top50" in pv and "edge" in pv:
        pv["top/edge"] = pv["top50"] / pv["edge"]
    print(pv.to_string(float_format=lambda x: f"{x:.4g}"))


_ratio_table(rd("mass_controls"), ["case", "variant"], "\n=== mass_controls (Kane, m=0.9): iso vs swapped anisotropy")
et = rd("edge_T")
if et is not None:
    print("\n=== edge_T: edge R_tot vs interface transmission (Kane, tensor 0.9)")
    if "failed" in et:
        print("failed runs:", int(et["failed"].fillna(0).sum()), "of", len(et)); et = et[et["failed"].fillna(0) == 0]
    print(et.pivot_table(index="T_edge", columns="case", values="R_tot").to_string(float_format=lambda x: f"{x:.4g}"))
_ratio_table(rd("kane_map"), ["case", "variant", "m_par"], "\n=== kane_map: mass tensor, parabolic vs Kane vs Kane+oblique")
_ratio_table(rd("thickness_kane"), ["case", "layers"], "\n=== thickness_kane: Kane, tensor 0.9")
uk = rd("uq_kane")
if uk is not None:
    print("\n=== uq_kane: Monte Carlo with mass tensor and Kane dispersion")
    if "failed" in uk:
        print("failed runs:", int(uk["failed"].fillna(0).sum()), "of", len(uk)); uk = uk[uk["failed"].fillna(0) == 0]
    sh = uk[uk.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
    r = sh["edge"] / sh["top50"]
    print(f"shared: n={len(r)} edge<top {100 * (r < 1).mean():.1f}%  median edge/top={r.median():.3g}  5-95%: {r.quantile(0.05):.3g}-{r.quantile(0.95):.3g}")
    for var in ("indep_edge", "indep_cnl"):
        e = uk[uk.variant == var].set_index("sample").R_tot
        j = e.index.intersection(sh.index)
        print(f"{var}: n={len(j)} edge<top {100 * (e[j] / sh.loc[j, 'top50'] < 1).mean():.1f}%")
    mi = uk[uk.variant == "shared"].groupby("sample").m_in.first().reindex(r.index)
    for lo, hi in ((0.3, 0.5), (0.5, 0.8), (0.8, 1.5)):
        sel = (mi >= lo) & (mi < hi)
        print(f"  m_par {lo}-{hi}: n={int(sel.sum())} edge<top {100 * (r[sel] < 1).mean():.1f}%  median edge/top={r[sel].median():.3g}")

_mp = rd("mass_phase")
if _mp is not None:
    if "failed" in _mp:
        print("\n=== mass_phase failed runs:", int(_mp["failed"].fillna(0).sum()), "of", len(_mp)); _mp = _mp[_mp["failed"].fillna(0) == 0]
    pv = _mp.pivot_table(index=["case", "m_c", "m_par"], columns="geom", values="R_tot")
    pv["top/edge"] = pv["top50"] / pv["edge"]
    print("=== mass_phase: top/edge ratio (two-band, best straight path)")
    print(pv["top/edge"].unstack("m_par").to_string(float_format=lambda x: f"{x:.3g}"))
