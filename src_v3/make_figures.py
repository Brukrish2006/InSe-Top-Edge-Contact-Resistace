"""v3 figures and derived analyses (compact model, TLM fits, UQ statistics).
Reads ../data_v3, writes ../figures_v3 and ../data_v3/derived_*.csv"""
import math, os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
import tcad2d as m  # noqa

D = os.path.join(os.path.dirname(__file__), "..", "data_v3")
F = os.path.join(os.path.dirname(__file__), "..", "figures_v3")
os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "lines.linewidth": 1.6, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e4e0",
                     "grid.linewidth": 0.6, "savefig.dpi": 300, "axes.axisbelow": True})
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
STY = {"edge": dict(color=C[0], marker="o", ls="-", label="Edge"),
       "top50": dict(color=C[1], marker="s", ls="--", label="Top, $L_c$=50 nm"),
       "top10": dict(color=C[3], marker="D", ls=":", label="Top, $L_c$=10 nm")}
INK = "#3a3a38"
DER = {}


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(F, name + ".png")); fig.savefig(os.path.join(F, name + ".pdf"))
    plt.close(fig); print("saved", name)


def rd(n):
    p = os.path.join(D, n + ".csv")
    return pd.read_csv(p) if os.path.exists(p) else None


# ---------------------------------------------------------------- Fig 1 verification
mg, s3, s35, tr = rd("verification_mott_gurney"), rd("verification_schottky_300K"), rd("verification_schottky_350K"), rd("verification_transient")
if mg is not None:
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for k, st in enumerate(("boltzmann", "fd")):
        d = mg[mg.statistics == st]
        ax[0].loglog(d.V, d.J_sim, "os"[k], color=C[k], ms=5, mfc="none" if k else C[k], label=f"DEVSIM ({'Boltzmann' if k == 0 else 'Fermi–Dirac'})")
    d = mg[mg.statistics == "fd"]
    ax[0].loglog(d.V, d.J_MG, "-", color=INK, lw=1, label="Mott–Gurney")
    ax[0].set_xlabel("V (V)"); ax[0].set_ylabel("J (A cm$^{-2}$)"); ax[0].legend(frameon=False)
    ax[0].set_title("(a) SCLC, 1 µm film", fontsize=9, loc="left")
    Dd = m.derived(dict(m.DEFAULTS, mh_dos=0.17), 1.30)
    for k, (dd, T0) in enumerate(((s3, 300.0), (s35, 350.0))):
        if dd is None:
            continue
        Vt = m.kq * T0; A = 120.17 * 0.17
        J0 = A * T0 ** 2 * math.exp(-0.40 / Vt)
        dd = dd[dd.V.abs() > 1e-9]
        ax[1].semilogy(dd.V, dd.J.abs(), "o", color=C[k], ms=4, label=f"DEVSIM {T0:.0f} K")
        vv = np.linspace(dd.V.min(), dd.V.max(), 300)
        sg = 1 if abs(dd.J.values[-1]) > abs(dd.J.values[0]) else -1
        ax[1].semilogy(vv, J0 * np.abs(np.exp(sg * vv / Vt) - 1), "-", color=INK, lw=0.9)
    ax[1].set_xlabel("$V$ (V)"); ax[1].set_ylabel("|J| (A cm$^{-2}$)"); ax[1].legend(frameon=False)
    ax[1].set_title("(b) Thermionic emission, 0.40 eV", fontsize=9, loc="left")
    if tr is not None:
        ax[2].semilogx(tr.t_s, tr.Tmax_K - 300, color=C[0], label="transient")
        ax[2].axhline(tr.Tmax_K.values[-1] - 300, color=INK, lw=0.8, ls="--", label="DC")
        ax[2].set_xlabel("t (s)"); ax[2].set_ylabel("$\\Delta T_{\\max}$ (K)"); ax[2].legend(frameon=False)
        ax[2].set_title("(c) Transient → DC (edge, 1 V)", fontsize=9, loc="left")
    save(fig, "fig1_verification")

# ---------------------------------------------------------------- Fig 2 thickness / band
th = rd("thickness")
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
N = np.linspace(2, 20, 200)
Eg = [m.band(m.INSE, n * 0.83e-7)[0] for n in N]
chi = [m.band(m.INSE, n * 0.83e-7)[1] for n in N]
ax[0].plot(N, Eg, color=C[0], label="$E_g(N)$ model")
ax[0].plot(N, np.array(chi) + np.array(Eg), color=C[1], ls="--", label="$E_V(N)$ below vacuum")
ax[0].plot([3, 5, 2], [1.69, 1.44, 1.95], "k^", ms=5, label="PL (lit.)")
ax[0].axhline(1.26, color=INK, lw=0.6, ls=":")
ax[0].set_xlabel("layers N"); ax[0].set_ylabel("energy (eV)"); ax[0].legend(frameon=False)
ax[0].set_title("(a) Layer-dependent band edges", fontsize=9, loc="left")
if th is not None:
    for g in ("top50", "edge"):
        for Phi, ls in ((4.8, ":"), (5.0, "-")):
            d = th[(th.geom == g) & np.isclose(th.Phi_M, Phi) & (th.study == "ungated")].sort_values("layers")
            ax[1].semilogy(d.layers, d.R_tot, color=STY[g]["color"], marker=STY[g]["marker"], ls=ls,
                           label=f"{STY[g]['label'].split(',')[0]}, {Phi} eV")
        d = th[(th.geom == g) & (th.study == "gated_on")].sort_values("layers")
        ax[1].semilogy(d.layers, d.R_tot, color=STY[g]["color"], marker="x", ls="-.", label=f"{STY[g]['label'].split(',')[0]}, gated on")
    ax[1].set_xlabel("layers N"); ax[1].set_ylabel("$R_{\\rm tot}$ at 0.05 V ($\\Omega\\cdot$cm)")
    ax[1].legend(frameon=False, fontsize=6); ax[1].set_title("(b) Thickness dependence", fontsize=9, loc="left")
save(fig, "fig2_thickness")

# ---------------------------------------------------------------- Fig 3 WF sweep, FD, tunnelling
wf = rd("wf_sweep")
if wf is not None:
    base = wf[wf.get("stat", pd.Series(index=wf.index)).isna() & wf.get("tunneling", pd.Series(index=wf.index)).isna()]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    for g in ("top10", "top50", "edge"):
        d = base[base.geom == g].sort_values("Phi_M"); s = STY[g]
        ax[0].semilogy(d.Phi_M, d.R_tot, color=s["color"], marker=s["marker"], ls=s["ls"], ms=4, label=s["label"] + " (FD)")
        b = wf[(wf.geom == g) & (wf.get("stat") == "boltzmann")].sort_values("Phi_M")
        ax[0].semilogy(b.Phi_M, b.R_tot, color=s["color"], marker="x", ls="none", ms=6)
    ax[0].axvline(5.0, color=INK, lw=0.6, ls=":")
    ax[0].set_xlabel("$\\Phi_M$ (eV)"); ax[0].set_ylabel("$R_{\\rm tot}$ ($\\Omega\\cdot$cm)")
    ax[0].legend(frameon=False, fontsize=6); ax[0].set_title("(a) Ungated; × = Boltzmann (clipped)", fontsize=9, loc="left")
    for k, g in enumerate(("top50", "edge")):
        for tun, ls in (("off", ":"), ("local", "--")):
            d = wf[(wf.geom == g) & (wf.get("tunneling") == tun)].sort_values("Phi_M")
            b0 = base[(base.geom == g) & base.Phi_M.isin(d.Phi_M)].sort_values("Phi_M")
            ax[1].plot(d.Phi_M.values, b0.R_tot.values / d.R_tot.values, color=STY[g]["color"], ls=ls, marker=STY[g]["marker"],
                       label=f"{STY[g]['label'].split(',')[0]}: nonlocal / {tun}")
    ax[1].axhline(1, color=INK, lw=0.6)
    ax[1].set_xlabel("$\\Phi_M$ (eV)"); ax[1].set_ylabel("$R_{\\rm nonlocal}/R_{\\rm ref}$")
    ax[1].legend(frameon=False, fontsize=6); ax[1].set_title("(b) Tunnelling model, ungated", fontsize=9, loc="left")
    save(fig, "fig3_wf_statistics_tunnelling")

# ---------------------------------------------------------------- Fig 4 anisotropy
an = rd("anisotropy")
if an is not None:
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    a1 = an[an.get("study").isna()]
    for g in ("top10", "top50", "edge"):
        for Phi, ls in ((4.8, "-"), (5.0, ":")):
            d = a1[(a1.geom == g) & np.isclose(a1.Phi_M, Phi)].sort_values("rmu")
            ax[0].loglog(d.rmu, d.R_tot, color=STY[g]["color"], marker=STY[g]["marker"], ls=ls,
                         label=f"{STY[g]['label'].split(',')[0]}{' '+STY[g]['label'].split(',')[1] if ',' in STY[g]['label'] else ''}, {Phi}")
    ax[0].set_xlabel("$\\mu_\\perp/\\mu_\\parallel$"); ax[0].set_ylabel("$R_{\\rm tot}$ ($\\Omega\\cdot$cm)")
    ax[0].legend(frameon=False, fontsize=5.5); ax[0].set_title("(a) Mobility anisotropy", fontsize=9, loc="left")
    a2 = an[an.get("study") == "Lc"]
    for rmu, col in ((1.0, C[0]), (0.01, C[1])):
        for Phi, ls in ((4.8, "-"), (5.0, ":")):
            d = a2[np.isclose(a2.rmu, rmu) & np.isclose(a2.Phi_M, Phi)].sort_values("Lc_nm")
            ax[1].loglog(d.Lc_nm, 2 * d.R_c, color=col, ls=ls, marker="s", label=f"$\\mu_\\perp/\\mu_\\parallel$={rmu:g}, {Phi} eV")
    ax[1].set_xlabel("$L_c$ (nm)"); ax[1].set_ylabel("$2R_c$ ($\\Omega\\cdot$cm)")
    ax[1].legend(frameon=False, fontsize=5.5); ax[1].set_title("(b) Top-contact length", fontsize=9, loc="left")
    a3 = an[an.get("study") == "thermal"]
    gl = ["top10", "top50", "edge"]; x = np.arange(3)
    for k, kp in enumerate((0.085, 0.0076)):
        d = a3[np.isclose(a3.kperp, kp)].set_index("geom").reindex(gl)
        ax[2].bar(x + (k - 0.5) * 0.38, d.Tmax - 300, width=0.36, color=C[k], label=f"$\\kappa_\\perp$={kp} W/cmK")
    ax[2].set_xticks(x); ax[2].set_xticklabels(["Top 10", "Top 50", "Edge"]); ax[2].set_ylabel("$\\Delta T_{\\max}$ at 1 V (K)")
    ax[2].legend(frameon=False, fontsize=6); ax[2].set_title("(c) Thermal anisotropy", fontsize=9, loc="left")
    save(fig, "fig4_anisotropy")

# ---------------------------------------------------------------- Fig 5 pinning top vs edge
pn = rd("pinning_top_edge")
if pn is not None:
    fig, ax = plt.subplots(1, 1, figsize=(3.6, 2.8))
    for k, Se in enumerate((1.0, 0.7, 0.5, 0.3)):
        d = pn[(pn.geom == "edge") & np.isclose(pn.get("S_edge"), Se)].sort_values("Phi_M")
        ax.semilogy(d.Phi_M, d.R_tot, color=C[0], alpha=1 - 0.2 * k, marker="o", ls="-", label=f"edge, $S_{{edge}}$={Se}")
    for k, St in enumerate((1.0, 0.5)):
        d = pn[(pn.geom == "top50") & np.isclose(pn.get("S_top"), St)].sort_values("Phi_M")
        ax.semilogy(d.Phi_M, d.R_tot, color=C[1], alpha=1 - 0.4 * k, marker="s", ls="--", label=f"top 50, $S_{{top}}$={St}")
    ax.set_xlabel("$\\Phi_M$ (eV)"); ax.set_ylabel("$R_{\\rm tot}$ ($\\Omega\\cdot$cm)"); ax.legend(frameon=False, fontsize=6)
    save(fig, "fig5_pinning_top_vs_edge")

# ---------------------------------------------------------------- Fig 6 gated transfer
ga = rd("gated")
if ga is not None:
    tf = ga[(ga.study == "transfer") & ga.get("tunneling", pd.Series(index=ga.index)).isna()]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
    for g in ("top50", "edge"):
        for Phi, ls in ((4.8, ":"), (5.0, "-"), (5.2, "--")):
            d = tf[(tf.geom == g) & np.isclose(tf.Phi_M, Phi)].sort_values("VG")
            ax[0].semilogy(d.VG, d.I.abs(), color=STY[g]["color"], ls=ls, marker=STY[g]["marker"], ms=2.5,
                           label=f"{STY[g]['label'].split(',')[0]}, {Phi} eV")
            if Phi == 5.0:
                d = d[d.VG <= 0.0]
                ax[1].plot(d.VG, 100 * d.tun_share, color=STY[g]["color"], ls=ls, marker=STY[g]["marker"], ms=2.5, label=f"{STY[g]['label'].split(',')[0]}, {Phi} eV")
                ax[2].semilogy(d.VG, 2 * d.R_c, color=STY[g]["color"], ls=ls, marker=STY[g]["marker"], ms=2.5, label=f"{STY[g]['label'].split(',')[0]}, 2$R_c$")
                ax[2].semilogy(d.VG, d.R_ch, color=STY[g]["color"], ls=":", marker="x", ms=2.5, label=f"{STY[g]['label'].split(',')[0]}, $R_{{ch}}$")
    for tun, mk in (("off", "+"), ("local", "1")):
        for g in ("top50", "edge"):
            d = ga[(ga.study == "transfer") & (ga.get("tunneling") == tun) & (ga.geom == g)].sort_values("VG")
            ax[0].semilogy(d.VG, d.I.abs(), color=STY[g]["color"], ls="none", marker=mk, ms=5)
    ax[0].set_xlabel("$V_G$ (V)"); ax[0].set_ylabel("$I_D$ (A/cm), $V_D$=0.05 V"); ax[0].legend(frameon=False, fontsize=5.5)
    ax[0].set_title("(a) Transfer (+ off, ⅄ local: 4.8 eV)", fontsize=8.5, loc="left")
    ax[1].set_xlabel("$V_G$ (V)"); ax[1].set_ylabel("tunnelling share (%)"); ax[1].legend(frameon=False, fontsize=6)
    ax[1].set_title("(b) Injection by tunnelling", fontsize=9, loc="left")
    ax[2].set_xlabel("$V_G$ (V)"); ax[2].set_ylabel("($\\Omega\\cdot$cm)"); ax[2].legend(frameon=False, fontsize=5.5)
    ax[2].set_title("(c) Contact vs channel, 5.0 eV", fontsize=9, loc="left")
    save(fig, "fig6_gated_transfer")

# ---------------------------------------------------------------- Fig 7 thermal (output, Peltier, transient)
tg, tt = rd("thermal_gated"), rd("thermal_transient")
if ga is not None and tg is not None:
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
    oc = ga[ga.study == "output"]
    for g in ("top50", "edge"):
        for Phi, ls in ((4.8, ":"), (5.0, "-")):
            d = oc[(oc.geom == g) & np.isclose(oc.Phi_M, Phi)].sort_values("V")
            ax[0].plot(d.P_el, d.Tmax - 300, color=STY[g]["color"], ls=ls, marker=STY[g]["marker"], ms=3, label=f"{STY[g]['label'].split(',')[0]}, {Phi} eV")
    ax[0].set_xlabel("$IV$ (W/cm)"); ax[0].set_ylabel("$\\Delta T_{\\max}$ (K)"); ax[0].legend(frameon=False, fontsize=6)
    ax[0].set_title("(a) Gated, $V_G$=−3 V", fontsize=9, loc="left")
    cases = ["full", "no_peltier", "bulk_joule_only", "isotropic_kappa", "contact_sink", "R_TBR_1e-9", "isothermal_stats"]
    y = np.arange(len(cases))
    for k, g in enumerate(("top50", "edge")):
        d = tg[(tg.geom == g) & np.isclose(tg.V, 2.0)].set_index("case").reindex(cases)
        ax[1].barh(y + (k - 0.5) * 0.38, d.Tmax - 300, height=0.36, color=STY[g]["color"], label=STY[g]["label"].split(",")[0])
    ax[1].set_yticks(y); ax[1].set_yticklabels(cases, fontsize=6); ax[1].set_xlabel("$\\Delta T_{\\max}$ at 2 V (K)")
    ax[1].legend(frameon=False, fontsize=6); ax[1].set_title("(b) Heat-model terms", fontsize=9, loc="left")
    if tt is not None:
        for g in ("top50", "edge"):
            d = tt[tt.geom == g]
            ax[2].semilogx(d.t, d.Tmax - 300, color=STY[g]["color"], ls=STY[g]["ls"], label=STY[g]["label"].split(",")[0])
            fin = d.Tmax.values[-1] - 300
            t63 = float(np.interp(0.632 * fin, d.Tmax.values - 300, d.t.values))
            DER[f"tau63_{g}_s"] = t63
        ax[2].set_xlabel("t (s)"); ax[2].set_ylabel("$\\Delta T_{\\max}$ (K)"); ax[2].legend(frameon=False, fontsize=6)
        ax[2].set_title("(c) 0→2 V step", fontsize=9, loc="left")
    save(fig, "fig7_thermal")

# ---------------------------------------------------------------- Fig 8 MoS2 calibration
ms = rd("mos2_calibration")
if ms is not None and (ms.step == "Lc_sweep").any():
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    out = {}
    for k, T0 in enumerate((300.0, 400.0)):
        d = ms[(ms.step == "Lc_sweep") & np.isclose(ms.T0, T0)].sort_values("Lc_nm")
        if d.empty:
            continue
        Rsh = float(d.Rsh.median())
        f = lambda Lc, rc, LT: np.sqrt(rc * Rsh) / np.tanh(Lc / LT)  # noqa
        try:
            (rc, LT), _ = curve_fit(f, d.Lc_nm.values * 1e-7, d.R_c.values, p0=(3e-7, 35e-7), maxfev=20000)
        except Exception:
            rc, LT = np.nan, np.nan
        out[T0] = dict(rho_c=rc, L_T_nm=LT * 1e7, Rsh=Rsh, Rc_long_ohm_um=float(d.R_c_ohm_um.values[-1]))
        ax[0].plot(d.Lc_nm, d.R_c_ohm_um, "o", color=C[k], label=f"TCAD {T0:.0f} K")
        ll = np.linspace(10, 210, 200)
        if np.isfinite(rc):
            ax[0].plot(ll, f(ll * 1e-7, rc, LT) * 1e4, "-", color=C[k], lw=1, label=f"TLM fit: $L_T$={LT*1e7:.0f} nm")
    ax[0].axhline(740, color=INK, lw=0.8, ls="--"); ax[0].text(150, 790, "exp. 740 Ω·µm (300 K)", fontsize=6.5)
    ax[0].set_xlabel("$L_c$ (nm)"); ax[0].set_ylabel("$R_c$ (Ω·µm)"); ax[0].legend(frameon=False, fontsize=6)
    ax[0].set_title("(a) Au/MoS$_2$, $n$=10$^{13}$ cm$^{-2}$", fontsize=9, loc="left")
    lab_ = ["$\\rho_c$ 300 K\n(×10$^{-7}$ Ω cm²)", "$L_T$ 300 K\n(nm)", "$\\rho_c$ 400 K\n(×10$^{-7}$)", "$L_T$ 400 K\n(nm)"]
    expv = [3.0, 35.0, 1.5, 20.0]; expe = [1.0, 10.0, 0.5, 7.0]
    simv = [out.get(300.0, {}).get("rho_c", np.nan) * 1e7, out.get(300.0, {}).get("L_T_nm", np.nan),
            out.get(400.0, {}).get("rho_c", np.nan) * 1e7, out.get(400.0, {}).get("L_T_nm", np.nan)]
    xx = np.arange(4)
    ax[1].bar(xx - 0.2, expv, width=0.38, yerr=expe, color="#c9c8c2", label="experiment", capsize=2)
    ax[1].bar(xx + 0.2, simv, width=0.38, color=C[0], label="TCAD (calibrated at 300 K $R_c$ only)")
    ax[1].set_xticks(xx); ax[1].set_xticklabels(lab_, fontsize=6.5); ax[1].legend(frameon=False, fontsize=6)
    ax[1].set_title("(b) Prediction vs measurement", fontsize=9, loc="left")
    save(fig, "fig8_mos2_calibration")
    cal = ms[ms.step == "calibration"]
    DER["mos2"] = dict(fits=out, T_vdW_cal=float(cal.T_vdW.values[0]) if len(cal) else np.nan,
                       mode=str(cal["mode"].values[0]) if "mode" in cal and len(cal) else "")

# ---------------------------------------------------------------- Fig 9 compact model
cm = rd("compact_model_data")
if cm is not None:
    Pm = m.DEFAULTS
    Eg, chi, _ = m.band(Pm, 5e-7)
    Dd = m.derived(Pm, Eg)
    q, Vt, t = m.q, Dd["Vt"], 5e-7
    rows = []
    for _, r in cm.iterrows():
        peq = Dd["Nv"] * m.F12(-r.phiBp / Vt)
        rho_if = Vt / (q * Dd["vR_p"] * r.T_vdW * peq)
        for mode, pbar in (("p_eq", peq), ("p_under", r.p_under)):
            rv = t / (2 * q * Pm["mu_p"] * r.rmu * pbar)
            rho = rho_if + rv
            Rsh = 1 / (q * Pm["mu_p"] * pbar * t)
            LT = math.sqrt(rho / Rsh)
            Rc = math.sqrt(rho * Rsh) / math.tanh(r.Lc_nm * 1e-7 / LT)
            rows.append(dict(mode=mode, Rc_tcad=r.R_c, Rc_cm=Rc, LT_nm=LT * 1e7, Lc_nm=r.Lc_nm, rmu=r.rmu,
                             d_vdW_nm=r.d_vdW_nm, phiBp=r.phiBp))
    cdf = pd.DataFrame(rows); cdf.to_csv(os.path.join(D, "derived_compact_model.csv"), index=False)
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9))
    for k, mode in enumerate(("p_eq", "p_under")):
        d = cdf[cdf["mode"] == mode]
        err = np.log10(d.Rc_cm / d.Rc_tcad)
        DER[f"compact_{mode}"] = dict(median_abs_pct=float(np.median(np.abs(10 ** err - 1)) * 100),
                                     p90_abs_pct=float(np.percentile(np.abs(10 ** err - 1), 90) * 100),
                                     max_abs_pct=float(np.max(np.abs(10 ** err - 1)) * 100), n=int(len(d)))
        for rmu, mk in ((1.0, "o"), (0.01, "s")):
            dd = d[np.isclose(d.rmu, rmu)]
            ax[k].loglog(dd.Rc_tcad, dd.Rc_cm, mk, color=C[k], ms=3.5, mfc="none" if rmu == 1 else C[k],
                         label=f"$\\mu_\\perp/\\mu_\\parallel$={rmu:g}")
        lim = [min(d.Rc_tcad.min(), d.Rc_cm.min()) * 0.7, max(d.Rc_tcad.max(), d.Rc_cm.max()) * 1.4]
        ax[k].plot(lim, lim, color=INK, lw=0.8)
        ax[k].set_xlabel("$R_c$ TCAD ($\\Omega\\cdot$cm)"); ax[k].set_ylabel("$R_c$ compact ($\\Omega\\cdot$cm)")
        ttl = "(a) $\\bar p=p_{\\rm eq}$ (predictive)" if mode == "p_eq" else "(b) $\\bar p$ from TCAD"
        ax[k].set_title(f"{ttl}; median err {DER[f'compact_{mode}']['median_abs_pct']:.0f}%", fontsize=8.5, loc="left")
        ax[k].legend(frameon=False, fontsize=6)
    save(fig, "fig9_compact_model")

# ---------------------------------------------------------------- Fig 10 UQ
uq = rd("uq_samples")
if uq is not None:
    uq = uq[uq.get("failed", pd.Series(0, index=uq.index)).fillna(0) == 0]
    pars = ["dEV", "S", "d_vdW", "Phi_vdW", "mu_p", "rmu", "NA", "eps_r", "mh"]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
    stats = {}
    for k, (g, Phi) in enumerate((("top50", 4.8), ("edge", 4.8), ("top50", 5.0), ("edge", 5.0))):
        d = uq[(uq.geom == g) & np.isclose(uq.Phi_M, Phi)]
        if len(d) < 10:
            continue
        lr = np.log10(d.R_tot)
        ax[0].hist(lr, bins=25, histtype="step", color=STY[g]["color"], ls="-" if Phi == 4.8 else "--",
                   label=f"{STY[g]['label'].split(',')[0]}, {Phi} eV")
        rho = {p: spearmanr(np.log10(d[p]) if p in ("mu_p", "rmu", "NA") else d[p], lr)[0] for p in pars}
        # first-order Sobol estimate by binning (correlation ratio)
        S1 = {}
        for p in pars:
            xb = pd.qcut(d[p].rank(method="first"), 8, labels=False)
            S1[p] = float(lr.groupby(xb).mean().var() / lr.var())
        stats[f"{g}_{Phi}"] = dict(median_R=float(10 ** lr.median()), p5=float(10 ** np.percentile(lr, 5)),
                                  p95=float(10 ** np.percentile(lr, 95)), spearman=rho, S1=S1, n=int(len(d)))
    ax[0].set_xlabel("log$_{10}$ $R_{\\rm tot}$ ($\\Omega\\cdot$cm)"); ax[0].set_ylabel("samples")
    ax[0].legend(frameon=False, fontsize=6); ax[0].set_title("(a) Monte Carlo output", fontsize=9, loc="left")
    x = np.arange(len(pars))
    for k, key in enumerate(("top50_4.8", "edge_4.8", "top50_5.0", "edge_5.0")):
        if key in stats:
            ax[1].bar(x + (k - 1.5) * 0.2, [stats[key]["S1"][p] for p in pars], width=0.19, color=C[k], label=key.replace("_", ", ") + " eV")
    ax[1].set_xticks(x); ax[1].set_xticklabels(["$\\Delta E_V$", "S", "$d_{vdW}$", "$\\Phi_{vdW}$", "$\\mu_p$", "$\\mu_\\perp/\\mu_\\parallel$", "$N_A$", "$\\varepsilon_r$", "$m_h$"], fontsize=6.5, rotation=30)
    ax[1].set_ylabel("first-order index $S_1$"); ax[1].legend(frameon=False, fontsize=6)
    ax[1].set_title("(b) Variance decomposition", fontsize=9, loc="left")
    save(fig, "fig10_uncertainty")
    DER["uq"] = stats

json.dump(DER, open(os.path.join(D, "derived_summary.json"), "w"), indent=1, default=float)
print(json.dumps(DER, indent=1, default=float)[:3000])
