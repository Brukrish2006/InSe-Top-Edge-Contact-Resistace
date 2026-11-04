"""Make all v2 figures from ../data_v2 into ../figures_v2 (PNG 300 dpi + PDF)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(__file__))
import inse2d as m  # noqa

D = os.path.join(os.path.dirname(__file__), "..", "data_v2")
F = os.path.join(os.path.dirname(__file__), "..", "figures_v2")
os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e4e0", "grid.linewidth": 0.6,
                     "savefig.dpi": 300, "axes.axisbelow": True, "font.family": "DejaVu Sans"})
# validated categorical palette (fixed order) + markers/linestyles as secondary encoding
STY = {"edge": dict(color="#2a78d6", marker="o", ls="-", label="Edge"),
       "top50": dict(color="#eb6834", marker="s", ls="--", label="Top, $L_c$=50 nm"),
       "top20": dict(color="#1baf7a", marker="^", ls="-.", label="Top, $L_c$=20 nm"),
       "top10": dict(color="#eda100", marker="D", ls=":", label="Top, $L_c$=10 nm")}
C3 = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#3a3a38"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(F, name + ".png"))
    fig.savefig(os.path.join(F, name + ".pdf"))
    plt.close(fig)
    print("saved", name)


def rd(n):
    return pd.read_csv(os.path.join(D, n + ".csv"))


# ---- Fig 1 verification
mg = rd("verification_mott_gurney"); sd = rd("verification_schottky_diode")
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
ax[0].loglog(mg.V, mg.J_sim_A_cm2, "o", color=C3[0], ms=6, label="DEVSIM v2 (1 µm, $N_A$=10$^8$)")
ax[0].loglog(mg.V, mg.J_MottGurney_A_cm2, "-", color=INK, lw=1.2, label=r"Mott–Gurney $\frac{9}{8}\varepsilon\mu V^2/L^3$")
ax[0].set_xlabel("V (V)"); ax[0].set_ylabel("J (A cm$^{-2}$)"); ax[0].legend(frameon=False)
ax[0].set_title("(a) Poisson–carrier coupling", fontsize=9, loc="left")
Dd = m.derived(m.DEFAULTS)
J0 = Dd["A_p"] * 300 ** 2 * np.exp(-0.40 / Dd["Vt"])
sd = sd[sd.V_drain.abs() > 1e-9]
V = sd.V_drain.values
ax[1].semilogy(V, np.abs(sd.J_A_cm2), "o", color=C3[1], ms=5, label="DEVSIM v2 Robin contact")
Vf = np.linspace(V.min(), V.max(), 200)
# identify forward polarity
fwdsign = 1 if abs(sd.J_A_cm2.values[-1]) > abs(sd.J_A_cm2.values[0]) else -1
ax[1].semilogy(Vf, J0 * np.abs(np.exp(fwdsign * Vf / Dd["Vt"]) - 1), "-", color=INK, lw=1.2,
               label=r"$A^*T^2e^{-\phi_B/kT}(e^{qV/kT}-1)$")
ax[1].set_xlabel("$V_{\\rm drain}$ (V)"); ax[1].set_ylabel("|J| (A cm$^{-2}$)")
ax[1].set_title("(b) Thermionic emission, $\\phi_B$=0.40 eV", fontsize=9, loc="left")
ax[1].legend(frameon=False)
save(fig, "fig1_verification")

# ---- Fig 2 barrier and R vs WF
wf = rd("wf_sweep")
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
phis = np.linspace(4.2, 5.5, 300)
for S, c, ls in ((1.0, C3[0], "-"), (0.5, C3[1], "--"), (0.3, C3[2], "-.")):
    pb = [m.barrier_heights(dict(m.DEFAULTS, S_pin=S), p)["phiBp_raw"] for p in phis]
    ax[0].plot(phis, pb, color=c, ls=ls, label=f"S = {S}")
ax[0].axhline(0, color=INK, lw=0.8)
ax[0].axvspan(5.0, 5.5, color="#ecebe6", zorder=0)
ax[0].text(5.03, 0.55, "$\\phi_{Bp}<0$ for S=1:\nclipped to 0\n(degenerate)", fontsize=7, color=INK)
ax[0].set_xlabel("$\\Phi_M$ (eV)"); ax[0].set_ylabel("Hole barrier $\\phi_{Bp}$ (eV)")
ax[0].legend(frameon=False, loc="lower left"); ax[0].set_title("(a) Input band alignment", fontsize=9, loc="left")
for g in ("top10", "top20", "top50", "edge"):
    d = wf[wf.geom == g].sort_values("Phi_M")
    s = STY[g]
    ax[1].semilogy(d.Phi_M, d.R_tot, color=s["color"], marker=s["marker"], ls=s["ls"], ms=4, label=s["label"])
ax[1].axvspan(5.0, 5.2, color="#ecebe6", zorder=0)
ax[1].set_xlabel("$\\Phi_M$ (eV)"); ax[1].set_ylabel("$R_{\\rm tot}$ at 0.05 V ($\\Omega\\cdot$cm)")
ax[1].legend(frameon=False); ax[1].set_title("(b) Total resistance, S = 1", fontsize=9, loc="left")
save(fig, "fig2_barrier_and_resistance_vs_wf")

# ---- Fig 3 decomposition
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
order = ["top10", "top20", "top50", "edge"]
for k, Phi in enumerate((4.8, 5.0)):
    d = wf[np.isclose(wf.Phi_M, Phi)].set_index("geom").loc[order]
    Rif = d.R_if.values; Rc2 = 2 * d.R_c.values; acc = Rc2 - Rif; Rch = d.R_ch.values
    x = np.arange(len(order))
    ax[k].bar(x, Rif, color=C3[1], width=0.6, label="Interface (thermionic/tunnel)", edgecolor="white", lw=1)
    ax[k].bar(x, acc, bottom=Rif, color="#f6b89c", width=0.6, label="Under-contact access", edgecolor="white", lw=1)
    ax[k].bar(x, Rch, bottom=Rc2, color=C3[0], width=0.6, label="Channel", edgecolor="white", lw=1)
    for xi, f in zip(x, d.frac_contact.values):
        ax[k].text(xi, d.R_tot.values[xi] * 1.02, f"{100*f:.0f}%", ha="center", va="bottom", fontsize=7, color=INK)
    ax[k].set_xticks(x); ax[k].set_xticklabels(["Top 10", "Top 20", "Top 50", "Edge"])
    ax[k].set_ylabel("$R$ at 0.05 V ($\\Omega\\cdot$cm)")
    ax[k].set_title(f"({'ab'[k]}) $\\Phi_M$ = {Phi} eV  (label: contact share)", fontsize=9, loc="left")
ax[0].legend(frameon=False, fontsize=7)
save(fig, "fig3_resistance_decomposition")

# ---- Fig 4 electro-thermal I-V
iv = rd("iv_electrothermal")
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
for Phi, ls in ((4.8, ":"), (5.0, "-")):
    for g in ("top10", "top50", "edge"):
        d = iv[(iv.geom == g) & np.isclose(iv.Phi_M, Phi)]
        s = STY[g]
        ax[0].semilogy(d.V, d.I_A_per_cm, color=s["color"], ls=ls, marker=s["marker"], ms=3,
                       label=f"{s['label']}, {Phi} eV")
        if Phi == 5.0:
            ax[1].plot(d.V, d.R_diff, color=s["color"], ls=s["ls"], marker=s["marker"], ms=3, label=s["label"])
            ax[2].plot(d.P_el, d.Tmax - 300, color=s["color"], ls=s["ls"], marker=s["marker"], ms=3, label=s["label"])
ax[0].set_xlabel("$V_{DS}$ (V)"); ax[0].set_ylabel("I (A/cm)"); ax[0].set_ylim(1e-6, 1e4)
ax[0].legend(frameon=False, fontsize=5.2, ncol=2, loc="upper center", columnspacing=0.6, handlelength=1.6)
ax[0].set_title("(a) I–V", fontsize=9, loc="left")
ax[1].set_xlabel("$V_{DS}$ (V)"); ax[1].set_ylabel("d$V$/d$I$ ($\\Omega\\cdot$cm)"); ax[1].set_yscale("log")
ax[1].set_title("(b) Differential R, 5.0 eV", fontsize=9, loc="left"); ax[1].legend(frameon=False, fontsize=6)
ax[2].set_xlabel("Dissipated power $IV$ (W/cm)"); ax[2].set_ylabel("$T_{\\max}-300$ K (K)")
ax[2].set_title("(c) Peak rise vs power", fontsize=9, loc="left"); ax[2].legend(frameon=False, fontsize=6)
save(fig, "fig4_electrothermal_iv")

# ---- Fig 5 SCLC / channel length
li = rd("lch_iv_isothermal"); ls_ = rd("lch_sweep")
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
for j, L in enumerate((25, 50, 100, 200)):
    d = li[(li.geom == "edge") & (li.Lch_nm == L)].sort_values("V")
    mloc = np.gradient(np.log(d.I_A_per_cm), np.log(d.V))
    ax[0].semilogx(d.V, mloc, color=C3[j], marker="os^D"[j], ms=3, label=f"$L_{{ch}}$={L} nm")
ax[0].axhline(1, color=INK, lw=0.7); ax[0].axhline(2, color=INK, lw=0.7, ls="--")
ax[0].set_xlabel("$V_{DS}$ (V)"); ax[0].set_ylabel("$m$ = dln$I$/dln$V$")
ax[0].set_title("(a) Local exponent, edge, 5.0 eV", fontsize=9, loc="left"); ax[0].legend(frameon=False, fontsize=6)
P = np.load(os.path.join(D, "profiles.npz"))
for j, V in enumerate((0.05, 0.5, 1.0, 2.0)):
    k = f"edge_Phi5.0_L100_V{V}"
    ax[1].semilogy(P[k + "_x"] * 1e7, P[k + "_p"] / 1e16, color=C3[j], label=f"{V} V")
ax[1].axhline(1, color=INK, lw=0.7, ls="--")
ax[1].set_xlabel("x (nm)"); ax[1].set_ylabel("$p/N_A$ (centre line)")
ax[1].set_title("(b) Hole profile, edge, 100 nm", fontsize=9, loc="left"); ax[1].legend(frameon=False, fontsize=6)
for g, Phi in (("edge", 5.0), ("top50", 5.0), ("edge", 4.8), ("top50", 4.8)):
    d = ls_[(ls_.geom == g) & np.isclose(ls_.Phi_M, Phi)].sort_values("Lch_nm")
    s = STY[g]
    ax[2].loglog(d.Lch_nm, d.R_tot, color=s["color"], marker=s["marker"], ls="-" if Phi == 5.0 else ":",
                 ms=4, label=f"{s['label'].split(',')[0]}, {Phi} eV")
ax[2].set_xticks([25, 50, 100, 200]); ax[2].set_xticklabels(["25", "50", "100", "200"]); ax[2].minorticks_off()
ax[2].set_xlabel("$L_{ch}$ (nm)"); ax[2].set_ylabel("$R_{\\rm tot}$ at 0.05 V ($\\Omega\\cdot$cm)")
ax[2].set_title("(c) Channel-length scaling", fontsize=9, loc="left"); ax[2].legend(frameon=False, fontsize=6)
save(fig, "fig5_sclc_channel_length")

# ---- Fig 6 tunneling & vdW
tv = rd("tunneling_vdw")
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
t = tv[tv.study == "tunneling"]
for g in ("top50", "edge"):
    a = t[t.geom == g].pivot_table(index="Phi_M", columns="tunneling", values="R_tot")
    s = STY[g]
    ax[0].plot(a.index, 100 * (1 - a[1.0] / a[0.0]), color=s["color"], marker=s["marker"], ls=s["ls"], label=s["label"])
ax[0].set_xlabel("$\\Phi_M$ (eV)"); ax[0].set_ylabel("R reduction by TFE (%)")
ax[0].set_title("(a) Field-dependent $\\Gamma(F)$", fontsize=9, loc="left"); ax[0].legend(frameon=False, fontsize=6)
c = tv[tv.study == "const_gamma"]
for g in ("top50", "edge"):
    for Phi, lsx in ((4.8, "-"), (4.6, ":")):
        d = c[(c.geom == g) & np.isclose(c.Phi_M, Phi)]
        s = STY[g]
        ax[1].semilogx(d.Gamma_const, d.R_tot / d.R_tot.iloc[0], color=s["color"], marker=s["marker"], ls=lsx,
                       label=f"{s['label'].split(',')[0]}, {Phi} eV")
ax[1].set_xlabel("constant $\\Gamma$"); ax[1].set_ylabel("$R(\\Gamma)/R(\\Gamma=1)$")
ax[1].set_title("(b) v1-style $\\Gamma$ sweep, corrected BC", fontsize=9, loc="left"); ax[1].legend(frameon=False, fontsize=6)
v = tv[tv.study == "vdw"]
edge = wf[wf.geom == "edge"].set_index("Phi_M").R_tot
for Phi, col in ((4.8, C3[1]), (5.0, C3[0])):
    d = v[np.isclose(v.Phi_M, Phi)].sort_values("d_vdW_nm")
    ax[2].semilogy(d.d_vdW_nm, d.R_tot, color=col, marker="s", label=f"Top 50, {Phi} eV")
    ax[2].axhline(edge.loc[Phi], color=col, lw=0.8, ls="--")
ax[2].text(0.02, edge.loc[4.8] * 0.6, "edge (dashed)", fontsize=7, color=INK)
ax[2].set_xlabel("vdW gap $d_{\\rm vdW}$ (nm)"); ax[2].set_ylabel("$R_{\\rm tot}$ ($\\Omega\\cdot$cm)")
ax[2].set_title("(c) Top-contact vdW gap", fontsize=9, loc="left"); ax[2].legend(frameon=False, fontsize=6)
save(fig, "fig6_tunneling_vdw")

# ---- Fig 7 thermal maps + BC sensitivity
fig = plt.figure(figsize=(7.2, 3.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1])
for r, g in enumerate(("top", "edge")):
    k = f"{g}_Phi5.0_L50_V2.0"
    axm = fig.add_subplot(gs[r, 0])
    x, y, T = P[k + "_2d_x"] * 1e7, P[k + "_2d_y"] * 1e7, P[k + "_2d_T"]
    tc = axm.tricontourf(x, y, T - 300, levels=20, cmap="Oranges")
    axm.invert_yaxis(); axm.set_aspect("auto"); axm.grid(False)
    axm.set_xlabel("x (nm)"); axm.set_ylabel("y (nm)")
    axm.set_title(f"({'ab'[r]}) {g} contacts, 5.0 eV, 2 V: $\\Delta T_{{\\max}}$={T.max()-300:.1f} K", fontsize=8.5, loc="left")
    fig.colorbar(tc, ax=axm, label="$\\Delta T$ (K)")
tb = rd("thermal_bc_sensitivity")
axb = fig.add_subplot(gs[:, 1])
cases = list(dict.fromkeys(tb.case + tb.drop(columns=["case"]).apply(
    lambda r: "", axis=1)))
labels = []
rowsT, rowsE = [], []
for _, rr in tb[tb.geom == "top50"].iterrows():
    lab = rr.case if rr.case == "baseline" else f"{rr.case}={[rr[c] for c in ('kappa','h_sub','h_contact') if c in rr and not pd.isna(rr[c])][0]:g}"
    labels.append(lab); rowsT.append(rr.dT)
for _, rr in tb[tb.geom == "edge"].iterrows():
    rowsE.append(rr.dT)
yy = np.arange(len(labels))
axb.barh(yy - 0.2, rowsT, height=0.38, color=STY["top50"]["color"], label="Top 50")
axb.barh(yy + 0.2, rowsE, height=0.38, color=STY["edge"]["color"], label="Edge")
axb.set_yticks(yy); axb.set_yticklabels(labels, fontsize=7); axb.set_xscale("log")
axb.set_xlabel("$\\Delta T_{\\max}$ at 2 V (K)"); axb.legend(frameon=False)
axb.set_title("(c) Thermal model sensitivity", fontsize=9, loc="left")
save(fig, "fig7_thermal")

# ---- Fig 8 sensitivity (material + doping)
ms = rd("material_sensitivity"); dp = rd("doping_sweep")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for k, Phi in enumerate((4.8, 5.0)):
    lab, lo, hi = [], [], []
    for g in ("top50", "edge"):
        base = wf[(wf.geom == g) & np.isclose(wf.Phi_M, Phi)].R_tot.values[0]
        for par in ("mu_p", "E_V", "E_g", "eps_r", "mh_dos"):
            d = ms[(ms.param == par) & (ms.geom == g) & np.isclose(ms.Phi_M, Phi)]
            r = d.R_tot / base
            nm_ = {"E_g": "E_g (χ fixed)", "E_V": "E_V (via χ)", "mh_dos": "m_h"}.get(par, par)
            lab.append(f"{g}: {nm_} {d.value.min():g}–{d.value.max():g}")
            lo.append(r.min()); hi.append(r.max())
        d = dp[(dp.geom == g) & np.isclose(dp.Phi_M, Phi)]
        r = d.R_tot / base
        lab.append(f"{g}: N_A 1e15–1e17"); lo.append(r.min()); hi.append(r.max())
    yy = np.arange(len(lab))
    for i in yy:
        col = STY["top50"]["color"] if lab[i].startswith("top") else STY["edge"]["color"]
        ax[k].plot([lo[i], hi[i]], [i, i], color=col, lw=5, solid_capstyle="butt")
    ax[k].axvline(1, color=INK, lw=0.8)
    ax[k].set_xscale("log"); ax[k].set_yticks(yy); ax[k].set_yticklabels(lab, fontsize=6)
    ax[k].set_xlabel("$R/R_{\\rm baseline}$"); ax[k].set_title(f"({'ab'[k]}) $\\Phi_M$={Phi} eV", fontsize=9, loc="left")
save(fig, "fig8_parameter_sensitivity")

# ---- Fig S1 mesh
me = rd("mesh_convergence")
fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.6))
for (g, Phi), d in me.groupby(["geom", "Phi_M"]):
    s = STY[g]
    d = d.sort_values("nodes")
    ax.semilogx(d.nodes, d.R_tot / d.R_tot.iloc[-1], color=s["color"], marker=s["marker"],
                ls="-" if Phi == 5.0 else ":", label=f"{s['label'].split(',')[0]}, {Phi} eV")
ax.set_xlabel("nodes"); ax.set_ylabel("$R/R_{\\rm finest}$"); ax.legend(frameon=False, fontsize=6)
save(fig, "figS1_mesh_convergence")
