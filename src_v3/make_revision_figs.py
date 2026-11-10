"""Figures for the revision studies (run_revision.py). Reads ../data_v3, writes ../figures_v3/fig13-15."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data_v3")
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures_v3")
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "lines.linewidth": 1.6, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e4e0",
                     "grid.linewidth": 0.6, "savefig.dpi": 300, "axes.axisbelow": True})
BLUE, ORANGE, GREEN, YELLOW, INK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#3a3a38", "#9a9a96"
rd = lambda n: pd.read_csv(os.path.join(D, n + ".csv"))
ok = lambda d: d[d["failed"].fillna(0) == 0] if "failed" in d else d


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(F, f"{name}.{ext}"))
    plt.close(fig); print("saved", name)


def first_order(x, y, nb=8):
    def idx(xx):
        q = pd.qcut(pd.Series(xx).rank(method="first"), nb, labels=False)
        return pd.Series(y).groupby(q.values).mean().var(ddof=0) / np.var(y)
    rng = np.random.default_rng(0)
    return idx(x), np.mean([idx(rng.permutation(x)) for _ in range(50)])


# ------------------------------------------------------------------ Fig 13: gated Monte Carlo
u = ok(rd("uq_gated")); u = u[np.isclose(u.VG, -3.0)]
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
for (geom, var), col, ls, lab in [(("edge", "shared"), BLUE, "-", "Edge"), (("top50", "shared"), ORANGE, "-", "Top, 50 nm"),
                                  (("edge", "indep_edge"), BLUE, ":", r"Edge, independent $S_{\rm edge}$")]:
    y = np.sort(np.log10(u[(u.geom == geom) & (u.variant == var)].R_tot.values))
    ax[0].plot(y, np.arange(1, len(y) + 1) / len(y), color=col, ls=ls, label=lab)
ax[0].set_xlabel(r"$\log_{10}R_{\rm tot}$ ($\Omega$cm)"); ax[0].set_ylabel("Cumulative fraction")
ax[0].legend(frameon=False, loc="lower right", fontsize=6.5); ax[0].set_title("(a)", loc="left", fontsize=9)
sh = u[u.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
ie = u[u.variant == "indep_edge"].set_index("sample").R_tot
ax[1].loglog(sh["top50"], sh["edge"], "o", color=BLUE, ms=2.8, alpha=0.7, label="shared $S$")
j = sh.index.intersection(ie.index)
ax[1].loglog(sh.loc[j, "top50"], ie.loc[j], "^", color=GREEN, ms=2.8, alpha=0.6, mfc="none", label=r"independent $S_{\rm edge}$")
lim = [min(sh.min().min(), ie.min()) / 3, sh["top50"].max() * 3]
ax[1].plot(lim, lim, color=INK, lw=0.8); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel(r"$R_{\rm tot}^{\rm top}$ ($\Omega$cm)"); ax[1].set_ylabel(r"$R_{\rm tot}^{\rm edge}$ ($\Omega$cm)")
ax[1].legend(frameon=False, loc="upper left", fontsize=6.5); ax[1].set_title("(b)", loc="left", fontsize=9)
par = ["cnl_off", "mh", "S", "eps_r", "dEV", "Phi_vdW", "d_vdW", "rmu", "mu_p", "NA"]
lab = {"cnl_off": "pin. level", "mh": "$m_h$", "S": "$S$", "eps_r": r"$\varepsilon_r$", "dEV": r"$\Delta E_V$",
       "Phi_vdW": r"$\Phi_{\rm vdW}$", "d_vdW": r"$d_{\rm vdW}$", "rmu": r"$\mu_\perp/\mu_\parallel$", "mu_p": r"$\mu_p$", "NA": "$N_A$"}
xs = np.arange(len(par)); floors = []
for k, (geom, col) in enumerate((("edge", BLUE), ("top50", ORANGE))):
    g = u[(u.geom == geom) & (u.variant == "shared")]; y = np.log10(g.R_tot.values)
    vals = []
    for p in par:
        si, fl = first_order(g[p].values, y); vals.append(si); floors.append(fl)
    ax[2].bar(xs + (k - 0.5) * 0.38, vals, 0.38, color=col, label="Edge" if geom == "edge" else "Top, 50 nm")
ax[2].axhline(np.mean(floors), color=GREY, ls="--", lw=0.9, label="noise floor")
ax[2].set_xticks(xs); ax[2].set_xticklabels([lab[p] for p in par], rotation=60, fontsize=7)
ax[2].set_ylabel("First-order index"); ax[2].legend(frameon=False, fontsize=6.5)
ax[2].set_title("(c)", loc="left", fontsize=9)
save(fig, "fig13_uq_gated")

# ------------------------------------------------------------------ Fig 14: measured alignment, heat and thickness
t = ok(rd("thermal_lit")); t = t[np.isclose(t.V, 1.0)].assign(dT=lambda d: d.Tmax - 300)
cases = ["direct (measured)", "ideal MIGS (DFT)", "interlayer (measured)", "unpinned"]
clab = ["direct\n$S$=0.05", "DFT\n$S$=0.18", "oxide\n$S$=0.5", "none\n$S$=1"]
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
for k, (geom, col, mk) in enumerate((("edge", BLUE, "o"), ("top50", ORANGE, "s"))):
    for i, c in enumerate(cases):
        d = t[(t["case"] == c) & (t.geom == geom)].set_index("variant").dT
        if "full" not in d:
            continue
        x = i + (k - 0.5) * 0.25
        lo, hi = d.get("contact_sink", d["full"]), d.get("R_TBR_3e-8", d["full"])
        ax[0].errorbar(x, d["full"], yerr=[[d["full"] - lo], [hi - d["full"]]], fmt=mk, color=col, ms=4.5, capsize=2,
                       mfc=col if geom == "edge" else "white", lw=1)
ax[0].set_yscale("log"); ax[0].set_xticks(range(4)); ax[0].set_xticklabels(clab, fontsize=7)
ax[0].set_ylabel(r"$\Delta T_{\max}$ at $V_D$=1 V (K)"); ax[0].set_title("(a) Pt, $V_G=-3$ V", loc="left", fontsize=9)
ax[0].legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="none", label="Edge"),
                      Line2D([], [], color=ORANGE, marker="s", mfc="white", ls="none", label="Top, 50 nm")],
             frameon=False, fontsize=6.5, loc="upper left")
th = ok(rd("thickness_lit")); th = th[np.isclose(th.VG, -3.0)]
lt = th[th.study == "literature"]
for c, ls in (("interlayer (measured)", "-"), ("direct (measured)", "--")):
    for geom, col, mk in (("edge", BLUE, "o"), ("top50", ORANGE, "s")):
        d = lt[(lt["case"] == c) & (lt.geom == geom)].sort_values("layers")
        ax[1].semilogy(d.layers, d.R_tot, color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white")
ax[1].set_xticks([3, 4, 6, 8, 12]); ax[1].set_xlabel("Layers"); ax[1].set_ylabel(r"On-state $R_{\rm tot}$ ($\Omega$cm)")
ax[1].legend(handles=[Line2D([], [], color=INK, ls="-", label="oxide interlayer"), Line2D([], [], color=INK, ls="--", label="direct")],
             frameon=False, fontsize=6.5, loc="upper left", ncol=2)
ax[1].set_ylim(1, 1e11)
ax[1].set_title("(b) measured alignment, Pt", loc="left", fontsize=9)
rf = th[th.study == "reference_S1"]
for frac, col in ((0.2, GREEN), (0.32, INK), (0.5, YELLOW)):
    for geom, mk, ls in (("edge", "o", "-"), ("top50", "s", "--")):
        d = rf[(np.isclose(rf.dEv_frac, frac)) & (rf.geom == geom)].sort_values("layers")
        ax[2].plot(d.layers, d.R_tot, color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white",
                   label=f"VB share {frac}" if geom == "top50" else None)
ax[2].set_xticks([3, 4, 6, 8, 12]); ax[2].set_xlabel("Layers"); ax[2].set_ylabel(r"On-state $R_{\rm tot}$ ($\Omega$cm)")
ax[2].legend(frameon=False, fontsize=6.5, loc="upper left")
ax[2].set_title("(c) reference, $\\Phi_M$=5.0 eV", loc="left", fontsize=9)
save(fig, "fig14_measured_alignment")

# ------------------------------------------------------------------ Fig 15: electron injection check and mesh study
kB = 8.617e-5
e = ok(rd("inse_electron"))
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
mc = {"In": BLUE, "Ti": GREEN, "Pd": ORANGE}
for metal, col in mc.items():
    for geom, mk, ls in (("edge", "o", "-"), ("top50", "s", "--")):
        d = e[(e.metal == metal) & (e.geom == geom) & np.isclose(e.VG, 0.0)].sort_values("T0")
        ax[0].plot(1 / (kB * d.T0), np.log(d.I / d.T0 ** 2), color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white")
        rows = []
        for vg, g in e[(e.metal == metal) & (e.geom == geom)].groupby("VG"):
            if len(g) >= 3:
                rows.append((vg, -np.polyfit(1 / (kB * g.T0.values), np.log(g.I.values / g.T0.values ** 2), 1)[0]))
        rows = np.array(rows)
        ax[1].plot(rows[:, 0], rows[:, 1], color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white",
                   label=metal if geom == "edge" else None)
ax[0].set_xlabel(r"$1/k_BT$ (eV$^{-1}$)"); ax[0].set_ylabel(r"$\ln(I/T^2)$"); ax[0].set_title("(a) $V_G$=0", loc="left", fontsize=9)
ax[1].axhline(0.30, color=GREY, ls=":", lw=1); ax[1].text(2.2, 0.31, "input 0.30 eV", fontsize=6.5, color=GREY)
ax[1].set_xlabel("$V_G$ (V)"); ax[1].set_ylabel("Arrhenius barrier (eV)")
ax[1].legend(handles=ax[1].get_legend_handles_labels()[0] + [Line2D([], [], color=INK, ls="-", marker="o", ms=3.5, label="edge"),
             Line2D([], [], color=INK, ls="--", marker="s", mfc="white", ms=3.5, label="top")], frameon=False, fontsize=6.5, ncol=2)
ax[1].set_title("(b) electrons, $S$=0.05", loc="left", fontsize=9)
ms = ok(rd("mesh_convergence"))
keys = [("table3_reference", -2.0, 1.0, "Table 3, 1 V"), ("interlayer (measured)", -3.0, 0.05, "oxide, on"),
        ("direct (measured)", -3.0, 0.05, "direct, on")]
for (case, vg, v, lab), col in zip(keys, (INK, GREEN, YELLOW)):
    for geom, mk, ls in (("edge", "o", "-"), ("top50", "s", "--")):
        d = ms[(ms["case"] == case) & (ms.geom == geom) & np.isclose(ms.VG, vg) & np.isclose(ms.V, v)].sort_values("mf")
        if len(d):
            ax[2].plot(d.mf, 100 * (d.R_tot / d.R_tot.iloc[0] - 1), color=col, ls=ls, marker=mk, ms=3.5,
                       mfc=col if geom == "edge" else "white", label=lab if geom == "edge" else None)
ax[2].axhline(0, color=GREY, lw=0.8)
ax[2].set_xlabel("Mesh factor"); ax[2].set_ylabel(r"$\Delta R_{\rm tot}$ vs factor 1 (%)"); ax[2].set_xticks([1, 2, 3])
ax[2].set_ylim(-9, 14); ax[2].legend(frameon=False, fontsize=6, loc="lower left", ncol=2, columnspacing=0.8); ax[2].set_title("(c) mesh", loc="left", fontsize=9)
save(fig, "fig15_electron_mesh")

# ------------------------------------------------------------------ Fig 16: comparison with Chen et al. 2019
try:
    c = pd.concat([rd("chen2019"), rd("chen2019_300K")], ignore_index=True)
    c = ok(c)
    A = c[c.part == "A_arrhenius"]
    rows = []
    for (case, metal, vg), g in A.groupby(["case", "metal", "VG"]):
        if len(g) >= 3:
            rows.append((case, metal, g.Phi_M.iloc[0], vg,
                         -np.polyfit(1 / (kB * g.T0.values), np.log(g.I.values / g.T0.values ** 2), 1)[0]))
    R = pd.DataFrame(rows, columns=["case", "metal", "Phi", "VG", "phiB"])
    R0 = R[R.VG == 0.0]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1.35, 1]))
    mk = dict(In="o", Ti="v", Cr="^", Pd="s")
    for case, col in (("control", BLUE), ("oxidized monolayer", ORANGE)):
        g = R0[R0["case"] == case]
        for _, r in g.iterrows():
            ax[0].plot(r.Phi, r.phiB, mk[r.metal], color=col, ms=5)
        gg = g[g.metal != "In"] if case != "control" else g
        sl, ic = np.polyfit(gg.Phi, gg.phiB, 1)
        xs = np.array([4.0, 5.25]); ax[0].plot(xs, sl * xs + ic, color=col, lw=1, ls="--")
        ax[0].text(4.05 if case == "control" else 4.75, 0.34 if case == "control" else 0.05,
                   f"$S_{{\\rm app}}$ = {sl:.2f}" + (" (excl. In)" if case != "control" else ""), color=col, fontsize=7)
    ax[0].plot(4.12, 0.065, "*", color=INK, ms=9, label="Ref. [11], In + oxide")
    ax[0].plot(5.12, 0.57, "D", color=INK, mfc="none", ms=6, label="Ref. [11], Pd + oxide")
    ax[0].axhline(0, color=GREY, lw=0.8)
    ax[0].set_xlabel(r"Metal work function $\Phi_M$ (eV)"); ax[0].set_ylabel("Arrhenius barrier at $V_G$=0 (eV)")
    h = [Line2D([], [], color=BLUE, marker="o", ls="none", label="control, $S$=0.05"),
         Line2D([], [], color=ORANGE, marker="o", ls="none", label="oxide, $S$=0.5")] + ax[0].get_legend_handles_labels()[0]
    ax[0].legend(handles=h, frameon=False, fontsize=6.3, loc="upper left")
    ax[0].set_title("(a) barriers, 300--400 K", loc="left", fontsize=9)
    b = c[(c.part == "B_oncurrent") & np.isclose(c.VG, 80.0)].set_index("case").I * 100   # A/cm -> uA/um
    xs = np.arange(2)
    ax[1].bar(xs - 0.2, [b["control"], b["oxidized monolayer"]], 0.38, color=[BLUE, ORANGE], label="this work ($L$=1 $\\mu$m)")
    ax[1].bar(xs + 0.2, [103 / 6, 103], 0.38, color="none", edgecolor=INK, hatch="//", label="Ref. [11]")
    ax[1].errorbar([0.2], [103 / 6], yerr=[[103 / 6 - 103 / 6], [103 / 3 - 103 / 6]], color=INK, capsize=2, lw=1)
    ax[1].set_yscale("log"); ax[1].set_xticks(xs); ax[1].set_xticklabels(["control", "oxide\ninterlayer"])
    ax[1].set_ylabel(r"$I_D$ at $V_G$=80 V, $V_D$=1 V ($\mu$A/$\mu$m)"); ax[1].legend(frameon=False, fontsize=6.3, loc="upper left")
    ax[1].set_ylim(5, 800); ax[1].set_title("(b) on-current, In contact", loc="left", fontsize=9)
    save(fig, "fig16_chen2019")
except Exception as e:  # keep the other figures if the Chen data are missing
    print("fig16 skipped:", e)
