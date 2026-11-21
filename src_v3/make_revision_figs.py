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
fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.96), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
for (geom, var), col, ls, lab in [(("edge", "shared"), BLUE, "-", "Edge"), (("top50", "shared"), ORANGE, "-", "Top, 50 nm"),
                                  (("edge", "indep_edge"), BLUE, ":", r"Edge, indep. $S$")]:
    y = np.sort(np.log10(u[(u.geom == geom) & (u.variant == var)].R_tot.values))
    ax[0].plot(y, np.arange(1, len(y) + 1) / len(y), color=col, ls=ls, label=lab)
ax[0].set_xlabel(r"$\log_{10}R_{\rm tot}$ ($\Omega$cm)"); ax[0].set_ylabel("Cumulative fraction")
ax[0].legend(frameon=False, loc="lower right", fontsize=6.5, handlelength=1.6); ax[0].set_title("(a)", loc="left", fontsize=9)
sh = u[u.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
ie = u[u.variant == "indep_edge"].set_index("sample").R_tot
ax[1].loglog(sh["top50"], sh["edge"], "o", color=BLUE, ms=2.8, alpha=0.7, label="shared $S$")
j = sh.index.intersection(ie.index)
ax[1].loglog(sh.loc[j, "top50"], ie.loc[j], "^", color=GREEN, ms=2.8, alpha=0.6, mfc="none", label=r"independent $S_{\rm edge}$")
lim = [min(sh.min().min(), ie.min()) / 3, sh["top50"].max() * 3]
ax[1].plot(lim, lim, color=INK, lw=0.8); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel(r"$R_{\rm tot}^{\rm top}$ ($\Omega$cm)"); ax[1].set_ylabel(r"$R_{\rm tot}^{\rm edge}$ ($\Omega$cm)")
ax[1].legend(frameon=False, loc="upper left", fontsize=7); ax[1].set_title("(b)", loc="left", fontsize=9)
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
ax[2].set_xticks(xs); ax[2].set_xticklabels([lab[p] for p in par], rotation=90, fontsize=7.5)
ax[2].set_ylabel("First-order index"); ax[2].legend(frameon=False, fontsize=7)
ax[2].set_title("(c)", loc="left", fontsize=9)
save(fig, "fig13_uq_gated")

# ------------------------------------------------------------------ Fig 14: measured alignment, heat and thickness
t = ok(rd("thermal_lit")); t = t[np.isclose(t.V, 1.0)].assign(dT=lambda d: d.Tmax - 300)
cases = ["direct (measured)", "ideal MIGS (DFT)", "interlayer (measured)", "unpinned"]
clab = ["direct\n0.05", "DFT\n0.18", "oxide\n0.5", "none\n1"]
fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.96))
for k, (geom, col, mk) in enumerate((("edge", BLUE, "o"), ("top50", ORANGE, "s"))):
    for i, c in enumerate(cases):
        d = t[(t["case"] == c) & (t.geom == geom)].set_index("variant").dT
        if "full" not in d:
            continue
        x = i + (k - 0.5) * 0.25
        lo, hi = d.get("contact_sink", d["full"]), d.get("R_TBR_3e-8", d["full"])
        ax[0].errorbar(x, d["full"], yerr=[[d["full"] - lo], [hi - d["full"]]], fmt=mk, color=col, ms=4.5, capsize=2,
                       mfc=col if geom == "edge" else "white", lw=1)
ax[0].set_yscale("log"); ax[0].set_ylim(1e-4, 1e4); ax[0].set_xlabel("Pinning factor $S$"); ax[0].set_xticks(range(4)); ax[0].set_xticklabels(clab, fontsize=7)
ax[0].set_ylabel(r"$\Delta T_{\max}$ at $V_D$=1 V (K)"); ax[0].set_title("(a) Pt, $V_G=-3$ V", loc="left", fontsize=9)
ax[0].legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="none", label="Edge"),
                      Line2D([], [], color=ORANGE, marker="s", mfc="white", ls="none", label="Top, 50 nm")],
             frameon=False, fontsize=7, loc="upper left")
th = ok(rd("thickness_lit")); th = th[np.isclose(th.VG, -3.0)]
lt = th[th.study == "literature"]
for c, ls in (("interlayer (measured)", "-"), ("direct (measured)", "--")):
    for geom, col, mk in (("edge", BLUE, "o"), ("top50", ORANGE, "s")):
        d = lt[(lt["case"] == c) & (lt.geom == geom)].sort_values("layers")
        ax[1].semilogy(d.layers, d.R_tot, color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white")
ax[1].set_xticks([3, 4, 6, 8, 12]); ax[1].set_xlabel("Layers"); ax[1].set_ylabel(r"On-state $R_{\rm tot}$ ($\Omega$cm)")
ax[1].legend(handles=[Line2D([], [], color=INK, ls="-", label="oxide interlayer"), Line2D([], [], color=INK, ls="--", label="direct")],
             frameon=False, fontsize=7, loc="upper left", ncol=1)
ax[1].set_ylim(1, 1e12)
ax[1].set_title("(b) Measured, Pt", loc="left", fontsize=9)
rf = th[th.study == "reference_S1"]
for frac, col in ((0.2, GREEN), (0.32, INK), (0.5, YELLOW)):
    for geom, mk, ls in (("edge", "o", "-"), ("top50", "s", "--")):
        d = rf[(np.isclose(rf.dEv_frac, frac)) & (rf.geom == geom)].sort_values("layers")
        ax[2].plot(d.layers, d.R_tot, color=col, ls=ls, marker=mk, ms=3.5, mfc=col if geom == "edge" else "white",
                   label=f"VB share {frac}" if geom == "top50" else None)
ax[2].set_xticks([3, 4, 6, 8, 12]); ax[2].set_xlabel("Layers"); ax[2].set_ylabel(r"On-state $R_{\rm tot}$ ($\Omega$cm)")
ax[2].set_ylim(0.1, 2.6); ax[2].legend(frameon=False, fontsize=6.5, loc="upper left")
ax[2].set_title("(c) Reference, 5 eV", loc="left", fontsize=9)
save(fig, "fig14_measured_alignment")

# ------------------------------------------------------------------ Fig 15: electron injection check and mesh study
kB = 8.617e-5
e = ok(rd("inse_electron"))
fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.96))
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
ax[1].axhline(0.30, color=GREY, ls=":", lw=1); ax[1].text(4.0, 0.31, "input 0.30 eV", fontsize=7, color=GREY, ha="right"); ax[1].set_ylim(-0.1, 0.46)
ax[1].set_xlabel("$V_G$ (V)"); ax[1].set_ylabel("Arrhenius barrier (eV)")
ax[1].legend(handles=ax[1].get_legend_handles_labels()[0] + [Line2D([], [], color=INK, ls="-", marker="o", ms=3.5, label="edge"),
             Line2D([], [], color=INK, ls="--", marker="s", mfc="white", ms=3.5, label="top")], frameon=False, fontsize=7, ncol=2, loc="upper left", columnspacing=0.8, handlelength=1.6)
ax[1].set_title("(b) electrons, $S$=0.05", loc="left", fontsize=9)
ms = ok(rd("mesh_convergence"))
keys = [("table3_reference", -2.0, 1.0, "Table 3"), ("interlayer (measured)", -3.0, 0.05, "oxide"),
        ("direct (measured)", -3.0, 0.05, "direct")]
for (case, vg, v, lab), col in zip(keys, (INK, GREEN, YELLOW)):
    for geom, mk, ls in (("edge", "o", "-"), ("top50", "s", "--")):
        d = ms[(ms["case"] == case) & (ms.geom == geom) & np.isclose(ms.VG, vg) & np.isclose(ms.V, v)].sort_values("mf")
        if len(d):
            ax[2].plot(d.mf, 100 * (d.R_tot / d.R_tot.iloc[0] - 1), color=col, ls=ls, marker=mk, ms=3.5,
                       mfc=col if geom == "edge" else "white", label=lab if geom == "edge" else None)
ax[2].axhline(0, color=GREY, lw=0.8)
ax[2].set_xlabel("Mesh factor"); ax[2].set_ylabel(r"$\Delta R_{\rm tot}$ vs factor 1 (%)"); ax[2].set_xticks([1, 2, 3])
ax[2].set_ylim(-13, 14); ax[2].legend(frameon=False, fontsize=6.5, loc="lower left", ncol=1, handlelength=1.6); ax[2].set_title("(c) mesh", loc="left", fontsize=9)
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
    fig, ax0 = plt.subplots(1, 1, figsize=(4.4, 3.2)); ax = [ax0]
    mk = dict(In="o", Ti="v", Cr="^", Pd="s")
    for case, col in (("control", BLUE), ("oxidized monolayer", ORANGE)):
        g = R0[R0["case"] == case]
        for _, r in g.iterrows():
            ax[0].plot(r.Phi, r.phiB, mk[r.metal], color=col, ms=5)
        gg = g[g.metal != "In"] if case != "control" else g
        sl, ic = np.polyfit(gg.Phi, gg.phiB, 1)
        xs = np.array([4.0, 5.25]); ax[0].plot(xs, sl * xs + ic, color=col, lw=1, ls="--")
        ax[0].text(4.05 if case == "control" else 4.75, 0.34 if case == "control" else 0.05,
                   f"$S_{{\\rm app}}$ = {sl:.2f}" + (" (excl. In)" if case != "control" else ""), color=col, fontsize=8)
    ax[0].plot(4.12, 0.065, "*", color=INK, ms=9, label="Chen et al., In + oxide")
    ax[0].plot(5.12, 0.57, "D", color=INK, mfc="none", ms=6, label="Chen et al., Pd + oxide")
    ax[0].axhline(0, color=GREY, lw=0.8)
    ax[0].set_xlabel(r"Metal work function $\Phi_M$ (eV)"); ax[0].set_ylabel("Arrhenius barrier at $V_G$=0 (eV)")
    h = [Line2D([], [], color=BLUE, marker="o", ls="none", label="control, $S$=0.05"),
         Line2D([], [], color=ORANGE, marker="o", ls="none", label="oxide, $S$=0.5")] + ax[0].get_legend_handles_labels()[0]
    ax[0].legend(handles=h, frameon=False, fontsize=7, loc="upper left")
    ax[0].set_title("Arrhenius barriers, 300–400 K", loc="left", fontsize=9)
    save(fig, "fig16_chen2019")
except Exception as e:  # keep the other figures if the Chen data are missing
    print("fig16 skipped:", e)

# ------------------------------------------------------------------ Fig 17: hole-mass anisotropy decides the geometry
def _crossover(r):
    lr = np.log(r.values)
    k = np.where((lr[:-1] > 0) & (lr[1:] <= 0))[0]
    if not len(k):
        return None
    i = k[0]
    return r.index[i] + (r.index[i + 1] - r.index[i]) * lr[i] / (lr[i] - lr[i + 1])


try:
    km = ok(rd("kane_map"))
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1]))
    order = [("phiB=0.2", "0.2 eV", GREEN), ("unpinned", "Pt, unpinned", BLUE), ("interlayer (measured)", "Pt, interlayer", YELLOW),
             ("phiB=0.6", "0.6 eV", GREY), ("phiB=1.0", "1.0 eV", INK), ("direct (measured)", "Pt, direct", ORANGE)]
    variants = [("parabolic", "s", "parabolic, normal path"), ("kane", "^", "two-band, normal path"),
                ("kane_oblique", "o", "two-band, best straight path")]
    xo = {}
    for case, lab_, col in order:
        d = km[km["case"] == case]
        pb = float(d[d.geom == "edge"].phiBp.iloc[0])
        for var, mk_, _ in variants:
            g = d[d.variant == var].pivot_table(index="m_par", columns="geom", values="R_tot").dropna()
            r = (g["top50"] / g["edge"]).sort_index()
            if var == "kane_oblique":
                ax[0].semilogy(r.index, r.values, color=col, ls="--" if case.startswith("phiB") else "-", marker="o", ms=3.5, label=lab_)
            xo.setdefault(var, []).append((pb, _crossover(r)))
    ax[0].axhline(1, color=GREY, lw=0.9)
    ax[0].axvspan(0.85, 0.95, color=GREY, alpha=0.2, lw=0)
    ax[0].text(0.9, 40, "calc.", ha="center", fontsize=7, color=INK)
    ax[0].text(0.5, 1.12, "edge better", fontsize=7, color=INK, va="bottom")
    ax[0].text(0.5, 0.89, "top better", fontsize=7, color=INK, va="top")
    ax[0].set_xlabel(r"In-plane hole mass $m_{\parallel}$ ($m_0$)")
    ax[0].set_ylabel(r"$R_{\rm tot}^{\rm top}/R_{\rm tot}^{\rm edge}$ (on state)")
    ax[0].legend(frameon=False, fontsize=6.5, loc="lower left", ncol=2, columnspacing=0.6, handlelength=1.6)
    ax[0].set_ylim(5e-3, 1e2)
    ax[0].set_title("(a) two-band, best straight path", loc="left", fontsize=9)
    for var, mk_, lab_ in variants:
        pts = sorted(xo[var])
        xs1 = [p for p, c in pts if c is not None]; ys1 = [c for p, c in pts if c is not None]
        ax[1].plot(xs1, ys1, marker=mk_, color=INK if var == "kane_oblique" else GREY, ls="-" if var == "kane_oblique" else ":",
                   ms=4.5, mfc="white" if var == "parabolic" else None, label=lab_)
        for p, c in pts:
            if c is None and var == "kane_oblique":
                ax[1].annotate("", xy=(p, 1.33), xytext=(p, 1.15), arrowprops=dict(arrowstyle="->", color=INK, lw=1.0))
    ax[1].axhline(0.9, color=GREY, lw=0.9, ls="--"); ax[1].text(0.12, 0.92, "calculated, 6 layers", fontsize=7, color=GREY)
    ax[1].text(0.62, 1.3, "no crossover\nup to 1.2 $m_0$", fontsize=7, color=INK, va="center")
    ax[1].set_xlim(0.1, 1.15); ax[1].set_ylim(0.25, 1.4)
    ax[1].set_xlabel(r"Hole barrier $\phi_{Bp}$ (eV)"); ax[1].set_ylabel(r"Crossover $m_{\parallel}$ ($m_0$)")
    ax[1].legend(frameon=False, fontsize=6.5, loc="lower left")
    ax[1].set_title("(b) crossover mass", loc="left", fontsize=9)
    save(fig, "fig17_mass_anisotropy")
    for var in xo:
        print("crossover", var, [(round(p, 2), None if c is None else round(c, 2)) for p, c in sorted(xo[var])])
except Exception as e:
    print("fig17 skipped:", e)

# ------------------------------------------------------------------ Fig 18: gated Monte Carlo with direction-dependent hole mass
try:
    ud = ok(rd("uq_kane"))
    ud = ud[np.isclose(ud.VG, -3.0)] if "VG" in ud else ud
    sh = ud[ud.variant == "shared"].pivot_table(index="sample", columns="geom", values="R_tot").dropna()
    mi = ud.groupby("sample").m_in.first()
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.96), gridspec_kw=dict(width_ratios=[1, 1.15, 1]))
    for col, c_, lab_ in (("edge", BLUE, "Edge"), ("top50", ORANGE, "Top, 50 nm")):
        y = np.sort(np.log10(sh[col].values))
        ax[0].plot(y, np.arange(1, len(y) + 1) / len(y), color=c_, label=lab_)
    ax[0].set_xlabel(r"$\log_{10}R_{\rm tot}$ ($\Omega$cm)"); ax[0].set_ylabel("Cumulative fraction")
    ax[0].legend(frameon=False, fontsize=7, loc="lower right"); ax[0].set_title("(a)", loc="left", fontsize=9)
    sc = ax[1].scatter(sh["top50"], sh["edge"], c=mi.reindex(sh.index), cmap="viridis", s=8, vmin=0.3, vmax=1.5)
    lim = [sh.min().min() / 3, sh.max().max() * 3]
    ax[1].plot(lim, lim, color=INK, lw=0.8); ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"$R_{\rm tot}^{\rm top}$ ($\Omega$cm)"); ax[1].set_ylabel(r"$R_{\rm tot}^{\rm edge}$ ($\Omega$cm)")
    cb = fig.colorbar(sc, ax=ax[1], pad=0.02, fraction=0.06); cb.set_label(r"$m_\parallel$ ($m_0$)", fontsize=8); cb.ax.tick_params(labelsize=7)
    ax[1].set_title("(b)", loc="left", fontsize=9)
    r = (sh["edge"] / sh["top50"])
    edges_ = np.array([0.3, 0.45, 0.6, 0.8, 1.05, 1.5])
    xb, fb, nb = [], [], []
    for lo, hi in zip(edges_[:-1], edges_[1:]):
        sel = (mi.reindex(r.index) >= lo) & (mi.reindex(r.index) < hi)
        if sel.sum():
            xb.append(np.sqrt(lo * hi)); fb.append(100 * (r[sel] < 1).mean()); nb.append(int(sel.sum()))
    ax[2].plot(xb, fb, "o-", color=INK, ms=4)
    ax[2].axhline(50, color=GREY, lw=0.8, ls="--"); ax[2].axvspan(0.85, 0.95, color=GREY, alpha=0.2, lw=0)
    ax[2].set_xscale("log"); ax[2].set_xticks([0.3, 0.5, 1.0, 1.5]); ax[2].set_xticklabels(["0.3", "0.5", "1", "1.5"]); ax[2].xaxis.set_minor_formatter(plt.NullFormatter())
    ax[2].set_ylim(0, 100); ax[2].set_xlabel(r"$m_\parallel$ ($m_0$)"); ax[2].set_ylabel("Samples with edge < top (%)")
    ax[2].set_title("(c)", loc="left", fontsize=9)
    save(fig, "fig18_uq_kane")
    print("uq_kane n=", len(r), "edge<top %", round(100 * (r < 1).mean(), 1), "bins", list(zip([round(x, 2) for x in xb], [round(f, 1) for f in fb], nb)))
except Exception as e:
    print("fig18 skipped:", e)

# ------------------------------------------------------------------ Fig 19: geometry phase diagram in (m_par, m_c)
try:
    mp = ok(rd("mass_phase_fine"))
    if mp is None or not len(mp):
        mp = ok(rd("mass_phase"))
    cases = [("unpinned", "(a) Pt, unpinned"), ("interlayer (measured)", "(b) Pt, interlayer"),
             ("phiB=0.8", "(c) 0.8 eV"), ("direct (measured)", "(d) Pt, direct")]
    fig, ax = plt.subplots(1, 4, figsize=(6.5, 2.3), sharey=True)
    from matplotlib.colors import TwoSlopeNorm
    norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    for k, (case, ttl) in enumerate(cases):
        d = mp[mp["case"] == case].pivot_table(index=["m_c", "m_par"], columns="geom", values="R_tot").dropna()
        r = np.log10(d["top50"] / d["edge"]).unstack("m_par")
        mc = r.index.values; mpar = r.columns.values
        pb = float(mp[(mp["case"] == case) & (mp.geom == "edge")].phiBp.iloc[0])
        im = ax[k].pcolormesh(np.log10(mpar), np.log10(mc), np.clip(r.values, -3, 3), cmap="RdBu", norm=norm, shading="nearest")
        if np.nanmin(r.values) < 0 < np.nanmax(r.values):
            ax[k].contour(np.log10(mpar), np.log10(mc), r.values, levels=[0], colors=INK, linewidths=1.2)
        for (mc_, mp_), mk_ in (((0.17, 0.9), "*"), ((0.078, 0.9), "o")):
            ax[k].plot(np.log10(mp_), np.log10(mc_), mk_, color="k", mfc="white", ms=7 if mk_ == "*" else 5)
        ax[k].plot(np.log10([0.3, 1.9]), np.log10([0.3, 1.9]), color=GREY, lw=0.8, ls=":")
        qthr = {"interlayer (measured)": 11.5, "phiB=0.8": 2.85, "direct (measured)": 2.45}.get(case)
        if qthr:  # constant mass ratio m_par/m_c at the threshold
            xs_ = np.array([0.3, 1.9]); ax[k].plot(np.log10(xs_), np.log10(xs_ / qthr), color=INK, lw=0.9, ls="--")
            ax[k].text(0.97, 0.03, f"{qthr:g}", transform=ax[k].transAxes, fontsize=6, color=INK,
                       ha="right", va="bottom", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
        ax[k].set_xlim(np.log10(0.27), np.log10(2.1)); ax[k].set_ylim(np.log10(mc.min()) - 0.08, np.log10(mc.max()) + 0.08)
        ax[k].set_xticks(np.log10([0.3, 0.6, 0.9, 1.9])); ax[k].set_xticklabels(["0.3", "0.6", "0.9", "1.9"], fontsize=7); ax[k].grid(False)
        ax[k].set_title(f"{ttl}\n$\\phi_{{Bp}}$={pb:.2f} eV", loc="left", fontsize=8)
        ax[k].set_xlabel(r"$m_\parallel$ ($m_0$)", fontsize=8)
    yt = [v for v in (0.078, 0.12, 0.17, 0.25, 0.4, 0.6, 0.9, 2.5) if v <= mc.max() * 1.01]
    ax[0].set_yticks(np.log10(yt)); ax[0].set_yticklabels([f"{v:g}" for v in yt], fontsize=7)
    ax[0].set_ylabel(r"$m_c$ ($m_0$)", fontsize=8)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03); cb.set_label(r"$\log_{10}(R^{\rm top}/R^{\rm edge})$", fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.savefig(os.path.join(F, "fig19_mass_phase.png"), bbox_inches="tight"); fig.savefig(os.path.join(F, "fig19_mass_phase.pdf"), bbox_inches="tight")
    plt.close(fig); print("saved fig19_mass_phase")
except Exception as e:
    print("fig19 skipped:", e)
