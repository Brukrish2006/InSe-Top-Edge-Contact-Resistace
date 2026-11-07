"""Figure: contact resistance with the measured / calculated pinning of InSe and the
literature band alignment (chi = 4.6 eV). Reads data_v3/barrier_sweep.csv and
data_v3/metals_pinning.csv (written by run_parallel.py); writes figures_v3/fig12_*."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = os.path.join(os.path.dirname(__file__), "..", "data_v3")
F = os.path.join(os.path.dirname(__file__), "..", "figures_v3")
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "lines.linewidth": 1.6, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e4e0",
                     "grid.linewidth": 0.6, "savefig.dpi": 300, "axes.axisbelow": True})
BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#3a3a38"

b = pd.read_csv(os.path.join(D, "barrier_sweep.csv"))
mt = pd.read_csv(os.path.join(D, "metals_pinning.csv"))

fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.1), gridspec_kw=dict(width_ratios=[1, 0.75, 1.25]))

# (a) on-state and ungated resistance versus hole barrier
on = b[(b.study == "gated") & np.isclose(b.VG, -3.0)]
for geom, col, mk, lab in (("edge", BLUE, "o", "Edge"), ("top50", ORANGE, "s", "Top, 50 nm")):
    d = on[(on.geom == geom) & on.tunneling.isna()].sort_values("phiB_set")
    ax[0].semilogy(d.phiB_set, d.R_tot, color=col, marker=mk, ls="-", label=lab)
    d = on[(on.geom == geom) & (on.tunneling == "off")].sort_values("phiB_set")
    ax[0].semilogy(d.phiB_set, d.R_tot, color=col, marker=mk, ls=":", mfc="white")
ax[0].set_xlabel(r"Hole barrier $\phi_{Bp}$ (eV)")
ax[0].set_ylabel(r"$R_{\rm tot}$ at $V_G=-3$ V ($\Omega$cm)")
h = ax[0].get_legend_handles_labels()[0] + [Line2D([], [], color=INK, ls="-", label="nonlocal tunnelling"),
                                            Line2D([], [], color=INK, ls=":", mfc="white", marker="o", label="no tunnelling")]
ax[0].legend(handles=h, loc="upper right", frameon=False, fontsize=6.3, ncol=2, columnspacing=0.6, handlelength=1.6, borderaxespad=0.1)
ax[0].set_ylim(1e-2, 1e16); ax[0].set_yticks([1e0, 1e4, 1e8, 1e12])
ax[0].set_title("(a)", loc="left", fontsize=9)

# (b) edge/top ratio of on-state resistance
e = on[(on.geom == "edge") & on.tunneling.isna()].set_index("phiB_set").R_tot
t = on[(on.geom == "top50") & on.tunneling.isna()].set_index("phiB_set").R_tot
r = (t / e).dropna()
ax[1].semilogy(r.index, r.values, color=INK, marker="D")
ax[1].set_xlabel(r"Hole barrier $\phi_{Bp}$ (eV)")
ax[1].set_ylabel(r"$R_{\rm tot}^{\rm top}/R_{\rm tot}^{\rm edge}$ (on state)")
ax[1].set_title("(b)", loc="left", fontsize=9)

# (c) real metals, on-state, four pinning cases
cases = ["unpinned", "interlayer (measured)", "ideal MIGS (DFT)", "direct (measured)"]
lab = {"unpinned": "1\nnone", "interlayer (measured)": "0.5\noxide (M)",
       "ideal MIGS (DFT)": "0.18\nDFT", "direct (measured)": "0.05\ndirect (M)"}
metals = ["Ti", "Cr", "Au", "Pd", "Pt"]
mk = dict(Ti="v", Cr="^", Au="o", Pd="s", Pt="D")
g = mt[mt.study == "gated"]
for i, c in enumerate(cases):
    for j, mtl in enumerate(metals):
        x = i + (j - 2) * 0.13
        for geom, col in (("edge", BLUE), ("top50", ORANGE)):
            d = g[(g["case"] == c) & (g.metal == mtl) & (g.geom == geom)]
            if len(d):
                ax[2].semilogy(x, d.R_tot.iloc[0], marker=mk[mtl], color=col, ls="none", ms=4.5,
                               mfc=col if geom == "edge" else "white")
ax[2].set_xticks(range(len(cases)))
ax[2].set_xticklabels([lab[c] for c in cases], fontsize=7)
ax[2].set_xlabel("Pinning factor $S$", labelpad=1)
ax[2].set_xlim(-0.5, 3.5)
ax[2].set_ylabel(r"$R_{\rm tot}$ at $V_G=-3$ V ($\Omega$cm)")
hm = [Line2D([], [], color=INK, marker=mk[x], ls="none", ms=4.5, label=x) for x in metals]
hm += [Line2D([], [], color=BLUE, marker="o", ls="none", label="edge (filled)"),
       Line2D([], [], color=ORANGE, marker="o", mfc="white", ls="none", label="top (open)")]
ax[2].legend(handles=hm, loc="upper center", frameon=False, ncol=4, fontsize=6, columnspacing=0.6, handletextpad=0.2, bbox_to_anchor=(0.5, 1.02))
ax[2].set_ylim(1e-1, 1e10)
ax[2].set_title("(c)", loc="left", fontsize=9)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"fig12_measured_pinning.{ext}"))
print("saved fig12_measured_pinning")
print("ratio", r.round(1).to_dict())
