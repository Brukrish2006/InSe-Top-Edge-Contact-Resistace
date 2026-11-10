"""2D field maps (hole density, temperature rise) of gated top and edge devices at 1 V.
Runs the DEVSIM model, stores fields in ../data_v3/fields_*.csv and draws ../figures_v3/fig11_field_maps.*"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import tcad2d as m  # noqa

HERE = os.path.dirname(__file__)
D = os.path.join(HERE, "..", "data_v3"); F = os.path.join(HERE, "..", "figures_v3")
NM = 1e-7
GATE = dict(t_ox=10e-7)
plt.rcParams.update({"font.size": 9, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})

fields = {}
for kind in ("top", "edge"):
    f = os.path.join(D, f"fields_{kind}.csv")
    if not os.path.exists(f):
        dv = m.Device(kind, Phi_M=5.0, Lc=50 * NM, gate=GATE)
        dv.set_gate(-2.0); dv.set_bias(1.0)
        I, Ip, In, It = dv.current()
        df = pd.DataFrame(dict(x_nm=dv.x / NM, y_nm=dv.y / NM, p=dv.nodev("Holes"), T=dv.nodev("T")))
        df.attrs["geo"] = dv.geo
        df.to_csv(f, index=False)
        with open(os.path.join(D, f"fields_{kind}_meta.txt"), "w") as fh:
            fh.write(f"I_A_per_cm={I}\nLc_nm={dv.geo['Lc'] / NM}\nLch_nm={dv.geo['Lch'] / NM}\nx_ch0_nm={dv.geo['x_ch0'] / NM}\nx_ch1_nm={dv.geo['x_ch1'] / NM}\n")
    fields[kind] = pd.read_csv(f)

def meta(kind):
    return dict(l.strip().split("=") for l in open(os.path.join(D, f"fields_{kind}_meta.txt")))

fig, ax = plt.subplots(2, 2, figsize=(6.5, 4.45), sharex=False)
for r, (kind, name) in enumerate((("top", "Top contact"), ("edge", "Edge contact"))):
    d = fields[kind]; mt = meta(kind)
    x0, x1 = float(mt["x_ch0_nm"]), float(mt["x_ch1_nm"])
    for c, (col, lab, cmap) in enumerate((("p", "$\\log_{10}\\,p$ (cm$^{-3}$)", "viridis"), ("T", "$\\Delta T$ (K)", "inferno"))):
        v = np.log10(np.maximum(d.p, 1.0)) if col == "p" else d["T"] - 300.0
        cs = ax[r, c].tricontourf(d.x_nm, d.y_nm, v, 40, cmap=cmap)
        cb = fig.colorbar(cs, ax=ax[r, c], pad=0.02, fraction=0.05); cb.set_label(lab, fontsize=8); cb.ax.tick_params(labelsize=8)
        if kind == "top":
            ax[r, c].axvline(x0, color="w", lw=0.6, ls=":"); ax[r, c].axvline(x1, color="w", lw=0.6, ls=":")
        ax[r, c].invert_yaxis(); ax[r, c].set_ylabel("$y$ (nm), oxide at bottom", fontsize=8); ax[r, c].tick_params(labelsize=8)
        ax[r, c].set_title(f"({'abcd'[2 * r + c]}) {name}, " + ("hole density" if c == 0 else "temperature"), fontsize=8, loc="left")
for c in range(2):
    ax[0, c].set_xlabel("$x$ (nm)", fontsize=8); ax[1, c].set_xlabel("$x$ (nm)", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(F, "fig11_field_maps.png")); fig.savefig(os.path.join(F, "fig11_field_maps.pdf"))
print("saved fig11_field_maps")
