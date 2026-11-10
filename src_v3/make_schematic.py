"""Fig. 0: device geometries (not to scale). Writes ../figures_v3/fig0_device_schematic.*"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

F = os.path.join(os.path.dirname(__file__), "..", "figures_v3")
plt.rcParams.update({"font.size": 8.5, "savefig.dpi": 300})
MET, INSE, OX, SI, INK = "#d8b24a", "#6fbf9b", "#d9d9d6", "#8f8f8b", "#3a3a38"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def stack(ax, x0, x1):
    ax.add_patch(Rectangle((x0, -1.2), x1 - x0, 0.5, fc=OX, ec=INK, lw=0.6))
    ax.add_patch(Rectangle((x0, -2.2), x1 - x0, 1.0, fc=SI, ec=INK, lw=0.6))
    ax.text((x0 + x1) / 2, -1.7, "p$^{++}$ Si back gate ($V_G$)", ha="center", va="center", fontsize=8, color="w")
    ax.text(x1 + 0.1, -0.95, "SiO$_2$, 10 nm", ha="left", va="center", fontsize=8)


def dim(ax, xa, xb, y, txt):
    ax.add_patch(FancyArrowPatch((xa, y), (xb, y), arrowstyle="<->", mutation_scale=7, lw=0.7, color=INK))
    ax.text((xa + xb) / 2, y + 0.12, txt, ha="center", va="bottom", fontsize=8)


fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.08))
# (a) top contact -----------------------------------------------------------------
a = ax[0]
stack(a, 0, 10)
a.add_patch(Rectangle((0, -0.7), 10, 0.7, fc=INSE, ec=INK, lw=0.6))
a.text(5, -0.35, "few-layer InSe ($t$ ≈ 5 nm, $p$-type)", ha="center", va="center", fontsize=8)
for xs in ((0, 3.2), (6.8, 10)):
    a.add_patch(Rectangle((xs[0], 0.12), xs[1] - xs[0], 0.7, fc=MET, ec=INK, lw=0.6))
a.text(2.0, 0.47, "source", ha="center", va="center", fontsize=8); a.text(8.4, 0.47, "drain", ha="center", va="center", fontsize=8)
a.annotate("vdW gap", xy=(1.6, 0.06), xytext=(4.9, 0.55), fontsize=8, ha="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
a.add_patch(FancyArrowPatch((0.5, 0.7), (0.5, -0.5), arrowstyle="-|>", mutation_scale=7, lw=1.2, color=ORANGE))
a.add_patch(FancyArrowPatch((0.5, -0.5), (2.6, -0.5), arrowstyle="-|>", mutation_scale=7, lw=1.2, color=ORANGE))
dim(a, 0, 3.2, 1.25, "$L_c$"); dim(a, 3.2, 6.8, 1.25, "$L_{ch}$")
a.set_xlim(-0.3, 12.2); a.set_ylim(-2.4, 2.1); a.axis("off")
a.set_title("(a) Top contact", fontsize=9, loc="left")
# (b) edge contact -----------------------------------------------------------------
b = ax[1]
stack(b, 0, 10)
b.add_patch(Rectangle((2.6, -0.7), 4.8, 0.7, fc=INSE, ec=INK, lw=0.6))
b.text(5, -0.35, "few-layer InSe", ha="center", va="center", fontsize=8)
b.add_patch(Rectangle((1.0, -0.7), 1.6, 1.5, fc=MET, ec=INK, lw=0.6)); b.add_patch(Rectangle((7.4, -0.7), 1.6, 1.5, fc=MET, ec=INK, lw=0.6))
b.text(1.8, 0.4, "source", ha="center", va="center", fontsize=8); b.text(8.2, 0.4, "drain", ha="center", va="center", fontsize=8)
b.add_patch(FancyArrowPatch((1.3, -0.35), (3.5, -0.35), arrowstyle="-|>", mutation_scale=7, lw=1.2, color=BLUE))
b.annotate("edge interface", xy=(2.6, -0.1), xytext=(4.9, 0.6), fontsize=8, ha="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
dim(b, 2.6, 7.4, 1.25, "$L_{ch}$")
b.set_xlim(-0.3, 12.2); b.set_ylim(-2.4, 2.1); b.axis("off")
b.set_title("(b) Edge contact", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(F, "fig0_device_schematic.png")); fig.savefig(os.path.join(F, "fig0_device_schematic.pdf"))
print("saved fig0")
