"""
Regenerates fig_channel_length.png on log-log axes from the already-computed
channel_length_study.csv data. A log-log plot is the correct representation
here because the manuscript's SCLC discussion (Section 3.6) is built entirely
around power-law exponents (R ~ L^n); on linear axes a power law is hard to
read visually, while on log-log axes it becomes a straight line whose slope
is the exponent, directly supporting the text's numerical claims.

Does not re-run DEVSIM -- reads the existing CSV only.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

data_dir = r"c:\Users\harsh\ドキュメント\InSe"

def main():
    df = pd.read_csv(os.path.join(data_dir, "channel_length_study.csv"))
    df = df.dropna(subset=["L_gap_nm", "R_Top", "R_Edge"]).sort_values("L_gap_nm")

    L = df["L_gap_nm"].values
    R_top = df["R_Top"].values
    R_edge = df["R_Edge"].values

    fig, ax = plt.subplots(figsize=(7, 5.5))

    ax.plot(L, R_top, marker='^', markersize=9, linestyle='--', linewidth=2,
            color='red', label='Top Contact (50 nm)')
    ax.plot(L, R_edge, marker='o', markersize=8, linestyle='-', linewidth=2,
            color='blue', label='Edge Contact')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Channel Length $L_{ch}$ (nm)')
    ax.set_ylabel(r'Total Resistance ($\Omega\cdot$cm)')
    ax.set_title(r'Total Resistance vs Channel Length ($\Phi_M$ = 5.0 eV)')
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(frameon=False)

    # Annotate local exponent for edge contact between last two points as a
    # concrete visual anchor for the ~3.3 exponent quoted in the text.
    n_local = np.log(R_edge[-1] / R_edge[-2]) / np.log(L[-1] / L[-2])
    ax.annotate(f"local slope $n$ = {n_local:.2f}",
                xy=(L[-2], R_edge[-2]), xytext=(0.55, 0.15),
                textcoords='axes fraction', fontsize=10,
                arrowprops=dict(arrowstyle='->', color='gray'))

    fig.tight_layout()
    out_path = os.path.join(data_dir, "fig_channel_length.png")
    plt.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")

if __name__ == "__main__":
    main()
