import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Set global plot style for publication quality
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

data_dir = r"c:\Users\harsh\ドキュメント\InSe"

def plot_iv_curves():
    """Plot I-V curves for both Edge and Top (L_c = 50 nm) contacts at different Work Functions."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    wfs = ['4.33', '4.8', '5.0', '5.2']
    colors = ['#d62728', '#ff7f0e', '#1f77b4', '#2ca02c']
    markers = ['o', 's', '^', 'D']

    for wf, color, marker in zip(wfs, colors, markers):
        fname_edge = os.path.join(data_dir, f"iv_edge_WF{wf}eV.csv")
        if os.path.exists(fname_edge):
            df = pd.read_csv(fname_edge)
            I_drain = np.abs(df['I_drain (A/cm)'])
            ax.plot(df['V_drain (V)'], I_drain, label=f'Edge, WF = {wf} eV',
                    color=color, linestyle='-', linewidth=2,
                    marker=marker, markersize=5, markevery=4)

        fname_top = os.path.join(data_dir, f"iv_top_Lc50nm_WF{wf}eV.csv")
        if os.path.exists(fname_top):
            df = pd.read_csv(fname_top)
            I_drain = np.abs(df['I_drain (A/cm)'])
            ax.plot(df['V_drain (V)'], I_drain, label=f'Top ($L_c$=50nm), WF = {wf} eV',
                    color=color, linestyle='--', linewidth=2,
                    marker=marker, markersize=5, markevery=4, markerfacecolor='none')

    ax.set_yscale('log')
    ax.set_xlabel('Drain Voltage (V)')
    ax.set_ylabel('|Drain Current| (A/cm)')
    ax.set_title('I-V Characteristics: Edge vs. Top Contact ($L_c$ = 50 nm)')
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax.grid(True, which="both", ls="--", alpha=0.2)

    plt.savefig(os.path.join(data_dir, 'fig1_iv_curves.png'))
    plt.close()

def plot_resistance_vs_wf():
    """Plot Contact Resistance vs Work Function for Edge and Top contacts."""
    df = pd.read_csv(os.path.join(data_dir, 'summary_results.csv'))

    fig, ax = plt.subplots(figsize=(7, 5))

    # Filter edge
    edge_df = df[df['Geometry'] == 'edge']
    ax.plot(edge_df['Work Function (eV)'], edge_df['R_total (Ohm*cm)'],
            marker='o', markersize=8, linestyle='-', linewidth=2, label='Edge Contact', color='#d62728')

    # Filter top L_c = 50nm (5e-6 cm)
    top_df = df[(df['Geometry'] == 'top') & (np.isclose(df['L_c (cm)'].astype(float), 5e-6))]
    ax.plot(top_df['Work Function (eV)'], top_df['R_total (Ohm*cm)'],
            marker='s', markersize=8, linestyle='--', linewidth=2, label=r'Top Contact ($L_c = 50$ nm)', color='#1f77b4')

    ax.set_yscale('log')
    ax.set_xlabel('Metal Work Function (eV)')
    ax.set_ylabel(r'Total Resistance ($\Omega \cdot$ cm)')
    ax.set_title('Resistance vs. Metal Work Function')
    ax.legend(frameon=False)
    ax.grid(True, which="both", ls="--", alpha=0.2)

    plt.savefig(os.path.join(data_dir, 'fig2_resistance_vs_wf.png'))
    plt.close()

def plot_resistance_vs_lc():
    """Plot Contact Resistance vs Contact Length for Top Contact, across all four work functions."""
    df = pd.read_csv(os.path.join(data_dir, 'summary_results.csv'))

    fig, ax = plt.subplots(figsize=(8, 5.5))
    wfs = [4.33, 4.8, 5.0, 5.2]
    colors = ['#d62728', '#ff7f0e', '#1f77b4', '#2ca02c']
    markers = ['o', 's', '^', 'D']

    for wf, color, marker in zip(wfs, colors, markers):
        # Filter top contact for specific WF
        top_df = df[(df['Geometry'] == 'top') & (np.isclose(df['Work Function (eV)'], wf))]
        # Convert L_c to nm for plotting
        L_c_nm = top_df['L_c (cm)'].astype(float) * 1e7
        ax.plot(L_c_nm, top_df['R_total (Ohm*cm)'],
                marker=marker, markersize=8, linestyle='-', linewidth=2, label=f'Top (WF = {wf} eV)', color=color)

        # Add edge contact resistance as a horizontal reference line
        edge_df = df[(df['Geometry'] == 'edge') & (np.isclose(df['Work Function (eV)'], wf))]
        if not edge_df.empty:
            edge_r = edge_df['R_total (Ohm*cm)'].iloc[0]
            ax.axhline(y=edge_r, linestyle=':', color=color, alpha=0.7, label=f'Edge Limit (WF = {wf} eV)')

    ax.set_yscale('log')
    ax.set_xlabel('Contact Length, $L_c$ (nm)')
    ax.set_ylabel(r'Total Resistance ($\Omega \cdot$ cm)')
    ax.set_title('Top Contact Resistance vs. Contact Length (All Work Functions)')
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    ax.grid(True, which="both", ls="--", alpha=0.3)
    fig.tight_layout()

    plt.savefig(os.path.join(data_dir, 'fig3_resistance_vs_lc.png'))
    plt.close()

if __name__ == "__main__":
    print("Generating plots...")
    plot_iv_curves()
    plot_resistance_vs_wf()
    plot_resistance_vs_lc()
    print("Plots saved successfully to the directory.")
