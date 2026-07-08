import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_thermal_results(csv_file):
    df = pd.read_csv(csv_file)

    df['L_c (cm)'] = pd.to_numeric(df['L_c (cm)'], errors='coerce')
    # Filter for top contact at L_c = 50 nm (5e-6 cm) -- matches the headlined
    # 94.3 K self-heating-reduction figure quoted in the Abstract/Results/Conclusion.
    df_top = df[(df['Geometry'] == 'top') & (np.isclose(df['L_c (cm)'], 5e-06))]
    df_edge = df[df['Geometry'] == 'edge']

    # Also we want to compare with L_c variation for Top contact at WF=5.0
    df_top_5 = df[(df['Geometry'] == 'top') & (df['Work Function (eV)'] == 5.0)]

    plt.figure(figsize=(12, 5))

    # Subplot 1: Max Temperature vs Work Function (Top vs Edge)
    plt.subplot(1, 2, 1)
    wf_edge = df_edge['Work Function (eV)'].values
    t_edge = df_edge['Max_T at 2.0V (K)'].values

    wf_top = df_top['Work Function (eV)'].values
    t_top = df_top['Max_T at 2.0V (K)'].values

    # Make sure we only plot where we have both
    common_wfs = sorted(list(set(wf_top) & set(wf_edge)))

    t_top_c = [t_top[np.where(wf_top == w)[0][0]] for w in common_wfs]
    t_edge_c = [t_edge[np.where(wf_edge == w)[0][0]] for w in common_wfs]

    x = np.arange(len(common_wfs))
    width = 0.35

    plt.bar(x - width/2, t_top_c, width, label='Top Contact (L_c=50nm)', color='#ff9999')
    plt.bar(x + width/2, t_edge_c, width, label='Edge Contact', color='#66b3ff')

    plt.ylabel('Max Lattice Temperature (K)')
    plt.xlabel('Work Function (eV)')
    plt.title('Self-Heating: Top vs Edge Contact @ V_ds=2.0V')
    plt.xticks(x, common_wfs)
    plt.ylim(290, max(max(t_top_c), max(t_edge_c)) + 20)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Subplot 2: Max Temperature vs L_c for Top Contact
    plt.subplot(1, 2, 2)
    # Convert L_c from string back to float, then to nm
    try:
        lc_nm = [float(x)*1e7 for x in df_top_5['L_c (cm)'].values]
        t_lc = df_top_5['Max_T at 2.0V (K)'].values

        # Sort by L_c
        sort_idx = np.argsort(lc_nm)
        lc_nm = np.array(lc_nm)[sort_idx]
        t_lc = np.array(t_lc)[sort_idx]

        plt.plot(lc_nm, t_lc, marker='o', linestyle='-', color='#ff9999', linewidth=2, markersize=8)
        plt.ylabel('Max Lattice Temperature (K)')
        plt.xlabel('Contact Length L_c (nm)')
        plt.title('Top Contact Self-Heating vs L_c (WF=5.0 eV)')
        plt.grid(linestyle='--', alpha=0.7)
    except Exception as e:
        print("Could not plot L_c dependence:", e)

    plt.tight_layout()
    plt.savefig('thermal_comparison.png', dpi=300)
    print("Saved thermal_comparison.png")

if __name__ == "__main__":
    plot_thermal_results("summary_results_thermal.csv")
