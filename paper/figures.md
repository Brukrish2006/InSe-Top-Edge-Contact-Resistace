# Research Paper Figures

The following figures have been generated from our finite-element simulations and are ready to be included in the final manuscript.

## Figure 1: I-V Characteristics
This figure demonstrates the fundamental transport behavior of both contact geometries (edge, and top at $L_c=50$ nm) under four work functions. The transition from Schottky-limited (low current) to injection-barrier-free (high current, space-charge-limited bulk transport — see Section 3.6) is clearly visible.

![Figure 1: I-V Curves](fig1_iv_curves.png)

## Figure 2: Contact Resistance vs. Work Function
This compares the total device resistance for both top and edge geometries. Edge contacts consistently outperform top contacts at every work function tested, most substantially once the injection barrier collapses (WF $\ge$ 5.0 eV), because they bypass the out-of-plane current crowding bottleneck. (Note: bulk transport in this high-WF regime is space-charge-limited rather than Ohmic; see Section 3.6.)

![Figure 2: Resistance vs Work Function](fig2_resistance_vs_wf.png)

## Figure 3: Resistance vs. Contact Length (Top Contact, All Work Functions)
This highlights the classic transfer-length limitation of top contacts, now shown across all four simulated work functions rather than a subset. The length-dependence is itself work-function dependent: minimal in the strongly Schottky-limited regime (4.33, 4.8 eV) and more pronounced once the barrier collapses (5.0, 5.2 eV). The edge-contact limit (dotted lines) represents the theoretically superior baseline at each work function.

![Figure 3: Resistance vs Contact Length](fig3_resistance_vs_lc.png)

## Figure 4: Multiphysics Thermal Analysis
This plot captures one of the paper's central findings. It shows how current crowding in top contacts leads to extreme, localized Joule heating. The left panel compares top ($L_c=50$ nm, the length underlying the headline claim) against edge contacts, where edge geometry drops the peak lattice temperature by roughly 94 K at the same bias (2.0 V, $\Phi_M=5.0$ eV).

![Figure 4: Thermal Comparison](thermal_comparison.png)

## Figure 5: Tunneling Enhancement ($\Gamma$) Sensitivity
This study sweeps the tunneling enhancement factor $\Gamma \in [1, 10, 100]$ across varying work functions. The overlapping curves prove that the device is bottlenecked by macroscopic drift-diffusion (either barrier $n_{eq}$ or bulk resistance) rather than purely interface emission velocity.

![Figure 5: Tunneling Sensitivity](fig_gamma_sensitivity_final.png)

## Figure 6: Explicit Van der Waals Gap Modeling
By introducing a vdW gap transmission penalty ($T_{vdW}$) to the top contact model in the Schottky regime, we verify that reducing interfacial velocity by up to a factor of 100 ($T_{vdW}=0.01$) does not change total resistance, confirming the macroscopic electrostatics dominate.

![Figure 6: vdW Gap Analysis](fig_vdw_gap_final.png)

## Figure 7: Channel Length Scaling
This sweeps the channel length ($L_{ch}$) from 10 nm to 100 nm, plotted on log-log axes so the power-law scaling is directly visible as a straight-line slope. Resistance scales superlinearly (>100x resistance jump over a 10x length increase, local exponent $n\approx3.3$ between 50-100 nm), consistent with the space-char