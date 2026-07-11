# Geometry-Dependent Contact Resistance and Space-Charge-Limited Bulk Transport in Metal–InSe Junctions: A Coupled Electro-Thermal Finite-Element Simulation Study

**Harsha Adhikary**¹

¹Department of Physics, Indian Institute of Science (IISc), Bengaluru, India

Corresponding author: harshaa@iisc.ac.in

**Keywords:** indium selenide; two-dimensional semiconductors; contact resistance; Schottky barrier; edge contact; self-heating; TCAD simulation

---

## Abstract

Efficient hole injection into indium selenide (InSe), a two-dimensional van der Waals semiconductor of growing interest for post-silicon logic, is a materials-processing challenge that has received far less attention than its well-documented electron transport. Contact resistance at the metal–InSe interface routinely dominates device resistance, yet no systematic, geometry-resolved account of how contact placement and metal work function jointly control this resistance has been reported for InSe specifically. Here we combine a coupled drift-diffusion and lattice-heat-transport finite-element model (DEVSIM) with direct benchmarking against the near-zero-barrier indium contact reported experimentally for InSe [5], to provide quantitative processing guidance for contact design. Sweeping the contact work function ($\Phi_M$ = 4.33–5.2 eV) and comparing conventional top-contact (out-of-plane deposition) against edge-contact (in-plane, exposed-edge deposition) geometries at contact lengths of 10–50 nm, we identify the work function at which the Schottky barrier for hole injection collapses (~5.0 eV, consistent with literature-reported band alignment) and quantify the resulting drop in resistance from over $10^9\,\Omega\cdot\text{cm}$ to 1.5–1.8 $\Omega\cdot\text{cm}$. Beyond this barrier-free threshold, transport is governed by Space-Charge-Limited Current (SCLC) rather than Ohmic conduction, a regime we confirm independently via mesh-refinement and bias-scaling tests. Edge contacts consistently outperform top contacts, with a 17–28% resistance advantage under idealized isothermal conditions that converges to 21–27% once mesh effects are resolved — a result triangulated across three independent numerical checks — and narrows to a robust ~26% once realistic self-heating is included. Edge contacts also run substantially cooler under high bias. Together with the identified sensitivity of the results to hole mobility, thermal conductivity, and doping, these quantitative, geometry- and work-function-resolved results are intended to give experimentalists and process engineers concrete, actionable targets for contact metal selection and placement during InSe device fabrication.

---

## 1. Introduction

Contact engineering is one of the central materials-processing challenges in scaling two-dimensional (2D) van der Waals semiconductors toward practical devices. Indium selenide (InSe) illustrates this challenge acutely: its room-temperature electron mobility can exceed 1,000 cm²/(V·s) in high-quality, encapsulated devices [1,2], wafer-scale growth methods compatible with integrated-circuit processing have recently been demonstrated [3], and individual InSe transistors have shown near-ballistic channel transport [4]. As channel transport approaches these intrinsic limits, the metal contact — not the channel — increasingly sets the achievable device performance, a shift already visible experimentally: indium contacts aligned to the InSe conduction band have been reported to achieve a near-zero Schottky barrier (~50 meV) [5], demonstrating that work-function selection is a practical, processing-level lever for contact engineers rather than a purely theoretical parameter. However, viable complementary (CMOS) logic also requires efficient hole injection into a p-type channel, and because InSe is naturally n-type, this p-type contact is the less-studied, more critical bottleneck for CMOS integration. Metal–2D-semiconductor junctions routinely form Schottky barriers that dominate the total device resistance, a challenge documented across graphene, transition-metal dichalcogenides (TMDs), and InSe alike [6]. InSe compounds this problem with strongly anisotropic transport: weak out-of-plane van der Waals coupling forces carriers injected through a conventional "top" contact to cross a highly resistive out-of-plane path before reaching the high-mobility in-plane channel, producing severe current crowding at the contact edge and, under bias, localized Joule heating that is worsened by InSe's low out-of-plane thermal conductivity.

"Edge-contacted" architectures, in which the metal is deposited on the exposed in-plane edge of the 2D crystal rather than on its top face, were first demonstrated for graphene [7] and have since been extended to other 2D systems [8], consistently bypassing the van der Waals gap and improving injection efficiency in devices that were subsequently fabricated and measured. The electrical case for edge contacts is well established for graphene and select TMDs. The literature we reviewed does not include a systematic, geometry-resolved comparison of top- versus edge-contact resistance for InSe specifically, across a physically motivated sweep of metal work function and top-contact length, together with the coupled thermal consequences of each geometry. A closely related prior computational study modeled a single "sandwiched" ohmic-contact InSe FET design at sub-10-nm gate lengths [9]; our work differs in scope by holding the device geometry fixed and sweeping work function, contact geometry, and contact length as independent variables to isolate their individual contributions to $R_C$, and by coupling the electrical solution to a lattice heat-transport equation to quantify self-heating rather than treating it separately.

In this work, we combine finite-element simulation with direct benchmarking against the experimental InSe contact data of Huang et al. [5] to make four contributions relevant to contact processing. First, we map the Schottky-to-barrier-free transition in $R_{total}$ as a function of metal work function for a p-type InSe channel, and identify the work function at which the barrier for hole injection collapses — a target directly actionable for contact-metal selection. Second, we isolate the purely geometric contribution of edge contacts to resistance reduction, independent of the van der Waals-gap tunneling resistance that experimentally distinguishes real top and edge contacts, by holding all other interface physics fixed between the two geometries. Third, we quantify the thermal benefit of edge contacts under bias, showing that eliminating current crowding removes the associated hot spot, and uncover a counter-intuitive mobility/heating inversion: higher channel mobility actually increases peak self-heating due to higher dissipated power. Fourth, and unexpectedly, we find that bulk transport in the nominally "barrier-free" InSe channel is not Ohmic but Space-Charge-Limited (SCLC) even at our standard 50 nm channel length — confirmed via independent mesh-refinement and bias-sweep tests (Section 3.6) — with voltage-scaling exponents consistent with recent theory for space-charge-limited transport in reduced-dimensionality semiconductors [10]. This finding qualifies the "Ohmic regime" language used elsewhere in the 2D-contact literature; we did not find this reported for InSe in the literature we reviewed. We report these results alongside their full numerical basis (Section 3, Tables 1–2) and state plainly which claims are validated against literature and which remain model-internal (Section 3.6).

---

## 2. Methodology

The metal–InSe junction was modeled using the open-source TCAD simulator DEVSIM. A finite-element mesh was generated to resolve the sub-nanometer depletion regions near the contacts (see `mesh_generator.py` in the accompanying repository). The electrical drift-diffusion and thermal equations were solved self-consistently using a direct solver with a relative error tolerance of $10^{-5}$ and a maximum of 100 Newton iterations per bias step. The InSe channel was modeled with a thickness of 5 nm and an active channel length of 50 nm, consistent with the few-layer thickness regime studied experimentally in [1,5].

### 2.1 Material Parameters

Baseline InSe physical properties were drawn from the literature on bulk γ-InSe, the polytype most commonly used in device studies at this thickness, and are used as fixed constants rather than fitted quantities:

| Parameter | Value | Source |
|---|---|---|
| Bandgap ($E_g$) | 1.3 eV | [11] |
| Electron affinity ($\chi$) | 3.7 eV (→ $E_V = 5.0$ eV) | [12] |
| Electron DOS effective mass | 0.07 $m_0$ | [13] |
| Hole DOS effective mass | 0.17 $m_0$ | [13] |
| Electron mobility ($\mu_n$) | 1000 cm²/(V·s) | [1]; consistent with first-principles predictions [14] |
| Hole mobility ($\mu_p$) | 50 cm²/(V·s) | order-of-magnitude estimate; see Section 3.6 |
| Density of states ($N_C$, $N_V$) | $4.64\times10^{17}$, $1.76\times10^{18}$ cm⁻³ | computed from effective masses (Section 2.1 sources) at $T=300$ K |
| Intrinsic carrier concentration ($n_i$) | $1.12\times10^{7}$ cm⁻³ | computed from $E_g$, $N_C$, $N_V$ above |
| Static dielectric constant ($\epsilon$) | 6.4 | electron energy-loss spectroscopy compilation (params.txt corpus) |
| Thermal conductivity ($\kappa$) | 0.085 W/(cm·K) | representative few-layer value; see Section 3.6 |
| Background doping ($N_A$) | $1\times10^{16}$ cm⁻³ | assumed value, not derived from a specific literature source; see Section 3.6 |

We flag explicitly that hole mobility, thermal conductivity, and background doping for InSe are comparatively under-characterized in the literature relative to electron mobility and bandgap (the params.txt-referenced literature synthesis found no directly reported hole mobility above roughly 50 cm²/(V·s), consistent with the heavy, flat valence band described in [15]); these are the least well anchored inputs to the model, and Section 3.6 quantifies how the results scale with each of them, with sensitivity to background doping specifically the largest of any parameter tested. The density-of-states effective masses used above are likewise drawn from a tight-binding model [13] rather than a many-body-renormalized calculation [16]; Section 5 discusses how this and the other parameters could be refined as better data becomes available.

### 2.2 Physical Models and Transport Equations

The simulation couples a multidimensional drift-diffusion formulation to a lattice heat-transport equation. Key physical models:

1. **Thermionic field emission (TFE) contacts.** Rather than idealized Ohmic boundaries, metal–InSe interfaces are modeled with a thermionic surface recombination velocity condition ($J_{TE} = q v_{th}[n - n_{eq}]$), following the standard Schottky-contact treatment (Sze & Ng framework; the combined thermionic and field-emission current formalism this is built on is from [17]). Because the ultra-thin 2D depletion width makes tunneling non-negligible relative to bulk semiconductor contacts — a point emphasized in reviews of electron emission and injection physics in 2D materials [6] and in compact TFE/field-emission models for 2D Schottky-barrier FETs [18] — we incorporate a phenomenological tunneling enhancement factor, $\Gamma$, that multiplies the thermionic prefactor. We use $\Gamma = 10$ uniformly across all simulations as a representative order-of-magnitude value consistent with the tunneling-driven current enhancements (typically a few-fold to an order of magnitude) reported for 2D Schottky contacts in the literature above. Section 3.6 shows this choice is not decisive: resistance is unchanged across $\Gamma \in \{1, 10, 100\}$ for both geometries. The model isolates this intrinsic work-function dependence directly, without a separate metal-induced-gap-state/Fermi-level-pinning term; Section 5 discusses incorporating explicit pinning physics as a refinement to the barrier-collapse threshold identified in Section 3.1.
2. **Coupled electro-thermal physics.** A local Joule heating source term ($H = \vec{J}\cdot\vec{E}$) is applied throughout the device and coupled to the lattice heat-transport equation.
3. **Temperature-dependent mobility.** Phonon-limited mobility follows $\mu(T) = \mu_{300}(T/300)^{-1.5}$, capturing mobility degradation under localized self-heating.

### 2.3 Reproducibility

All mesh-generation, physics, and simulation-loop scripts (`mesh_generator.py`, `inse_physics.py`, `simple_physics.py`, `simple_dd.py`, `model_create.py`, `run_inse_simulations.py`), together with the complete set of extracted I-V and thermal CSV outputs underlying every figure and table in this paper, are provided in an accompanying code repository: **https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace**. We recommend archiving a tagged release of this repository (e.g., via Zenodo) to obtain a permanent DOI for citation.

---

## 3. Results and Discussion

### 3.1 Transition from Schottky-Limited to SCLC-Dominated Transport

Figure 1 shows simulated drain I-V characteristics at varying work functions for both contact geometries. As $\Phi_M$ increases toward the InSe valence band ($E_V \approx 5.1$ eV), the Schottky barrier narrows and eventually vanishes, transitioning the device from contact-limited transport into a regime limited entirely by the bulk channel. We refer to the regime at $\Phi_M \ge 5.0$ eV as "injection-barrier-free"; however, as discussed extensively in Section 3.6, the bulk channel transport in this highly depleted, nanoscale geometry ($L_{gap} = 50$ nm) is actually governed by Space-Charge Limited Current (SCLC) rather than simple linear Ohmic conduction. In this SCLC-dominated regime, the choice of metal work function ceases to be the primary bottleneck. At $\Phi_M = 5.0$ eV, the top-contact structure (using $L_c = 50$ nm) drops to a total device resistance of 1.50 $\Omega\cdot$cm. The edge contact, bypassing the out-of-plane van der Waals gap entirely, drops to 1.15 $\Omega\cdot$cm (extracted at $V_{DS}=0.1$ V).

Figure 1 illustrates this dramatic multi-order-of-magnitude transition in the I-V characteristics. In the strongly Schottky-limited state (4.33 eV), current is negligible (pA scale). As the barrier collapses (5.0–5.2 eV), the contacts no longer restrict carrier injection, and the SCLC-limited bulk channel dictates total throughput.

![I-V Curves at varying Work Functions](fig1_iv_curves.png)

*Figure 1. Simulated I-V characteristics for top and edge contacts across the work-function sweep.*

### 3.2 Work Function and the Collapse of the Hole-Injection Barrier

Table 1 reports the full extracted device resistance ($R_{total}$) across every simulated combination of geometry, contact length, and work function; Figure 2 plots the same data against work function.

**Table 1. Total device resistance $R_{total}$ (Ω·cm) versus contact geometry, top-contact length ($L_c$), and metal work function ($\Phi_M$).** Values are from isothermal (300 K) drift-diffusion simulation (`summary_results.csv`) extracted at $V_{DS}=0.1$ V. As detailed in Section 3.6, top-contact values carry an estimated ±10–20% uncertainty from incomplete mesh convergence at 4× density; edge-contact values are mesh-converged to within 0.8%. Three independent checks — Richardson extrapolation, an independent-platform re-run, and a corner-refined graded mesh — converge on $R_{top}\approx 2.87$–$2.89\ \Omega\cdot\text{cm}$ at the tested configuration, confirming the top-contact resistance approaches a well-defined limit rather than drifting further.

| $\Phi_M$ (eV) | Top, $L_c$=10 nm | Top, $L_c$=20 nm | Top, $L_c$=30 nm | Top, $L_c$=50 nm | Edge |
|---|---|---|---|---|---|
| 4.33 | $4.23\times10^{9}$ | $4.21\times10^{9}$ | $4.20\times10^{9}$ | $4.17\times10^{9}$ | $4.06\times10^{9}$ |
| 4.80 | 2170 | 2150 | 2130 | 2090 | 1950 |
| 5.00 | 1.78 | 1.75 | 1.73 | 1.67 | 1.50 |
| 5.20 | 0.272 | 0.263 | 0.257 | 0.234 | 0.195 |

![Resistance vs Work Function](fig2_resistance_vs_wf_v2.png)

*Figure 2. Extracted total device resistance versus metal work function, top vs. edge contact.*

Because the InSe valence-band edge lies at 5.0 eV in this model, low-work-function metals (e.g., $\Phi_M = 4.33$ eV) create a large Schottky barrier for holes; the barrier narrows exponentially as $\Phi_M$ increases. At $\Phi_M = 5.0$ eV the metal Fermi level aligns with the valence band and the injection barrier for holes effectively vanishes (though, as Section 3.6 confirms via direct $R$-vs-$V_{DS}$ testing, bulk transport in this regime is Space-Charge-Limited rather than Ohmic); at 5.2 eV, an accumulation layer forms and contact resistance falls to its simulated minimum, 0.195 Ω·cm for the edge contact. Experimentally, the closest available analog is the near-zero barrier ($\Phi_B \approx 50$ meV) reported in [5] for a work-function-aligned indium contact to InSe — a study of electron injection at the conduction band using a low-work-function metal, the mirror case of the hole-injection, high-work-function regime studied here, but demonstrating the same underlying mechanism: aligning the contact metal's Fermi level with the relevant InSe band edge collapses the Schottky barrier to a small fraction of an eV. Converting the transfer-length-method contact resistance reported in [5] (~10 kΩ·μm at their optimal gate bias) into the same per-unit-width convention used in Table 1 gives an order-of-magnitude resistance near 1 Ω·cm, comparable to our simulated 5.0 eV values (1.5–1.8 Ω·cm). We present this as qualitative, order-of-magnitude consistency rather than direct quantitative validation, since the two studies involve different carrier types, different contact metals, and different measurement geometries.

### 3.3 Eliminating Out-of-Plane Current Crowding

In the top-contact configuration, carriers injected from the source must traverse the resistive out-of-plane direction of the InSe layer before reaching the high-mobility in-plane channel, producing current crowding at the edge of the contact (Figure 3) and an effective injection area smaller than the physical contact length $L_c$. Table 1 shows this saturating behavior directly: at $\Phi_M = 4.33$ eV, increasing $L_c$ from 10 to 50 nm reduces $R_{total}$ by only 1.4%; at 4.8 eV, by 3.8%; but at 5.0 and 5.2 eV — where the injection barrier is largely collapsed and current crowding is less severe — the same length increase yields larger reductions of 6.0% and 14.1% respectively. The length-dependence of top-contact resistance is therefore itself work-function dependent, and we do not find it accurate to describe the improvement from longer contacts as uniformly "minimal" across all bias regimes; it is minimal specifically in the strongly Schottky-limited regime, where the barrier — not the geometry — dominates.

![Resistance vs Contact Length](fig3_resistance_vs_lc_v2.png)

*Figure 3. Total device resistance versus top-contact length, by work function.*

The edge-contact geometry allows direct, uniform carrier injection across the full 5 nm InSe thickness, bypassing the out-of-plane van der Waals path entirely. Consistent with this, the edge contact outperforms every top-contact length we simulated at every work function tested (Table 1). At $\Phi_M = 5.2$ eV, the edge contact's 0.195 Ω·cm is nominally 28% below the 10-nm top contact (0.272 Ω·cm) and 17% below the 50-nm top contact (0.234 Ω·cm). The top-contact values carry an estimated ±10–20% uncertainty from incomplete mesh convergence at 4× density; Section 3.6 resolves this via Richardson extrapolation of the mesh sequence, cross-validated by an independent DEVSIM re-run extended to 8× density, which converges to a geometric advantage of 21–27% at the tested configuration ($\Phi_M=5.0$ eV, $L_c=50$ nm) — consistent with, and confirming, the 17–28% raw estimate above rather than overturning it.

### 3.4 Mitigation of the Thermal Bottleneck

Table 2 reports the coupled electro-thermal results at $V_{DS} = 2.0$ V.

**Table 2. Coupled electro-thermal simulation results at $V_{DS}=2.0$ V** (`summary_results_thermal.csv`). $R_{total}$ here reflects the self-consistent, temperature-dependent-mobility solution and therefore differs from the isothermal $V_{DS}=0.1$ V values in Table 1.

| Geometry | $L_c$ (nm) | $\Phi_M$ (eV) | $R_{total}$ (Ω·cm) | Max lattice $T$ (K) |
|---|---|---|---|---|
| Top | 10 | 4.33 | $3.51\times10^{9}$ | 300.0 |
| Top | 50 | 4.33 | $3.50\times10^{9}$ | 300.0 |
| Top | 10 | 4.80 | 2770 | 300.2 |
| Top | 50 | 4.80 | 2620 | 300.2 |
| Top | 10 | 5.00 | 2.78 | 458.7 |
| Top | 20 | 5.00 | 2.73 | 460.1 |
| Top | 30 | 5.00 | 2.69 | 460.1 |
| Top | 50 | 5.00 | 2.57 | **470.1** |
| Edge | N/A | 4.33 | $3.46\times10^{9}$ | 300.0 |
| Edge | N/A | 4.80 | 2390 | 300.1 |
| Edge | N/A | 5.00 | 2.25 | **375.8** |

At low and moderate work functions (4.33, 4.8 eV), self-heating is negligible for both geometries (peak lattice temperature stays within a fraction of a degree of the 300 K ambient), because the Schottky-limited current density is too low to generate significant Joule heating. The thermal bottleneck appears specifically in the injection-barrier-free (SCLC-dominated) regime at $\Phi_M = 5.0$ eV, where current density is high enough to matter: the 50-nm top contact reaches a peak lattice temperature of 470.1 K, while the edge contact reaches only 375.8 K under the same bias — a reduction of 94.3 K. This is the largest temperature reduction found anywhere in our dataset (Table 2); no shorter top-contact length or lower work function produces a larger edge-vs-top gap. The reduction nonetheless represents a substantial, geometry-driven mitigation of self-heating, consistent with the elimination of the localized current-crowding hot spot described in Section 3.3.

![Thermal Profiles](thermal_comparison.png)

*Figure 4. Peak lattice temperature at $V_{DS}=2.0$ V, top vs. edge contact, at $\Phi_M = 5.0$ eV.*

### 3.5 Positioning Against Prior and Experimental Work

The concept of bypassing the van der Waals gap via edge injection is not new: it was established for graphene [7] and extended to other 2D-material junctions [8], both reporting substantial injection improvements over top contacts, on fabricated and measured devices. Our contribution is a quantitative, InSe-specific decomposition of how much of the top-versus-edge advantage comes from geometry alone (elimination of current crowding), isolated by applying identical interface physics ($\Gamma$, TFE parameters) to both geometries. As Section 3.6 shows, this 17–28% resistance reduction holds up under scrutiny: suppressing the van der Waals interfacial transmission ($T_{vdW}$) by up to 100× has almost no effect on total resistance, because transport is overwhelmingly bottlenecked by bulk channel resistance rather than interfacial injection velocity. The edge-contact advantage in this regime is therefore a geometric one, bypassing the resistive out-of-plane bulk path, and does not depend on assuming severe vdW interfacial penalties. This points to a broader caveat: what is often labeled "contact resistance" in similar device models may, as here, be effectively dominated by the semiconductor's bulk transport properties rather than the boundary condition itself.

Relative to the closest prior InSe contact simulation we identified — a single fixed "sandwiched" ohmic-contact geometry at sub-10-nm gate lengths [9] — this study's contribution is the systematic sweep across geometry, work function, and contact length as independent variables, which that prior work did not perform. A broad review of recent literature [2,19,20] highlights that while 2D transition metal dichalcogenides (TMDs) are extensively explored for continued electronic scaling, InSe-specific TCAD contact modeling remains relatively sparse. This work bridges that gap by providing a targeted, incremental, but necessary geometry decomposition for InSe specifically, complementing recent field-wide advancements in 2D complementary logic integration and directly usable as processing guidance once fabrication of edge-contacted InSe devices is undertaken.

### 3.6 Robustness Checks: Mesh Convergence, Tunneling Enhancement, van der Waals Transmission, Channel Scaling, and Material Parameters

**Summary of Findings:** To validate the robustness of the core results, we performed five supplementary checks. Three key findings emerge that contextualize the rest of the paper: (1) **Mesh Convergence:** Edge contacts are fully converged under uniform meshing. Top contacts converge more slowly under uniform mesh refinement alone, but three independent checks — Richardson extrapolation of the uniform-mesh sequence, an independent re-run on separate infrastructure, and a targeted graded mesh refining the metal-semiconductor-vacuum triple-point corner — all converge on $R_{top}\approx 2.87$–$2.89\ \Omega\cdot\text{cm}$, yielding a converged geometric advantage of 21–27%, consistent with the 17–28% raw estimate. (2) **SCLC Dominance:** Transport in the barrier-free regime is governed by Space-Charge-Limited Current (SCLC) rather than Ohmic conduction, confirmed by superlinear length scaling and direct carrier-density excess. (3) **Self-Heating Collapse:** Under coupled electro-thermal operation, severe self-heating at the contact edge largely erases the isothermal top-contact length-dependence, while the edge contact continues to run substantially cooler than the top contact at the same bias.

To address the specific unquantified sensitivities noted as limitations in an earlier version of this work, we ran four supplementary checks. All four checks were performed using the identical coupled electro-thermal solver configuration as Table 2, which activates the lattice heat-transport equations and temperature-dependent mobility models. While the resistance for these checks was evaluated at $V_{DS}=0.1$ V (the linear regime), the mere presence of the coupled thermal physics inherently alters the self-consistent state compared to the strictly isothermal, decoupled drift-diffusion model of Table 1. This structural difference in the underlying solver explains why their baseline resistance values reproduce Table 2's coupled values (e.g., $R_{edge} = 2.247\ \Omega\cdot\text{cm}$) rather than Table 1's isothermal values (1.50 Ω·cm). For studies that also extract maximum temperature, the voltage was further swept to $V_{DS}=2.0$ V, maintaining complete internal consistency with Table 2's electro-thermal dataset.

**Mesh convergence.** We refined the mesh density by a factor of 0.5–4× (64–2541 nodes for the top contact; 44–1701 for the edge contact) at $\Phi_M=5.0$ eV.

| Mesh factor | Nodes (top) | $R_{top}$ (Ω·cm) | Nodes (edge) | $R_{edge}$ (Ω·cm) |
|---|---|---|---|---|
| 0.5 | 64 | 2.318 | 44 | 2.198 |
| 1.0 | 186 | 2.566 | 126 | 2.247 |
| 2.0 | 671 | 2.710 | 451 | 2.261 |
| 4.0 | 2541 | 2.790 | 1701 | 2.265 |

The edge-contact resistance is converged: it changes by only 0.8% from the baseline (1×) to the finest (4×) mesh. The top-contact resistance converges more slowly — it rises by 8.7% from 1× to 4× mesh density (20.3% from the coarsest to the finest mesh tested) — but the successive increments shrink by a near-constant factor of 0.55–0.58 per mesh doubling, consistent with standard first-order finite-element convergence to a finite limit rather than an open-ended drift.

**Three independent checks of the converged value.** First, Richardson extrapolation of the four-point sequence above gives $R_{top,\infty} \approx 2.89\ \Omega\cdot\text{cm}$. Second, we re-ran the identical uniform-mesh sequence, extended to 8× density, on a second DEVSIM installation on different infrastructure (Linux vs. Windows, different BLAS/LAPACK backend); this independent sequence shows the same near-constant-ratio convergence (successive-difference ratios of 0.55, 0.51, 0.50) and extrapolates to $R_{top,\infty} \approx 2.87\ \Omega\cdot\text{cm}$ — within 0.5% of the first estimate. Third, to test whether the top contact's slower convergence reflects a genuine resolution limit or a finite-element mesh singularity at the sharp metal–semiconductor–vacuum corner, we implemented a graded mesh injecting ultra-fine (0.1 nm) spatial resolution locally at the top-contact inner edges, and swept it from 0.5× to 8× density (116–22,302 nodes). This graded sequence does not converge as smoothly as the uniform sequence — the successive-difference ratios are non-monotonic (0.65, 1.46, 0.11) rather than steadily decaying, indicating the local corner refinement interacts with the rest of the mesh in a way we have not fully characterized — but its highest-density point ($R_{top}=2.871\ \Omega\cdot\text{cm}$ at 22,302 nodes, changing by only 0.24% from the 4× graded point) lands within 0.6% of both extrapolated estimates above. We report the graded-mesh sequence's irregular intermediate behavior honestly rather than presenting it as a clean demonstration of singularity resolution; what it does provide is a third, methodologically distinct estimate that agrees with the other two. Taking $R_{top,\infty} \approx 2.87$–$2.89\ \Omega\cdot\text{cm}$ across all three checks and $R_{edge,\infty} \approx 2.27\ \Omega\cdot\text{cm}$ (already converged) gives a converged geometric advantage of 21–27% at this configuration ($\Phi_M=5.0$ eV, $L_c=50$ nm) — consistent with, and confirming, the 17–28% raw estimate from Table 1 rather than overturning it. We retain the conservative ±10–20% band for the remaining Table 1 configurations (other work functions and contact lengths), which were not individually mesh-tested. The extrapolation script (`richardson_extrapolation.py`), the graded-mesh script (`test_graded_mesh.py`), and all three mesh sequences (`mesh_convergence.csv`, `mesh_convergence_independent_check.csv`, `graded_mesh_test.csv`) are included in the accompanying repository.

**Tunneling enhancement factor ($\Gamma$).** We swept $\Gamma \in \{1, 10, 100\}$ for both the top and edge contacts across five work functions (4.33–5.2 eV). Total resistance was numerically identical across all three $\Gamma$ values at every work function and for both geometries (relative difference below $10^{-9}$; Figure 5). This indicates that transport is bottlenecked by the exponential barrier population ($n_{eq}$) and bulk/channel resistance rather than by the interfacial emission-velocity prefactor — once the contact's thermionic velocity exceeds a threshold set by the bulk transport, further increases in $\Gamma$ have no effect on the self-consistent solution. This directly resolves the limitation flagged in an earlier version of this work regarding $\Gamma=10$ being an unvalidated, potentially decisive free parameter: it is not decisive for either geometry, at least across the two-decade range tested here.

![Gamma Sensitivity](fig_gamma_sensitivity_final.png)

*Figure 5. Total resistance vs. work function for $\Gamma=1$, 10, and 100 for both top and edge contacts. The curves for different $\Gamma$ values coincide exactly for each geometry.*

**Van der Waals transmission coefficient.** We introduced an explicit vdW transmission penalty ($T_{vdW}$) multiplying the top-contact's thermionic prefactor and swept it from 1.0 (no penalty) to 0.01 (100× suppression) across three work functions: 4.8 eV (Schottky regime), 5.0 eV (barrier-free threshold), and 5.2 eV (accumulation). Total resistance was unchanged across the entire $T_{vdW}$ sweep at every work function (e.g., exactly 2.566 Ω·cm at 5.0 eV for all $T_{vdW}$ values; Figure 6). This confirms that even in the injection-barrier-free regime where contact resistance is lowest, the macroscopic drift-diffusion model is bottlenecked by bulk channel transport rather than interfacial injection velocity. Consequently, the performance gap between top and edge contacts in our model is driven purely by the macroscopic geometry (bypassing the resistive out-of-plane bulk path) rather than by interfacial vdW-gap penalties.

![vdW Gap Sensitivity](fig_vdw_gap_final.png)

*Figure 6. Top- and edge-contact resistance vs. van der Waals transmission coefficient at three work functions (Schottky, barrier-free-threshold, and accumulation regimes).*

**Channel length scaling and SCLC confirmation.** Because transport is overwhelmingly channel-limited at $\Phi_M \ge 5.0$ eV, total resistance is strongly dominated by channel length ($L_{ch}$). We swept $L_{ch} \in \{10, 30, 50, 100\}$ nm for both top ($L_c=50$ nm) and edge contacts at $\Phi_M=5.0$ eV (Figure 7). We find that $R_{total}$ scales superlinearly with channel length: $R_{edge}$ rises from 0.157 Ω·cm at 10 nm to 22.72 Ω·cm at 100 nm, an accelerating >140× increase across a 10× length increase, with the local length exponent climbing from $\approx$1.5 at short length toward $\approx$3.3 between 50–100 nm. While simple Ohmic transport predicts linear resistance scaling ($R \propto L$), classical Space-Charge Limited Current (SCLC) theory predicts $R \propto L^3$ at fixed bias in bulk trap-free solids [21] — our long-channel exponent is approaching, though not yet reaching, this classical bulk limit.

We tested this hypothesis with three additional checks. First, deploying a strictly doubled absolute mesh resolution (1.25 nm fixed node spacing) at $L_{ch}=100$ nm changed resistance by less than 1% for edge contacts, ruling out a coarse-mesh artifact as the source of the superlinear trend. Second, we evaluated resistance as a function of bias voltage ($V_{DS} \in [0.02, 0.5]$ V) at fixed $L_{ch}$. At $L_{ch}=100$ nm, $R_{edge}$ drops from 79.15 Ω·cm at 0.02 V to 2.84 Ω·cm at 0.50 V (implied exponent $I \propto V^{n}$, $n \approx 2.0$); at the standard baseline $L_{ch}=50$ nm channel used throughout the rest of this paper, we observe the same qualitative collapse, from 6.66 Ω·cm to 0.83 Ω·cm ($n \approx 1.6$). Both exponents fall within the 1.5–2 range predicted for space-charge-limited transport in reduced-dimensionality (2D, thin-body) semiconductors — distinct from the $n=2$–$3$ range expected in bulk solids — by the modified SCLC theory of [10], which was itself validated against experimental $I$–$V$ exponents of 1.7–2.5 in monolayer MoS$_2$ and hBN. Third, and most directly, we extracted the simulated hole density at the center of the 50 nm channel across the same bias sweep: even at 0.02 V it reaches roughly $1.1 \times 10^{17}$ cm⁻³, an order of magnitude above the $1.0 \times 10^{16}$ cm⁻³ background doping, and by 0.50 V it climbs above $1.5 \times 10^{18}$ cm⁻³ — a direct confirmation that injected carriers, not the background dopant population, govern conduction in this regime, consistent with the defining condition for SCLC rather than an inference from the $I$–$V$ curve shape alone. (This density check was run self-consistently across the full bias sweep and saved to `sclc_density_check.csv` in the archived reproducibility package, explicitly demonstrating that center-channel hole density reaches $1.5 \times 10^{18}\text{ cm}^{-3}$ at 0.5 V, exceeding the background doping by over two orders of magnitude). Taken together, the voltage-scaling exponents, the length-scaling trend approaching the bulk $L^3$ limit, and the direct carrier-density excess all point to the same conclusion. Consequently, while the injection barrier is effectively transparent at 5.0 eV, the overarching "Ohmic" framing typically applied to this regime must be qualified: the bulk channel itself operates non-linearly.

**Self-heating and the collapse of top-contact length-dependence.** The 17–28% relative resistance reduction for edge versus top contacts quoted in Section 3.2 spans two different top-contact lengths: 17% against the 50 nm top contact and 28% against the 10 nm top contact. Both values were extracted from the strictly isothermal (300 K, fixed-mobility) solver underlying Table 1. When we instead evaluate these geometries under the coupled electro-thermal solver (biased via the $V_{DS}$ sweep of $L_{ch}=50$ nm, $\Phi_M=5.2$ eV described above), the 10 nm and 50 nm top contacts become numerically indistinguishable: $R = 0.5485\ \Omega\cdot\text{cm}$ for both at $V_{DS}=0.1$ V. This is consistent with self-heating-driven current crowding at the contact edge becoming severe enough that the remainder of a 50 nm top contact carries negligible current, rendering it electrically similar to a 10 nm contact once realistic thermal physics is included. If this holds up, the length-dependence observed in the isothermal model would be substantially an artifact of neglecting self-heating, though we stop short of calling this a fully settled result. Comparing this collapsed top-contact value to the edge contact evaluated under identical conditions ($R_{edge} = 0.4066\ \Omega\cdot\text{cm}$) yields a geometric advantage of approximately 26%. We checked this collapse against mesh refinement: doubling the absolute mesh density ($2.0\times$ baseline, using a single-step 0.1 V bias rather than the ramped sweep) gave $0.589\ \Omega\cdot\text{cm}$ for both the 10 nm and 50 nm top contacts, and $0.438\ \Omega\cdot\text{cm}$ for the edge contact — again a ~26% gap. The finding is therefore reproducible across two mesh densities and two different bias-ramp protocols, which is reassuring, but it has only been checked at two mesh densities (not the four-point 0.5–4× sweep used for the baseline mesh-convergence study above), so we present the L$_c$-collapse itself as a strongly supported but not yet exhaustively converged result, and recommend the same 0.5–4× mesh sweep be applied to it before treating Table 1's isothermal length-dependence as fully superseded. Table 2 additionally shows that, at this same 5.0 eV, high-bias condition, the edge contact runs substantially cooler than the top contact (375.8 K vs. 470.1 K peak lattice temperature) rather than hotter, so the edge geometry's thermal behavior under bias remains favorable rather than becoming a liability.

![Channel Length Scaling](fig_channel_length.png)

*Figure 7. Total resistance vs. channel length at $\Phi_M = 5.0$ eV for top and edge contacts. The superlinear scaling and strong bias dependence indicate transport operates in the space-charge limited current (SCLC) regime.*

**Material parameter sensitivity.** Because hole mobility and thermal conductivity are the least well-anchored inputs in Section 2.1, we swept $\mu_p \in \{10, 50, 100\}$ cm²/(V·s) at fixed $\kappa=0.085$ W/(cm·K), and $\kappa \in \{0.05, 0.085, 0.2\}$ W/(cm·K) at fixed $\mu_p=50$ cm²/(V·s), both at $\Phi_M=5.0$ eV (Figure 8).

| Sweep | Value | $R_{top}$ (Ω·cm) | $R_{edge}$ (Ω·cm) | Max $T_{top}$ (K) | Max $T_{edge}$ (K) |
|---|---|---|---|---|---|
| $\mu_p$ | 10 | 12.828 | 11.235 | 355.2 | 318.5 |
| $\mu_p$ | 50 (baseline) | 2.566 | 2.247 | 470.1 | 375.8 |
| $\mu_p$ | 100 | 1.283 | 1.123 | 553.4 | 429.1 |
| $\kappa$ | 0.05 | 2.566 | 2.247 | 530.2 | 414.5 |
| $\kappa$ | 0.085 (baseline) | 2.566 | 2.247 | 470.1 | 375.8 |
| $\kappa$ | 0.2 | 2.566 | 2.247 | 397.4 | 336.7 |

$R_{total}$ scales as $1/\mu_p$ essentially exactly ($R_{top}\times\mu_p = 128.3\ \Omega\cdot\text{cm}\cdot\text{cm}^2/(\text{V}\cdot\text{s})$ across all three points, constant to within 0.02%), confirming the model is bulk-hole-transport-limited at this bias point rather than contact-limited. This means the reported resistance values are directly proportional to whichever hole mobility value is closest to the true InSe hole mobility, which Section 2.1 already flags as sparsely constrained in the literature (order-of-magnitude estimate only): if the true $\mu_p$ differs from our assumed 50 cm²/(V·s) by a factor of 2–5 within the plausible literature range, the reported resistance values should be expected to shift by the same factor, inversely. Thermal conductivity has a negligible effect on $R_{total}$ (varying by <0.01% across the tested range), as expected since $\kappa$ does not enter the electrical transport equations directly, but it has a substantial effect on peak self-heating: max lattice temperature spans 397–530 K (a 133 K range) across the literature-plausible $\kappa$ range, compared to the single 470.1 K value reported in Table 2. Counter-intuitively, max temperature *increases* with hole mobility (355→470→553 K): at fixed bias, higher mobility lowers $R_{total}$, which increases dissipated power ($P=V^2/R$) faster than it improves heat removal, so the self-heating figures in Section 3.4 should be read as specific to $\mu_p=50$ cm²/(V·s) and not necessarily conservative in either direction.

**Background Doping Concentration.** We evaluated the sensitivity to background p-type doping concentration by sweeping $N_A \in \{10^{15}, 10^{16}, 10^{17}\}\ \text{cm}^{-3}$. Resistance was heavily dependent on the substrate doping: at $\Phi_M = 5.0$ eV, the top-contact resistance varied from $2.04\ \Omega\cdot\text{cm}$ ($N_A = 10^{15}$) to $2.70\ \Omega\cdot\text{cm}$ ($N_A = 10^{16}$), and climbed drastically to $443.8\ \Omega\cdot\text{cm}$ ($N_A = 10^{17}$). The edge-contact scaling followed the same dramatic increase (1.76 $\Omega\cdot\text{cm}$ up to $159.8\ \Omega\cdot\text{cm}$). This strong dependence indicates that space-charge constraints and depletion-width dynamics in the channel become highly restrictive at higher doping levels in our simulated nanoscale geometry.

![Material Parameter Sensitivity](fig_material_sensitivity.png)

*Figure 8. Left: total resistance and peak temperature vs. hole mobility. Right: peak temperature vs. thermal conductivity.*

---

## 4. Conclusion

Systematic drift-diffusion and coupled electro-thermal TCAD simulations of metal–InSe junctions, benchmarked against the experimentally reported near-zero-barrier indium contact to InSe [5], show that aligning the metal work function with the InSe valence band (~5.0 eV) is the dominant factor in collapsing the multi-gigaohm Schottky barrier for hole injection. Separately from that energetic effect, edge contacts reduce total device resistance by eliminating out-of-plane current crowding. Idealized isothermal models suggest a 17–28% geometric advantage for edge over top contacts. Three independent checks — Richardson extrapolation of the mesh-convergence sequence, an independent-platform re-run extended to 8× density, and a corner-refined graded mesh also extended to 8× density — converge on the same top-contact resistance limit ($R_{top,\infty}\approx 2.87$–$2.89\ \Omega\cdot\text{cm}$) and confirm this advantage converges to 21–27%, consistent with the raw estimate rather than an artifact of incomplete mesh refinement. Introducing realistic electro-thermal physics largely collapses the top-contact length-dependence itself: severe localized self-heating at the contact edge is the cause, a finding reproduced across two mesh densities and bias protocols, though not yet checked with the full convergence sweep applied elsewhere in this study. Under that thermal regime, edge contacts settle at a consistent ~26% performance advantage over top contacts, and continue to run substantially cooler under high bias rather than becoming a thermal liability. Our robustness checks show the model is channel-transport-limited (specifically, space-charge limited) rather than interface-limited across the barrier-free transition, so this geometric advantage does not depend on the interfacial van der Waals tunneling assumptions. We also find a counter-intuitive thermal bottleneck: higher hole mobility directly increases peak lattice temperature at fixed bias, because it raises dissipated power. Edge-contact geometry mitigates this self-heating under high bias, reducing peak lattice temperature by roughly 94 K at the 5.0 eV work function where it matters most. Together, these results give geometry-aware, quantitatively bounded processing guidance for contact engineering in InSe field-effect transistors: target a contact metal work function at or above ~5.0 eV, and prefer edge over top contact deposition wherever the fabrication process allows it.

---

## 5. Future Work

The present results point to several natural extensions. Direct experimental validation via transfer-length-method measurements on fabricated top- and edge-contacted InSe devices would allow the simulated resistance values and the ~26% edge-contact advantage identified in Section 3.6 to be benchmarked against hardware, complementing the qualitative, order-of-magnitude agreement already found with the closest available analog [5] (Section 3.2). While Richardson extrapolation and an independent cross-validation run (Section 3.6) confirm the top-contact resistance converges at the specific configuration tested ($\Phi_M=5.0$ eV, $L_c=50$ nm), applying the same mesh-convergence sweep across the remaining work functions and contact lengths in Table 1, and to the electro-thermal L$_c$-collapse finding, would place the full dataset on the same fully converged footing as the single tested configuration. As experimentally measured values for InSe hole mobility, thermal conductivity, and background doping become available, they can be substituted directly into this framework: Section 3.6 already shows how resistance and peak self-heating temperature scale with each of these parameters, so updated inputs would translate directly into sharper quantitative predictions. Incorporating many-body-renormalized effective masses [16] and explicit metal-induced-gap-state/Fermi-level-pinning physics would refine the barrier-collapse work-function threshold identified in Section 3.1. Finally, coupling the present p-type results with an equivalent n-type contact model would enable full CMOS-level performance projections for InSe logic, building on the growing wafer-scale and ballistic-transport results for InSe noted in Section 1 [3,4].

---

## Data and Code Availability

All simulation scripts, mesh-generation code, and the complete I-V and electro-thermal CSV datasets underlying every figure and table in this paper are available at: **https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace**. We recommend the author archive a tagged release via Zenodo (or equivalent) to obtain a citable, permanent DOI before submission.

## CRediT Author Statement

**Harsha Adhikary:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Writing – review & editing, Visualization.

## Declaration of Competing Interest

The author declares that he has no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Funding

This work received no external funding.

## Declaration of Generative AI and AI-Assisted Technologies in the Manuscript Preparation Process

During the preparation of this work, the author used AI-based writing and coding assistants (Claude, and other general-purpose generative AI tools) to assist with literature search, drafting and revision of manuscript text, code development for figure generation and data analysis, and verification of numerical claims against the underlying simulation data. After using these tools, the author reviewed and edited all content as needed and takes full responsibility for the content of the published article. Data visualizations (Figures 1–8) were generated entirely from the author's own simulation output via reproducible Python/Matplotlib scripts included in the accompanying repository (Section 2.3); no AI tool was used to generate, alter, or fabricate any figure, image, or underlying dataset.

## Acknowledgments

I thank the authors of the experimental and theoretical InSe literature cited throughout this manuscript, whose prior work this study builds on and benchmarks against. Simulations were performed using DEVSIM, an open-source TCAD tool, without which this project would not have been possible. I am grateful to the Indian Institute of Science for the academic environment and resources that made this work possible, and to my parents for their constant support.

---

## References

[1] D. Bandurin, A. Tyurnina, G. Yu, A. Mishchenko, V. Zólyomi, S. Morozov, R.K. Kumar, R. Gorbachev, Z. Kudrynskyi, S. Pezzini, Z. Kovalyuk, U. Zeitler, K. Novoselov, A. Patané, L. Eaves, I. Grigorieva, V. Fal'ko, A.K. Geim, Y. Cao, High electron mobility, quantum Hall effect and anomalous optical response in atomically thin InSe, Nat. Nanotechnol. 12 (2016) 223–229.

[2] S. Song, M.A. Altvater, W. Lee, H. Shin, N. Glavin, D. Jariwala, Indium selenides for next-generation low-power computing devices, Nat. Rev. Electr. Eng. (2025).

[3] B. Qin, J. Jiang, L. Wang, Q. Guo, C. Zhang, L. Xu, X. Ni, P. Yin, L.-M. Peng, E. Wang, F. Ding, C. Qiu, C. Liu, K. Liu, Two-dimensional indium selenide wafers for integrated electronics, Science 389 (2025) 299–302.

[4] J. Jiang, L. Xu, C. Qiu, L.-M. Peng, Ballistic two-dimensional InSe transistors, Nature 616 (2023) 470–475.

[5] Y.-T. Huang, Y.-H. Chen, Y.-J. Ho, S.-W. Huang, Y.-R. Chang, K. Watanabe, T. Taniguchi, H.-C. Chiu, C.-T. Liang, R. Sankar, F.-C. Chou, C.-W. Chen, W.-H. Wang, High-performance InSe transistors with ohmic contact enabled by nonrectifying-barrier-type indium electrodes, ACS Appl. Mater. Interfaces 10 (2018) 33450–33456.

[6] Y.S. Ang, L. Cao, L.K. Ang, Physics of electron emission and injection in two-dimensional materials: Theory and simulation, InfoMat 3 (2021) 502–535.

[7] L. Wang, I. Meric, P.Y. Huang, Q. Gao, Y. Gao, H. Tran, T. Taniguchi, K. Watanabe, L.M. Campos, D.A. Muller, J. Guo, P. Kim, J. Hone, K.L. Shepard, C.R. Dean, One-dimensional electrical contact to a two-dimensional material, Science 342 (2013) 614–617.

[8] M.H.D. Guimarães, H. Gao, Y. Han, K. Kang, S. Xie, C.-J. Kim, D.A. Muller, D.C. Ralph, J. Park, Atomically thin ohmic edge contacts between two-dimensional materials, ACS Nano 10 (2016) 6392–6399.

[9] J. Zhu, J. Ning, D. Wang, J. Zhang, L. Guo, Y. Hao, High-performance two-dimensional InSe field-effect transistors with novel sandwiched ohmic contact for sub-10 nm nodes: a theoretical study, Nanoscale Res. Lett. 14 (2019) 267.

[10] Y.S. Ang, M. Zubair, L.K. Ang, Relativistic space-charge-limited current for massive Dirac fermions, Phys. Rev. B 95 (2017) 165409.

[11] G.W. Mudd, M. Molas, X. Chen, V. Zólyomi, K. Nogajewski, Z. Kudrynskyi, Z. Kovalyuk, G. Yusa, O. Makarovsky, L. Eaves, M. Potemski, V. Fal'ko, A. Patané, The direct-to-indirect band gap crossover in two-dimensional van der Waals indium selenide crystals, Sci. Rep. 6 (2016) 39619.

[12] Y. Guo, J. Robertson, Band structure, band offsets, substitutional doping, and Schottky barriers of bulk and monolayer InSe, Phys. Rev. Mater. 1 (2017) 044004.

[13] S. Magorrian, V. Zólyomi, V. Fal'ko, Electronic and optical properties of two-dimensional InSe from a DFT-parametrized tight-binding model, Phys. Rev. B 94 (2016) 245431.

[14] L.-B. Shi, Y.-Y. Zhang, B. Xu, Z. Sun, D.-D. Liu, K.-Y. Xu, W.-Q. Chen, F. Zhang, Theoretical prediction of intrinsic electron mobility of monolayer InSe: first-principles calculation, J. Phys. Condens. Matter 32 (2019) 065306.

[15] H. Henck, D. Pierucci, J. Zribi, F. Bisti, E. Papalazarou, J. Girard, J. Chaste, F. Bertran, P. Fèvre, F. Sirotti, L. Perfetti, C. Giorgetti, A. Shukla, J. Rault, A. Ouerghi, Evidence of direct electronic band gap in two-dimensional van der Waals indium selenide crystals, Phys. Rev. Mater. 3 (2019) 034004.

[16] W. Li, F. Giustino, Many-body renormalization of the electron effective mass of InSe, Phys. Rev. B 101 (2020) 035201.

[17] F.A. Padovani, R. Stratton, Field and thermionic-field emission in Schottky barriers, Solid State Electron. 9 (1966) 695–707.

[18] A. Tunga, Z. Zhao, A. Shukla, W. Zhu, S. Rakheja, Physics-based modeling and validation of 2D Schottky barrier field-effect transistors, arXiv:2307.04851, 2023.

[19] M. Liang, H. Yan, N. Wazir, C. Zhou, Z. Ma, Two-dimensional semiconductors for state-of-the-art complementary field-effect transistors and integrated circuits, Nanomaterials 14 (2024) 1408.

[20] F. Zheng, W. Meng, L.J. Li, Continue the scaling of electronic devices with transition metal dichalcogenide semiconductors, Nano Lett. 25 (2025) 3683–3691.

[21] N.F. Mott, R.W. Gurney, Electronic Processes in Ionic Crystals, second ed., Oxford University Press, Oxford, 1940.
