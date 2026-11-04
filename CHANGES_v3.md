# v3 additions (relative to v2)

| Addition | Where | Key result |
|---|---|---|
| Back gate (oxide + gate electrode, flat band at V_G=0) | `tcad2d.py` `_bottom` | SB-pFET transfer curves, on/off > 1e12 |
| Fermi-Dirac statistics (Nilsson eta, Bednarczyk F1/2) | `delta_expr`, `F12_expr` | Contacts beyond E_V allowed; R falls a further 2-2.6x up to 5.65 eV |
| Anisotropic mobility (mu_perp/mu_par) and kappa (8.5 / 0.76 W/mK) | `_build`, `_build_thermal` | Top-contact R doubles for mu ratio 1 -> 0.002; edge unaffected |
| Nonlocal Tsu-Esaki/WKB tunnelling (semi-implicit G = a - b n) | `update_tunnel`, `_tsu_esaki` | 80-100 % of gated on-current; local Gamma underestimates ~60x |
| T-dependent statistics/injection | `nonisothermal` | <1 % effect on temperatures |
| Interface, Peltier, tunnelling heat; energy closure | `_build_thermal`, `junction_heat` | heat out = IV within 0.3-1.2 %; v2 bulk-only misses up to 23 % |
| Transient heating | `transient_step` | tau63 ~ 0.20-0.25 ns |
| Layer-dependent band gap / affinity | `band()` | top-contact optimum ~4 layers when gated |
| Separate S_top / S_edge | `phi_eff` | S_edge = 0.7 already reverses the top/edge ranking |
| MoS2 calibration (English et al. 2016) | `studies.py s_mos2` | T_vdW = 1.0; L_T(300 K) 34 nm vs 35-40 nm measured |
| Compact TLM model | `make_figures.py` | median error 11 % vs TCAD |
| Monte Carlo UQ (120 samples x 4) | `studies.py s_uq` | pinning and E_V dominate; edge < top in 100 % of paired samples |
| Verification V6-V9, CI | `verification.py`, `.github/workflows/verify.yml` | all pass |

Fixes made during v3 development: FD-consistent equilibrium initialisation; carriers solved before
the heat equation is switched on; mobility temperature dependence in isothermal runs at T0 != 300 K.

## Figure review (v3.1)

All ten v3 figures were checked for axis labels and units, legend overlap, colour use, grayscale legibility and print-size text.

Fixed:
- **Fig. 7(b)** was empty: the plotting code filtered on 2 V, but the thermal ablation was run at 0.5 and 1 V. It now shows the 1 V ablation with readable labels. The panel (c) title said "0→2 V"; the transient is a 0→1 V step (tunnelling off).
- **Fig. 8**: panel (b) mixed nm and 10^-7 Ω cm^2 on one unlabelled axis; it is now two panels with their own axes and units. The annotation that overlapped data became a legend entry.
- **Fig. 10**: colours now follow the rest of the paper (blue = edge, orange = top; hatching = 5.0 eV). The overlapping histograms are replaced by cumulative distributions, and a new panel shows the paired top/edge samples (edge < top in 100% of pairs at both work functions).
- **Fig. 4**: legends no longer cover data; panel (b) no longer uses the edge colour for a top contact; panel (c) uses geometry colours with hatching for κ⊥ and gives units in W/mK.
- **Fig. 6**: the + and Y-shaped symbols in panel (a) are now explained in the legend.
- **Fig. 2, Fig. 1, Fig. 5**: legend overlap (Fig. 2b, Fig. 1a), integer layer ticks, an unlabelled bulk-gap line and the analytic line (Fig. 1b), and line styles so the pinning curves are distinguishable without colour (Fig. 5).

Added:
- **Fig. 0** (`make_schematic.py`): device geometries.
- **Fig. 11** (`make_maps.py`): 2D hole-density and temperature maps of the gated top and edge devices at 1 V, with fields in `data_v3/fields_*.csv`.

No simulation results changed.
