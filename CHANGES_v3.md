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
