# Changes from v1 to v2

## Defects in v1 (all fixed in `src_v2/inse2d.py`)
| # | v1 defect | Effect | Fix |
|---|---|---|---|
| 1 | Contact flux `q v (1+Gamma)(p - p_eq)` was the *residual* of the DEVSIM contact equation, which replaces the bulk equation, so `p_s = p_eq` exactly (Dirichlet) | Gamma and T_vdW had no effect (identical to 1e-9) | Robin flux BC balanced against the drift-diffusion current (`edge_model` + `node_model * ContactSurfaceArea/NodeVolume`); verified (V2, V4) |
| 2 | `CreatePE` never called: Poisson kept the equilibrium charge after drift-diffusion | No space-charge feedback, so SCLC could not be inferred | Full charge q(p - n + N_D - N_A); verified against Mott-Gurney (V3) |
| 3 | `NetDoping = +1e16` (n-type) | Contradicted the p-type text; wrong channel | `NetDoping = -N_A` |
| 4 | Contact potential used `-0.5 V_t ln(Nc/Nv)` | 34 meV inconsistency between potential BC and p_eq | `+0.5 V_t ln(Nc/Nv)`; p_eq derived from psi_c |
| 5 | R from 0.1 V paired with T from 2 V | Mixed operating points | All quantities at the same bias; T vs power |
| 6 | `mesh_generator.py` truncated inside the edge mesh | Edge case not reproducible | New meshes in `inse2d.py` |
| 7 | `summary_results.csv` not produced by any public script; differed from the thermal solver by 1.5-2x at 0.1 V | Table 1 not reproducible | One solver; isothermal/thermal agree to <2e-4 |
| 8 | Joule heat J.E; fixed-T substrate directly under the film | Heat source can be negative; unrealistic heat sink | Quasi-Fermi Joule heat (Wachutka); Robin h_sub; energy balance checked |
| 9 | p_eq up to 2000 N_V for Phi_M > E_V | Outside Boltzmann validity | Negative barriers clipped and flagged |
| 10 | "21-27%" advantage | Arithmetic gives 20.9-21.5% | Superseded |

## Conclusions withdrawn
* Resistance is invariant to the tunnelling factor (artifact of #1).
* The barrier-free channel is SCLC (no Poisson coupling, #2; the corrected model shows m = 1.3-1.5 < 2).
* The edge advantage is purely geometric (corrected model: it is mostly the vdW gap).
* Edge contacts run 94 K cooler (artifact of #5 and #8; corrected model: edge runs hotter at equal power).
* The doping sensitivity of 443 Ohm cm at 1e17 (wrong-sign doping, #3).

## New results (v2)
See `paper_v2/InSe_revised.pdf`, Sec. 3.
