# Geometry-Dependent Contact Resistance and Space-Charge-Limited Bulk Transport in Metal–InSe Junctions

Simulation code, data, and manuscript accompanying the paper *"Geometry-Dependent Contact
Resistance and Space-Charge-Limited Bulk Transport in Metal–InSe Junctions: A Coupled
Electro-Thermal Finite-Element Simulation Study"* (Harsha Adhikary, Department of Physics,
Indian Institute of Science).

Coupled drift-diffusion / lattice-heat-transport finite-element simulations (built on the
open-source TCAD tool [DEVSIM](https://devsim.org/)) of metal–InSe junctions, comparing
top-contact and edge-contact geometries across a sweep of metal work function and contact
length, and quantifying self-heating under bias.

## Repository structure

```
src/       DEVSIM mesh-generation, physics, simulation-loop, and plotting scripts
data/      CSV outputs (I-V curves, summary tables, robustness-check sweeps)
figures/   The 8 figures used in the manuscript
paper/     Final manuscript (Markdown source, Word, and PDF) and highlights
```

### `src/` — core simulation pipeline

| Script | Purpose |
|---|---|
| `mesh_generator.py` | 2D finite-element mesh generation for top- and edge-contact geometries |
| `inse_physics.py` | InSe material parameters, thermionic-field-emission (TFE) contact boundary conditions, and the coupled lattice-temperature equation |
| `simple_physics.py`, `simple_dd.py`, `model_create.py` | Core DEVSIM drift-diffusion/Poisson template physics |
| `run_inse_simulations.py` | Main sweep loop (geometry × work function × contact length), produces the `iv_*.csv` files and `summary_results*.csv` |

### `src/` — robustness-check studies (Section 3.6 of the manuscript)

| Script | Produces |
|---|---|
| `mesh_convergence_study.py` | `mesh_convergence.csv` — 0.5–4× mesh-density sweep |
| `richardson_extrapolation.py` | Richardson-extrapolates `mesh_convergence.csv` (and the independent cross-check `mesh_convergence_independent_check.csv`) to estimate the mesh-converged resistance and the resulting edge-vs-top advantage, addressing the top-contact mesh-convergence gap directly |
| `test_graded_mesh.py` | `graded_mesh_test.csv` — third, independent convergence check using a corner-refined (graded) mesh at the top contact's metal-semiconductor-vacuum triple point, 0.5×–8× density |
| `gamma_sensitivity_study.py` | `gamma_sensitivity.csv`, `fig_gamma_sensitivity_final.png` — tunneling-enhancement-factor sweep |
| `vdw_gap_study.py` | `vdw_gap_study.csv`, `fig_vdw_gap_final.png` — van der Waals transmission-coefficient sweep |
| `channel_length_study.py`, `replot_channel_length.py` |