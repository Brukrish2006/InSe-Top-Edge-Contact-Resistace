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
| `gamma_sensitivity_study.py` | `gamma_sensitivity.csv`, `fig_gamma_sensitivity_final.png` — tunneling-enhancement-factor sweep |
| `vdw_gap_study.py` | `vdw_gap_study.csv`, `fig_vdw_gap_final.png` — van der Waals transmission-coefficient sweep |
| `channel_length_study.py`, `replot_channel_length.py` | `channel_length_study.csv`, `fig_channel_length.png` — channel-length scaling / SCLC check |
| `material_sensitivity_study.py` | `material_sensitivity.csv`, `fig_material_sensitivity.png` — hole-mobility and thermal-conductivity sweeps |
| `doping_study.py` | `doping_study.csv` — background-doping sweep |
| `sclc_density_check.py` | `sclc_density_check.csv` — center-channel carrier-density vs. bias (direct SCLC confirmation) |
| `sclc_test.py` | `sclc_vds_sweep*.csv` — resistance vs. bias voltage at fixed channel length |
| `verify_lc_collapse.py` | `verify_lc_collapse.csv` — mesh-refinement check of the self-heating length-collapse result |

### `src/` — plotting

| Script | Produces |
|---|---|
| `plot_results.py` | `fig1_iv_curves.png`, `fig2_resistance_vs_wf_v2.png`, `fig3_resistance_vs_lc_v2.png` |
| `plot_thermal.py` | `thermal_comparison.png` |
| `run_sim.ps1` | Optional Windows/PowerShell convenience launcher for the main sweep |

### `data/`

All CSV outputs referenced by the manuscript's tables and figures: the per-configuration I-V
curves (`iv_edge_WF*.csv`, `iv_top_Lc*_WF*.csv`), the two summary tables (`summary_results.csv`,
`summary_results_thermal.csv` — Tables 1–2), and every robustness-check CSV listed above, including
`mesh_convergence_independent_check.csv` — a second, independently-run top-contact mesh sequence
(0.5×–8×, different DEVSIM installation) used to cross-validate the Richardson extrapolation in
Section 3.6.

### `figures/`

The 8 figures used in the manuscript, matching the filenames cited in `paper/draft_paper_FINAL.md`:
`fig1_iv_curves.png`, `fig2_resistance_vs_wf_v2.png`, `fig3_resistance_vs_lc_v2.png`,
`thermal_comparison.png`, `fig_gamma_sensitivity_final.png`, `fig_vdw_gap_final.png`,
`fig_channel_length.png`, `fig_material_sensitivity.png`.

### `paper/`

- `draft_paper_FINAL.md` — the manuscript source (Markdown, with LaTeX-style math).
- `InSe_Contact_Resistance_Manuscript_v12.docx` / `.pdf` — the current formatted manuscript.
- `highlights.md` — the 3–5 bullet "Highlights" required by some journal submission systems.

Earlier numbered manuscript drafts and internal review notes are intentionally excluded from
this repository (see `.gitignore`) — only the current version is tracked.

## Running the simulations

Requires [DEVSIM](https://devsim.org/) and a working Python 3 environment (numpy, pandas,
matplotlib).

```
cd src
python run_inse_simulations.py
```

This writes one `iv_*.csv` file per (geometry, work function, contact length) configuration and
aggregates the extracted contact resistance — and, for the thermal run, peak lattice
temperature — into `summary_results.csv` and `summary_results_thermal.csv`. Each robustness-check
script in the table above can be run independently once the main sweep has completed, and each
plotting script regenerates its corresponding figure(s) from the CSVs in `data/`.

## Citing this repository

If you use this code or data, please cite the accompanying manuscript. We recommend archiving a
tagged release of this repository via [Zenodo](https://zenodo.org/) to obtain a permanent,
citable DOI, and adding it here once available.

## License

Released under the [MIT License](LICENSE).
