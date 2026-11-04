# Gate-Controlled Hole Injection at Top and Edge Contacts to Few-Layer InSe (v3)

[![verification](https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace/actions/workflows/verify.yml/badge.svg)](https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace/actions/workflows/verify.yml)

**Current version: v3** (`src_v3/`, `data_v3/`, `figures_v3/`, `paper_v3/`). v3 extends the corrected
v2 model with:
- a back gate;
- Fermi-Dirac statistics;
- anisotropic mobility and thermal conductivity;
- nonlocal Tsu-Esaki/WKB tunnelling;
- temperature-dependent statistics;
- interface, Peltier and tunnelling heat, and transient heating;
- layer-dependent InSe band structure;
- separate top/edge pinning;
- a one-parameter calibration on Au/MoS2 with out-of-sample tests;
- a compact model and Monte Carlo uncertainty propagation.

v2 (`src_v2/` ...) and the original v1 (`src/` ...) are kept unchanged. See `CHANGES_v2.md` and `CHANGES_v3.md`.

## v3 quick start
```bash
pip install -r requirements.txt
export DEVSIM_MATH_LIBS=/usr/lib/x86_64-linux-gnu/liblapack.so.3:/usr/lib/x86_64-linux-gnu/libblas.so.3   # Linux
cd src_v3
python verification.py && python check_verification.py      # V1-V9 (also run by CI)
python studies.py wf pinning aniso compact thickness gated thermal uq mos2   # several hours on 2 cores
python make_figures.py
cd ../paper_v3 && pdflatex InSe_v3.tex && pdflatex InSe_v3.tex
```

| v3 file | Purpose |
|---|---|
| `src_v3/tcad2d.py` | Model: meshes (top/edge, optional oxide + gate), anisotropic FD drift-diffusion, flux contacts, nonlocal or local tunnelling, full heat model, transient, extraction |
| `src_v3/verification.py`, `check_verification.py` | V1-V9 benchmarks and CI tolerance check |
| `src_v3/studies.py` | All v3 studies (InSe and the MoS2 calibration) |
| `src_v3/make_figures.py` | Figures 1-10, compact-model fit, TLM fits, UQ statistics (`data_v3/derived_*`) |
| `paper_v3/` | `InSe_v3.tex/.pdf`, `Response_to_Reviewers_v3.tex/.pdf` |

---

# Hole Injection at Top and Edge Metal Contacts to Few-Layer InSe (v2, revised)

> **Revision note (November 2026).** The first version of this study (`src/`, `data/`, `figures/`,
> `paper/`) was rejected by the *Journal of Computational Electronics*. Reviewer 1 correctly
> identified implementation errors in the public code. The model has been rewritten and verified
> in `src_v2/`, every result was regenerated in `data_v2/` and `figures_v2/`, and the revised
> manuscript and the point-by-point response are in `paper_v2/`. See `CHANGES_v2.md` for the list of
> defects fixed and the conclusions withdrawn. The v1 folders are kept unchanged for transparency;
> **do not use v1 results.**

## v2 quick start

```bash
pip install devsim numpy scipy pandas matplotlib
# Linux: point DEVSIM at a BLAS/LAPACK (Windows/MKL users: set DEVSIM_MATH_LIBS to mkl_rt)
export DEVSIM_MATH_LIBS=/usr/lib/x86_64-linux-gnu/liblapack.so.3:/usr/lib/x86_64-linux-gnu/libblas.so.3
cd src_v2
python verification.py                  # V1-V5 benchmarks -> data_v2/verification_*.csv
python studies.py mesh wf wf_iso pinning tunnel doping material thermalbc lch iv profiles
python make_figures.py                  # -> figures_v2/
cd ../paper_v2 && pdflatex InSe_revised.tex && pdflatex Response_to_Reviewers.tex
```
The full set of studies takes about 15 min on two cores.

| v2 file | Purpose |
|---|---|
| `src_v2/inse2d.py` | Model: meshes, full Poisson coupling, p-type doping, flux (Robin) thermionic/TFE contacts with self-consistent field-dependent Gamma(F), WKB vdW-gap transmission, Fermi-level pinning, quasi-Fermi Joule heating, Robin thermal BCs, exact R decomposition, energy balance |
| `src_v2/verification.py` | Ohmic resistor, finite-velocity contact resistor, Mott-Gurney SCLC, Schottky diode, 2D energy balance |
| `src_v2/studies.py` | All parameter studies used in the paper |
| `src_v2/make_figures.py` | Figures 1-8 and S1 |
| `data_v2/` | All CSV outputs (+ `profiles.npz`: centre-line and 2D profiles) |
| `paper_v2/` | `InSe_revised.tex/.pdf`, `Response_to_Reviewers.tex/.pdf` |

---

## Original (v1) README

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