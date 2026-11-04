# Gate-Controlled Hole Injection at Top and Edge Contacts to Few-Layer InSe

[![verification](https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace/actions/workflows/verify.yml/badge.svg)](https://github.com/Brukrish2006/InSe-Top-Edge-Contact-Resistace/actions/workflows/verify.yml)

A 2D DEVSIM TCAD study of hole injection into few-layer InSe through top and edge metal contacts (`src_v3/`, `data_v3/`, `figures_v3/`, `paper_v3/`). The model includes:
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

## Quick start
```bash
pip install -r requirements.txt
export DEVSIM_MATH_LIBS=/usr/lib/x86_64-linux-gnu/liblapack.so.3:/usr/lib/x86_64-linux-gnu/libblas.so.3   # Linux
cd src_v3
python verification.py && python check_verification.py      # V1-V9 (also run by CI)
python studies.py wf pinning aniso compact thickness gated thermal uq mos2   # several hours on 2 cores
python make_figures.py && python make_schematic.py && python make_maps.py
cd ../paper_v3 && pdflatex InSe_v3.tex && pdflatex InSe_v3.tex
```

| File | Purpose |
|---|---|
| `src_v3/tcad2d.py` | Model: meshes (top/edge, optional oxide + gate), anisotropic FD drift-diffusion, flux contacts, nonlocal or local tunnelling, full heat model, transient, extraction |
| `src_v3/verification.py`, `check_verification.py` | V1-V9 benchmarks and CI tolerance check |
| `src_v3/studies.py` | All v3 studies (InSe and the MoS2 calibration) |
| `src_v3/make_figures.py` | Figures 1-10, compact-model fit, TLM fits, UQ statistics (`data_v3/derived_*`) |
| `src_v3/make_schematic.py`, `make_maps.py` | Fig. 0 (geometries) and Fig. 11 (2D hole-density and temperature maps; fields in `data_v3/fields_*.csv`) |
| `paper_v3/` | `InSe_v3.tex/.pdf` and its included tables and text |
