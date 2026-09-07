# Topological Phenomena Protected by Diabolical Textures

Figure-generation code for v2 of [arXiv:2605.30421](https://arxiv.org/abs/2605.30421).

## Setup

Requires Julia 1.12, Python 3.12 or newer, and a TeX installation with
`latexmk`, `pdflatex`, TikZ/PGFPlots, `standalone`, and Computer Modern fonts.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
julia --project=. -e 'using Pkg; ENV["PYTHON"]=Sys.which("python"); Pkg.instantiate(); Pkg.build("PyCall")'
```

## Generate figures

Run from this directory. Figures and numerical tables are written to `build/`.
All system sizes and parameter values are specified in the scripts.

```sh
python reproduce.py all
```

To generate a subset, replace `all` with a group in the table below.
Use `--output PATH` to choose another output directory.

| Figure | Group | Output PDF |
| --- | --- | --- |
| Main 1(a,b) | `main-rm` | `texture.pdf` |
| Main 1(c–e) | `main-rm` | `rm_phase_spectrum_mu.pdf` |
| Main 2 | `main-rm` | `Scaling.pdf` |
| Main 3 | `diagrams` | `classd_phase_2panel.pdf` |
| Supplemental 1 | `landau` | `Spectrum collapse and Wavefunctions.pdf` |
| Supplemental 2 | `diagrams` | `rice_mele_homogeneous_phases.pdf` |
| Supplemental 3 | `supp-numerics` | `rice_mele_droplet_profiles.pdf` |
| Supplemental 4(a,b) | `diagrams` | `rice_mele_unwinding_phases.pdf` |
| Supplemental 4(c) | `supp-numerics` | `rice_mele_numerical_scaling.pdf` |
| Supplemental 5 | `diagrams` | `classd_homogeneous_phases.pdf` |
| Supplemental 6 | `continuum` | `classd_droplets_edges.pdf` |
| Supplemental 7 | `supp-numerics` | `classd_numerical_spectral_flow.pdf` |
| Supplemental 8 | `continuum` | `classd_trap_modes.pdf` |
| Supplemental 9(a,b) | `diagrams` | `classd_unwinding_phases_norho.pdf` |
| Supplemental 9(c) | `supp-numerics` | `classd_numerical_scaling.pdf` |
| Supplemental 10 | `diagrams` | `suspension_construction_schematic.pdf` |

The generators in `scripts/` are:

- `main-rm`: `berry_textures.jl` and `make_main_rice_mele_figures.jl`.
- `landau`: `landau_levels.jl`.
- `supp-numerics`: `make_numerical_universality_figures.jl`.
- `continuum`: `make_classd_supp_figures.py`.
- `diagrams`: `make_rice_mele_texture.jl` and the sources in `figures/tex/`.

The `main-rm` group also compiles `figures/tex/rm_phase_spectrum_mu.tex`.
For Main 1, place `texture.pdf` and `rm_phase_spectrum_mu.pdf` side by side
at equal height. For Supplemental 4 and 9, place the phase diagram and
scaling plot side by side at 66% and 32% of the line width, respectively.
