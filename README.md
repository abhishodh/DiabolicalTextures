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

The table maps manuscript panels to the generated PDFs and their source panels.

| Figure | Group | Generated source |
| --- | --- | --- |
| Main 1(a,b) | `main-rm` | `texture.pdf` |
| Main 1(c) | `main-rm` | `rice_mele_spectral_panels.pdf`, upper panel (spectral flow in alpha) |
| Main 1(d) | `main-rm` | `Scaling.pdf`, panel (a) (neutral-gap scaling) |
| Main 1(e) | `main-rm` | `rice_mele_spectral_panels.pdf`, lower panel (spectral flow in beta) |
| Main 1(f,g) | `main-rm` | `rm_phase_spectrum_mu.pdf`, upper and lower right-hand phase diagrams |
| Main 2 | `diagrams` | `classd_phase_2panel.pdf` |
| End Matter 3(a,b) | `main-rm` | `Scaling.pdf`, panels (b,c) (local observable and entanglement) |
| Supplemental 1 | `landau` | `Spectrum collapse and Wavefunctions.pdf` |
| Supplemental 2 | `diagrams` | `rice_mele_homogeneous_phases.pdf` |
| Supplemental 3 | `supp-numerics` | `rice_mele_droplet_profiles.pdf` |
| Supplemental 4(a,b) | `diagrams` | `rice_mele_unwinding_phases.pdf` |
| Supplemental 4(c) | `supp-numerics` | `rice_mele_numerical_scaling.pdf` |
| Supplemental 5 | `diagrams` | `classd_homogeneous_phases.pdf` |
| Supplemental 6 | `continuum` | `classd_droplets_edges.pdf` |
| Supplemental 7 | `supp-numerics` | `classd_numerical_spectral_flow.pdf` |
| Supplemental 8(a,b) | `diagrams` | `classd_unwinding_phases_norho.pdf` |
| Supplemental 8(c) | `supp-numerics` | `classd_numerical_scaling.pdf` |

The generators in `scripts/` are:

- `main-rm`: `berry_textures.jl` and `make_main_rice_mele_figures.jl`.
- `landau`: `landau_levels.jl`.
- `supp-numerics`: `make_numerical_universality_figures.jl`.
- `continuum`: `make_classd_supp_figures.py`.
- `diagrams`: `make_rice_mele_texture.jl` and the sources in `figures/tex/`.

The `main-rm` group also compiles `figures/tex/rm_phase_spectrum_mu.tex`.
For Main 1, keep (a,b) stacked at left, place (c,d,e) at equal height in
the upper-right row, and place equally sized (f,g) in the row below.
Extract and relabel the source panels according to the table. The gap
panel 1(d) shows the trap-critical scaling $\Delta\sim1/\sqrt{L}$.
After assembling Main 1, prepare its embedded fonts for inclusion in LaTeX:

```sh
python scripts/prepare_pdf_fonts.py build/figure1_layout.pdf
```

For End Matter 3, place only the local-observable and entanglement panels
side by side, relabeled (a,b). For Supplemental 4 and 8, place the phase diagram and
scaling plot side by side at 66% and 32% of the line width, respectively.
