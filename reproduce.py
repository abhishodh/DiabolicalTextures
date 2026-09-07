"""Generate the main-text and supplemental figures."""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

GROUPS = {
    "main-rm": ["texture.pdf", "rm_phase_spectrum_mu.pdf", "Scaling.pdf"],
    "landau": ["Spectrum collapse and Wavefunctions.pdf"],
    "supp-numerics": ["rice_mele_droplet_profiles.pdf", "rice_mele_numerical_scaling.pdf",
                      "classd_numerical_scaling.pdf", "classd_numerical_spectral_flow.pdf"],
    "continuum": ["classd_droplets_edges.pdf", "classd_trap_modes.pdf"],
    "diagrams": ["classd_phase_2panel.pdf", "rice_mele_homogeneous_phases.pdf",
                 "classd_homogeneous_phases.pdf", "rice_mele_unwinding_phases.pdf",
                 "classd_unwinding_phases_norho.pdf", "suspension_construction_schematic.pdf"],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", choices=["all", *GROUPS])
    parser.add_argument("--output", type=Path, default=ROOT / "build")
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--latexmk", default="latexmk")
    args = parser.parse_args()
    out = args.output.resolve()
    if out == ROOT or (ROOT in out.parents and
                       out.parts[len(ROOT.parts)] in {"scripts", "figures"}):
        parser.error("Choose an output directory separate from the source files.")
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TEXTURE_OUTPUT=str(out), MPLBACKEND="Agg",
               MPLCONFIGDIR=str(out / ".matplotlib"), OPENBLAS_NUM_THREADS="2")

    def run(command):
        subprocess.run(list(map(str, command)), cwd=out, env=env, check=True)

    def julia(script):
        run([args.julia, f"--project={ROOT}", ROOT / "scripts" / script])

    def latex(stem):
        run([args.latexmk, "-pdf", "-interaction=nonstopmode", "-halt-on-error", stem + ".tex"])

    groups = list(GROUPS) if args.group == "all" else [args.group]
    for tex in (ROOT / "figures" / "tex").glob("*.tex"):
        shutil.copyfile(tex, out / tex.name)
    for group in groups:
        if group == "main-rm":
            julia("berry_textures.jl")
            julia("make_main_rice_mele_figures.jl")
            latex("rm_phase_spectrum_mu")
        elif group == "landau":
            julia("landau_levels.jl")
        elif group == "supp-numerics":
            julia("make_numerical_universality_figures.jl")
        elif group == "continuum":
            run([sys.executable, ROOT / "scripts" / "make_classd_supp_figures.py"])
        elif group == "diagrams":
            julia("make_rice_mele_texture.jl")
            for filename in GROUPS[group]:
                latex(Path(filename).stem)
        missing = [f for f in GROUPS[group] if not (out / f).is_file()]
        if missing:
            raise RuntimeError(f"Missing outputs for {group}: {missing}")
    print(f"Figures: {out}")


if __name__ == "__main__":
    main()
