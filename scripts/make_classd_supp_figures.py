"""Vector figures for the Class-D Supplement, from the stated continuum theory.

The edge envelopes and dispersions are the local linear-mass solutions, not
finite-lattice data. The quadratic-trap eigenfunctions are computed by a
converged finite-difference solution of the dimensionless squared Dirac operator.
"""

from pathlib import Path
import os


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(os.environ.get("TEXTURE_OUTPUT", Path(__file__).resolve().parents[1] / "build")).resolve()
ROOT.mkdir(parents=True, exist_ok=True)
DATA_BLUE = "#0057b8"
DATA_RED = "#c62828"
DATA_BLACK = "#000000"
P = {
    "trivial": "#a6d5e2", "kitaev": "#dc8b8e", "plus": "#b49fcd",
    "minus": "#e4c891", "boundary": "#4a3a7e", "ink": "#26262a",
    "lower": "#23698b", "upper": "#bd592b", "neutral": "#f0f2f3",
}
plt.rcParams.update({
    "text.usetex": True, "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 10, "legend.fontsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.7, "lines.linewidth": 1.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "savefig.transparent": False, "pdf.fonttype": 42,
})


def panel(ax, letter, detail=""):
    label = f"({letter})" + (f" {detail}" if detail else "")
    ax.text(0.5, 1.065, label, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=11)


def width_fraction(alpha, rho):
    a = np.abs(np.asarray(alpha, dtype=float))
    r = abs(rho)
    if not 0 <= r < 1:
        raise ValueError("The interior-droplet formula here assumes |rho| < 1.")
    inside = (a > 1-r) & (a < 1+r)
    q = (1+a*a-r*r)/(2*np.maximum(a, 1e-15))
    return np.where(inside, np.arccos(np.clip(q, -1, 1))/np.pi, 0.0)


def quadratic_trap(n=2400, extent=8.0, states=3):
    du = 2*extent/(n+1)
    u = np.linspace(-extent, extent, n+2)[1:-1]
    diag = 2/du**2 + u**4 - 2*u
    off = np.full(n-1, -1/du**2)
    matrix = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    vals, vecs = np.linalg.eigh(matrix)
    vals, vecs = vals[:states], vecs[:, :states]
    vecs /= np.sqrt(du)
    density = (vecs[:, 0]**2 + vecs[::-1, 0]**2)/2
    return u, np.sqrt(vals), density, du


def droplets(save):
    rho, alpha, ly, v = 0.5, 1.0, 512.0, 0.5
    fraction = float(width_fraction(alpha, rho))
    ym, yp = 0.75-fraction/2, 0.75+fraction/2
    gradient = 2*np.pi/ly * np.sqrt(1-rho*rho/4)
    ell = np.sqrt(v/gradient)
    ell_fraction = ell/ly
    fig, axes = plt.subplots(1, 4, figsize=(10.2, 3.15))
    for ax, chirality, letter in zip(axes[:2], [1, -1], "ab"):
        ax.set_facecolor(P["trivial"])
        ax.axhspan(ym, yp, color=P["plus"] if chirality > 0 else P["minus"])
        for y, sign, color in [(ym, -1, P["lower"]), (yp, 1, P["upper"])]:
            ax.axhspan(y-ell_fraction, y+ell_fraction, color=color,
                       alpha=0.22, lw=0, zorder=2)
            ax.axhline(y, color=color, lw=1.4)
            for x in [0.22, 0.66]:
                dx = 0.15*sign*chirality
                ax.annotate("", xy=(x+dx, y), xytext=(x-dx, y),
                            arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.7})
            ax.plot([0.075, 0.075], [y-ell_fraction, y+ell_fraction],
                    color=color, lw=1.0)
            ax.plot([0.064, 0.086], [y-ell_fraction, y-ell_fraction],
                    color=color, lw=1.0)
            ax.plot([0.064, 0.086], [y+ell_fraction, y+ell_fraction],
                    color=color, lw=1.0)
        ax.text(0.5, 0.75, r"$p_x+i p_y$" if chirality > 0 else r"$p_x-i p_y$",
                ha="center", va="center", fontsize=12)
        ax.text(0.5, 0.28, r"$C=0$", ha="center")
        ax.text(0.5, 0.94, r"$\rho=1/2$" if chirality > 0 else r"$\rho=-1/2$",
                ha="center")
        ax.annotate("", xy=(0.94, yp), xytext=(0.94, ym),
                    arrowprops={"arrowstyle": "|-|", "lw": 0.8, "color": P["ink"]})
        ax.text(0.915, 0.75, r"$w$", ha="right", va="center")
        ax.text(0.105, 0.635, r"$\ell_c$", ha="left", va="center", fontsize=9)
        ax.text(0.5, 0.53, r"$\ell_c\sim L_y^{1/2}$"+"\n"+r"$\vartheta=1/2$",
                ha="center", va="center", fontsize=9)
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel=r"$x/L_x$", ylabel=r"$y/L_y$")
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        panel(ax, letter)
    ax = axes[2]
    z = np.linspace(0.5, 1, 4000)
    ax.axvspan(ym, yp, color=P["plus"], alpha=0.33, lw=0)
    for yc, color, label in [(ym, DATA_BLUE, r"$y_-$"), (yp, DATA_RED, r"$y_+$")]:
        density = ly/(np.sqrt(np.pi)*ell)*np.exp(-((z-yc)*ly/ell)**2)
        ax.plot(z, density, color=color, label=label)
        ax.axvline(yc, color=color, ls=":", lw=0.7)
        ax.plot([yc-ell_fraction, yc+ell_fraction], [8.0, 8.0],
                color=color, lw=1.0)
        ax.plot([yc-ell_fraction, yc-ell_fraction], [7.35, 8.65],
                color=color, lw=1.0)
        ax.plot([yc+ell_fraction, yc+ell_fraction], [7.35, 8.65],
                color=color, lw=1.0)
        ax.text(yc, 10.2, r"$\ell_c$", color=color,
                ha="center", va="bottom", fontsize=9)
    ax.set(xlim=(0.52, 0.98), ylim=(0, 50), xlabel=r"$y/L_y$", ylabel=r"$L_y|f(y)|^2$")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=2, handlelength=1.4)
    panel(ax, "c")
    ax = axes[3]
    k = np.linspace(-0.45, 0.45, 500)
    for n in [1, 2, 3]:
        e = np.sqrt((v*k)**2+2*n*v*gradient)
        ax.plot(k, e, color=DATA_BLACK, lw=0.9)
        ax.plot(k, -e, color=DATA_BLACK, lw=0.9)
    ax.plot(k, -v*k, color=DATA_BLUE, label=r"$y_-$")
    ax.plot(k, v*k, color=DATA_RED, label=r"$y_+$")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.set(xlim=(-0.45, 0.45), ylim=(-0.3, 0.3), xlabel=r"$k_x$", ylabel=r"$E$")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28),
              ncol=2, handlelength=1.4)
    panel(ax, "d")
    for ax in axes:
        ax.yaxis.labelpad = 1.5
        ax.tick_params(axis="y", pad=2)
    fig.subplots_adjust(left=0.05, right=0.985, bottom=0.28, top=0.84, wspace=0.48)
    save(fig, "classd_droplets_edges")
    return {"rho": rho, "alpha": alpha, "Ly": ly, "width_fraction": fraction,
            "edge_positions_fraction": [ym, yp], "interface_gradient": gradient,
            "interface_length": ell, "interface_length_fraction": ell_fraction,
            "first_transverse_gap": float(np.sqrt(2*v*gradient))}


def trap_figure(save, u, energies, density):
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.2))
    ax = axes[0]
    ax.plot(u, np.exp(-u*u)/np.sqrt(np.pi), color=DATA_BLUE,
            label=r"$p=1,\ \vartheta=1/2$")
    ax.plot(u, density, color=DATA_RED, label=r"$p=2,\ \vartheta=2/3$")
    ax.set(xlim=(-3.5, 3.5), ylim=(0, 0.62), xlabel=r"$u=(y-y_*)/\ell$",
           ylabel=r"$\ell |f(y)|^2$")
    ax.legend(frameon=False, loc="upper right", fontsize=8.5)
    panel(ax, "a")
    ax = axes[1]
    k = np.linspace(-1.6, 1.6, 500)
    for sign in [-1, 1]:
        ax.plot(k, sign*k, color=DATA_BLUE, label=r"$p=1$" if sign == 1 else None)
        ax.plot(k, sign*np.sqrt(k*k+2), color=DATA_BLUE, ls="--", lw=1.0, alpha=0.6)
        ax.plot(k, sign*np.sqrt(k*k+energies[0]**2), color=DATA_RED,
                label=r"$p=2$" if sign == 1 else None)
        ax.plot(k, sign*np.sqrt(k*k+energies[1]**2), color=DATA_RED, ls="--", lw=1.0, alpha=0.6)
    ax.set(xlim=(-1.6, 1.6), ylim=(-2.6, 2.6), xlabel=r"$k_x\ell$", ylabel=r"$E\ell/v$")
    ax.legend(frameon=False, loc="upper center", ncol=2)
    panel(ax, "b")
    fig.subplots_adjust(left=0.095, right=0.98, bottom=0.2, top=0.85, wspace=0.28)
    save(fig, "classd_trap_modes")


def main():
    def save(fig, stem):
        fig.savefig(ROOT/f"{stem}.pdf", metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)

    u, energies, density, _ = quadratic_trap()
    droplets(save)
    trap_figure(save, u, energies, density)


if __name__ == "__main__":
    main()
