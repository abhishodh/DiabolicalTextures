using DelimitedFiles
using LaTeXStrings
using LinearAlgebra
using PyPlot
using SparseArrays
using Statistics

# Numerical figures for the Rice--Mele and two-dimensional class-D textures.

const ROOT = abspath(get(ENV, "TEXTURE_OUTPUT", joinpath(@__DIR__, "..", "build")))
const DATA_DIR = joinpath(ROOT, "numerical_data")
mkpath(DATA_DIR)

rc("text", usetex=true)
rc("font", family="serif", serif=["Computer Modern"], size=18)
rc("axes", labelsize=20, linewidth=1.0)
rc("xtick", labelsize=16, direction="out")
rc("ytick", labelsize=16, direction="out")
rc("legend", fontsize=13, frameon=false)
rc("lines", linewidth=2.0)

const BLUE = "#0057b8"
const RED = "#c62828"
const BLACK = "#000000"
const PURPLE = "#4a3a7e"
const INK = "#26262a"

function savefig_clean(fig, stem)
    fig.savefig(joinpath(ROOT, stem * ".pdf"), bbox_inches="tight",
                metadata=Dict("CreationDate" => nothing, "ModDate" => nothing))
    fig.savefig(joinpath(DATA_DIR, stem * ".png"), dpi=180, bbox_inches="tight")
    close(fig)
end

function panel!(ax, letter; fontsize=18)
    ax.text(0.5, 1.045, "(" * letter * ")", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=fontsize)
end

"""Rice--Mele single-particle Hamiltonian in the convention of the Supplement."""
function rice_mele_hamiltonian(L::Int, alpha, beta, mu; periodic::Bool)
    @assert iseven(L)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    xmax = L

    for x in 1:L
        theta = 2pi * (1 - beta) * x / xmax
        Jx = alpha + sin(theta)
        push!(rows, x); push!(cols, x); push!(vals, (-1)^x * Jx - mu)
    end

    lastbond = periodic ? L : L - 1
    for x in 1:lastbond
        theta = 2pi * (1 - beta) * x / xmax
        delta = cos(theta)
        hopping = (1 - (-1)^x * delta) / 2
        xp = x == L ? 1 : x + 1
        push!(rows, x); push!(cols, xp); push!(vals, hopping)
        push!(rows, xp); push!(cols, x); push!(vals, hopping)
    end
    return sparse(rows, cols, vals, L, L)
end

"""Low-lying spectrum of an open Rice--Mele chain via its tridiagonal form."""
function rice_mele_open_near_zero(L::Int, alpha, beta, mu; nev=18)
    @assert iseven(L)
    diagonal = zeros(L)
    offdiagonal = zeros(L - 1)
    for x in 1:L
        theta = 2pi * (1 - beta) * x / L
        diagonal[x] = (-1)^x * (alpha + sin(theta)) - mu
        if x < L
            delta = cos(theta)
            offdiagonal[x] = (1 - (-1)^x * delta) / 2
        end
    end
    allvals = eigvals(SymTridiagonal(diagonal, offdiagonal))
    order = sortperm(abs.(allvals))
    return sort(allvals[order[1:min(nev, L)]])
end

const sx = ComplexF64[0 1; 1 0]
const sy = ComplexF64[0 -im; im 0]
const sz = ComplexF64[1 0; 0 -1]
const s0 = Matrix{ComplexF64}(I, 2, 2)
# Majorana basis (a_even, b_even, a_odd, b_odd) of the Supplement.
const G1 = -kron(s0, sx)
const G2 = -kron(sx, sz)
const G3 = -kron(s0, sy)
const G4 = -kron(sy, sz)
const GD = -im * G1 * G2

function add_block!(rows, cols, vals, block, i0, j0)
    for a in axes(block, 1), b in axes(block, 2)
        z = block[a, b]
        abs(z) < 1e-15 && continue
        push!(rows, i0 + a)
        push!(cols, j0 + b)
        push!(vals, z)
    end
end

"""
    classd_strip_hamiltonian(ncell, alpha, beta, rho, kx; periodic_y)

Exact Wilson--Dirac lattice representative of the class-D family at fixed
longitudinal momentum. The texture varies over `ncell` two-layer cells.
"""
function classd_strip_hamiltonian(ncell::Int, alpha, beta, rho, kx;
                                  periodic_y::Bool)
    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]

    B_rho = (rho / 2) * GD
    C_rho = (-im * rho / 2) * G1 * G4
    T_rho = (B_rho - im * C_rho) / 2

    denom = periodic_y ? ncell : ncell - 1
    for j in 1:ncell
        theta = 2pi * (1 - beta) * (j - 1) / denom
        m = alpha + sin(theta)
        delta = cos(theta)
        tx = (1 - m) / 2
        d1 = tx * sin(kx)
        d3 = (1 + m) / 2 - tx * cos(kx)
        ay = (1 + delta) / 2
        onsite = d1 * G1 + d3 * G3 + ay * G4 + (rho / 2) * GD
        add_block!(rows, cols, vals, onsite, 4(j - 1), 4(j - 1))
    end

    lastlink = periodic_y ? ncell : ncell - 1
    for j in 1:lastlink
        theta = 2pi * (1 - beta) * (j - 1) / denom
        delta = cos(theta)
        ty = (1 - delta) / 2
        T = -(ty / 2) * G4 - im * (ty / 2) * G2 + T_rho
        jp = j == ncell ? 1 : j + 1
        add_block!(rows, cols, vals, T, 4(j - 1), 4(jp - 1))
        add_block!(rows, cols, vals, T', 4(jp - 1), 4(j - 1))
    end
    H = sparse(rows, cols, vals, 4ncell, 4ncell)
    return H
end

function near_zero(H; nev=18)
    n = size(H, 1)
    # Dense Hermitian diagonalization is deterministic and robust at the exact
    # zero modes. The largest class-D matrix below is 1024 by 1024.
    allvals = eigvals(Hermitian(Matrix(H)))
    order = sortperm(abs.(allvals))
    return sort(real.(allvals[order[1:min(nev, n)]]))
end

function rm_zero_mode_spacing(vals)
    levels = sort(abs.(vals))
    length(levels) >= 3 || error("Need the shifted zero mode and first transverse pair")
    return mean(levels[2:3])
end

function classd_zero_mode_spacing(vals)
    levels = sort(abs.(vals))
    length(levels) >= 4 || error("Need the zero-mode pair and first transverse pair")
    return mean(levels[3:4]) - mean(levels[1:2])
end

function adjacent_positive_spacing(vals)
    positive = sort(vals[vals .> 1e-10])
    length(positive) >= 2 || error("Need two positive levels")
    return positive[2] - positive[1]
end

function power_fit(Ls, gaps)
    x, y = log.(Float64.(Ls)), log.(gaps)
    A = hcat(ones(length(x)), x)
    coeff = A \ y
    residual = y - A * coeff
    r2 = 1 - sum(abs2, residual) / sum(abs2, y .- mean(y))
    return coeff[2], exp(coeff[1]), r2
end

function scaling_data(; rice_mele_only=false)
    Lrm_half = [128, 160, 192, 256, 320, 384, 512, 640, 768, 1024]
    Lrm_side = [256, 512, 1024, 2048, 4096, 8192]
    Lrm_cap = [256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
    mu = 0.5
    beta_cap = 0.30
    alpha_cap = sin(2pi * beta_cap) + sqrt(mu^2 - cos(2pi * beta_cap)^2)

    rm_half = Float64[]
    rm_side = Float64[]
    rm_cap = Float64[]
    for L in Lrm_half
        vals = near_zero(rice_mele_hamiltonian(L, 1.0, 0.0, 0.0; periodic=true); nev=16)
        push!(rm_half, rm_zero_mode_spacing(vals))
    end
    for L in Lrm_side
        vals = rice_mele_open_near_zero(L, 1 + mu, 0.0, mu; nev=16)
        push!(rm_side, adjacent_positive_spacing(vals))
    end
    for L in Lrm_cap
        vals = rice_mele_open_near_zero(L, alpha_cap, beta_cap, mu; nev=20)
        push!(rm_cap, adjacent_positive_spacing(vals))
    end

    writedlm(joinpath(DATA_DIR, "rice_mele_gap_scaling.csv"),
             hcat(Lrm_half, rm_half), ',')
    writedlm(joinpath(DATA_DIR, "rice_mele_side_gap_scaling.csv"),
             hcat(Lrm_side, rm_side), ',')
    writedlm(joinpath(DATA_DIR, "rice_mele_cap_gap_scaling.csv"),
             hcat(Lrm_cap, rm_cap), ',')
    if rice_mele_only
        return (; Lrm_half, Lrm_side, Lrm_cap, rm_half, rm_side, rm_cap,
                  alpha_cap, beta_cap)
    end

    Ld = [32, 40, 48, 64, 80, 96, 128, 160, 192, 256]
    rho = 0.5
    classd_rho0 = Float64[]
    classd_side = Float64[]
    classd_cap = Float64[]
    alpha_d_cap = sin(2pi * beta_cap) + sqrt(rho^2 - cos(2pi * beta_cap)^2)
    for L in Ld
        vals = near_zero(classd_strip_hamiltonian(L, 1.0, 0.0, 0.0, 0.0;
                                                  periodic_y=true); nev=24)
        push!(classd_rho0, classd_zero_mode_spacing(vals))
        vals = near_zero(classd_strip_hamiltonian(L, 1 + rho, 0.0, rho, 0.0;
                                                  periodic_y=true); nev=24)
        push!(classd_side, adjacent_positive_spacing(vals))
        vals = near_zero(classd_strip_hamiltonian(L, alpha_d_cap, beta_cap, rho, 0.0;
                                                  periodic_y=false); nev=24)
        push!(classd_cap, adjacent_positive_spacing(vals))
    end

    writedlm(joinpath(DATA_DIR, "classd_gap_scaling.csv"),
             hcat(Ld, classd_rho0, classd_side, classd_cap), ',')
    return (; Lrm_half, Lrm_side, Lrm_cap, rm_half, rm_side, rm_cap, Ld, classd_rho0,
            classd_side, classd_cap, alpha_cap, alpha_d_cap, beta_cap)
end

function plot_scaling_panel(datasets, expected, xlabel, stem;
                            panel_letter=nothing, figsize=(6.4, 4.5),
                            fitted_lines=false, data_fit_legend=false,
                            fit_labels=nothing, label_x_fractions=nothing,
                            label_offsets=nothing)
    fig, ax = subplots(1, 1, figsize=figsize)
    marker_styles = ["o", "s", "D"]
    fit_curves = NamedTuple[]
    for (i, ((Ls, gaps, color, label), exponent)) in enumerate(zip(datasets, expected))
        firstfit = max(3, length(Ls) - 5)
        slope, amp_fit, r2 = power_fit(Ls[firstfit:end], gaps[firstfit:end])
        plot_label = data_fit_legend ? "_nolegend_" : label
        ax.loglog(Ls, gaps, marker_styles[i], ms=7.0,
                  mfc="none", mec=color, mew=1.6, color=color,
                  linestyle="none", label=plot_label)
        xfit = range(minimum(Ls), maximum(Ls), length=200)
        amp_ref = exp(mean(log.(gaps[firstfit:end]) .- exponent .* log.(Ls[firstfit:end])))
        fit_amp = fitted_lines ? amp_fit : amp_ref
        fit_exponent = fitted_lines ? slope : exponent
        ax.loglog(xfit, fit_amp .* xfit .^ fit_exponent,
                  color=color, lw=1.5, linestyle="--")
        push!(fit_curves, (; xmin=minimum(Ls), xmax=maximum(Ls),
                            amp=fit_amp, exponent=fit_exponent, color))
        println(stem, " ", i, ": amplitude=", amp_fit, " slope=", slope,
                " expected=", exponent, " R2=", r2)
    end
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(L"\Delta", fontsize=15)
    ax.xaxis.set_minor_formatter(PyPlot.matplotlib.ticker.NullFormatter())
    if data_fit_legend
        marker_handles = [ax.plot(Float64[], Float64[], marker=marker_styles[i],
                                  ms=6.5, color=datasets[i][3], markerfacecolor="none",
                                  linestyle="none")[1] for i in eachindex(datasets)]
        data_handle = PyPlot.pybuiltin("tuple")(marker_handles)
        fit_handle = ax.plot(Float64[], Float64[], color=INK, lw=1.5,
                             linestyle="--")[1]
        tuple_handler = PyPlot.matplotlib.legend_handler.HandlerTuple(ndivide=nothing)
        ax.legend([data_handle, fit_handle], [L"\mathrm{data}", L"\mathrm{fit}"],
                  handler_map=Dict(PyPlot.pybuiltin("tuple") => tuple_handler),
                  loc="lower left", handlelength=1.8, labelspacing=0.35,
                  frameon=true, facecolor="white", edgecolor="none", framealpha=1.0,
                  fontsize=11)
    else
        ax.legend(loc="lower left", handlelength=1.2, labelspacing=0.35,
                  frameon=true, facecolor="white", edgecolor="none", framealpha=1.0,
                  fontsize=11)
    end
    if fit_labels !== nothing
        fractions = label_x_fractions === nothing ? fill(0.55, length(fit_curves)) : label_x_fractions
        offsets = label_offsets === nothing ? fill(1.12, length(fit_curves)) : label_offsets
        fig.canvas.draw()
        for (i, (curve, fit_label)) in enumerate(zip(fit_curves, fit_labels))
            lx = exp(log(curve.xmin) + fractions[i] * log(curve.xmax / curve.xmin))
            ly = curve.amp * lx^curve.exponent
            x1, x2 = lx / 1.08, lx * 1.08
            y1 = curve.amp * x1^curve.exponent
            y2 = curve.amp * x2^curve.exponent
            p1 = Vector{Float64}(ax.transData.transform((x1, y1)))
            p2 = Vector{Float64}(ax.transData.transform((x2, y2)))
            angle = atan(p2[2] - p1[2], p2[1] - p1[1]) * 180 / pi
            ax.text(lx, ly * offsets[i], fit_label, color=curve.color,
                    fontsize=9.2, ha="center", va="bottom", rotation=angle,
                    rotation_mode="anchor")
        end
    end
    if panel_letter !== nothing
        panel!(ax, panel_letter; fontsize=14)
    end
    ax.tick_params(which="both", width=1.0, labelsize=13)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.18, top=0.97)
    savefig_clean(fig, stem)
end

function plot_scaling(data)
    datasets_rm = [
        (data.Lrm_half, data.rm_half, BLUE, L"\mathrm{Dirac\ trap}:\ L^{-1/2}"),
        (data.Lrm_side, data.rm_side, BLACK, L"\mathrm{lobe\ side}:\ L^{-1}"),
        (data.Lrm_cap, data.rm_cap, RED, L"\mathrm{lobe\ cap}:\ L^{-2/3}"),
    ]
    rm_fit_labels = [latexstring("z\\vartheta=" * string(round(-power_fit(
        Ls[max(3, length(Ls)-5):end], gaps[max(3, length(Ls)-5):end])[1]; digits=3)) *
        "\\ [" * theory * "]") for ((Ls, gaps, _, _), theory) in
        zip(datasets_rm, ["1/2", "1", "2/3"])]
    plot_scaling_panel(datasets_rm, [-1 / 2, -1.0, -2 / 3],
                       L"L", "rice_mele_numerical_scaling";
                       panel_letter="c", figsize=(4.7, 4.5),
                       fitted_lines=true, data_fit_legend=true,
                       fit_labels=rm_fit_labels,
                       label_x_fractions=[0.54, 0.59, 0.47],
                       label_offsets=[1.16, 1.18, 1.17])

    hasproperty(data, :Ld) || return
    datasets_d = [
        (data.Ld, data.classd_rho0, BLUE, L"\rho=0:\ L_y^{-1/2}"),
        (data.Ld, data.classd_side, BLACK, L"\mathrm{lobe\ side}:\ L_y^{-2/3}"),
        (data.Ld, data.classd_cap, RED, L"\mathrm{lobe\ cap}:\ L_y^{-1/2}"),
    ]
    plot_scaling_panel(datasets_d, [-1 / 2, -2 / 3, -1 / 2],
                       L"L_y", "classd_numerical_scaling";
                       panel_letter="c", figsize=(4.7, 4.5),
                       fitted_lines=true, data_fit_legend=true,
                       fit_labels=[latexstring("z\\vartheta=" * string(round(-power_fit(
                           Ls[max(3, length(Ls)-5):end], gaps[max(3, length(Ls)-5):end])[1]; digits=3)) *
                           "\\ [" * theory * "]") for ((Ls, gaps, _, _), theory) in
                           zip(datasets_d, ["1/2", "2/3", "1/2"])],
                       label_x_fractions=[0.62, 0.65, 0.34],
                       label_offsets=[1.04, 1.04, 1.04])
end

function plot_rice_mele_trap_profiles()
    fig, axes = subplots(1, 3, figsize=(15.2, 4.25))

    x = collect(range(0, 1, length=1001))
    mu = 0.5
    edge = sqrt.(2 .+ 2 .* sin.(2pi .* x)) .- mu
    inside = edge .< 0
    xc = x[inside]
    axes[1].plot(x, edge, color=BLUE, lw=2.2)
    axes[1].axhline(0, color=INK, lw=1.0)
    axes[1].fill_between(xc, edge[inside], 0, color="#b49fcd", alpha=0.9)
    width = acos((2-mu^2)/2) / pi
    xminus, xplus = 0.75-width/2, 0.75+width/2
    axes[1].axvline(xminus, color=RED, lw=1.3, linestyle=":")
    axes[1].axvline(xplus, color=RED, lw=1.3, linestyle=":")
    axes[1].annotate("", xy=(xplus, -0.46), xytext=(xminus, -0.46),
                     arrowprops=Dict("arrowstyle" => "<->", "color" => INK,
                                     "linewidth" => 1.2))
    axes[1].text((xminus + xplus) / 2, -0.50, L"w", ha="center", va="top",
                 fontsize=15)
    axes[1].text((xminus + xplus) / 2, -0.16, "metal", ha="center", va="center",
                 fontsize=14)
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(-0.62, 1.55)
    axes[1].set_xticks([0, 0.5, 1])
    axes[1].set_xlabel(L"x/L", fontsize=15)
    axes[1].set_ylabel(L"R(x)-\mu", fontsize=15)
    panel!(axes[1], "a"; fontsize=14)

    ax = axes[2]
    interface_color = "#bd592b"
    xleft, xright = 0.26, 0.74
    edgewidth = 0.055
    ax.axvspan(0, xleft, color="#a6d5e2", alpha=0.85, lw=0)
    ax.axvspan(xleft, xright, color="#b49fcd", alpha=0.90, lw=0)
    ax.axvspan(xright, 1, color="#a6d5e2", alpha=0.85, lw=0)
    ax.axvspan(xleft-edgewidth, xleft+edgewidth, color=interface_color, alpha=0.23, lw=0)
    ax.axvspan(xright-edgewidth, xright+edgewidth, color=interface_color, alpha=0.23, lw=0)
    ax.axvline(xleft, color=interface_color, lw=1.35, linestyle=":")
    ax.axvline(xright, color=interface_color, lw=1.35, linestyle=":")
    ax.text(0.12, 0.37, "insulator", ha="center", va="center", fontsize=12)
    ax.text(0.50, 0.72, "metal", ha="center", va="center", fontsize=13)
    ax.text(0.50, 0.61, L"\xi_{\rm min}\sim L^{1/2}",
            ha="center", va="center", fontsize=10.5)
    ax.text(0.50, 0.50, L"(z,p,\vartheta)=(2,2,1/2)",
            ha="center", va="center", fontsize=10)
    ax.text(0.88, 0.37, "insulator", ha="center", va="center", fontsize=12)
    ax.text(xleft, 0.42, L"x_-", ha="center", va="top", fontsize=12)
    ax.text(xright, 0.42, L"x_+", ha="center", va="top", fontsize=12)
    ax.annotate("", xy=(xright, 0.24), xytext=(xleft, 0.24),
                arrowprops=Dict("arrowstyle" => "<->", "color" => INK,
                                "linewidth" => 1.2))
    ax.text(0.50, 0.205, L"w\sim L", ha="center", va="top", fontsize=13)
    ax.annotate("", xy=(xleft+edgewidth, 0.92), xytext=(xleft-edgewidth, 0.92),
                arrowprops=Dict("arrowstyle" => "<->", "color" => interface_color,
                                "linewidth" => 1.15))
    ax.text(xleft, 0.84, L"\xi_{\rm int}\sim L^{1/3}",
            ha="center", va="center", fontsize=9.5)
    ax.text(xleft, 0.76, L"(z,p,\vartheta)=(2,1,1/3)",
            ha="center", va="center", fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ["top", "right", "bottom", "left"]
        ax.spines[side].set_visible(false)
    end
    panel!(ax, "b"; fontsize=14)

    u2 = collect(range(-3.2, 4.6, length=900))
    dens2 = exp.(-u2 .^ 2) ./ sqrt(pi)
    umax, ngrid = 10.0, 900
    h = umax / ngrid
    u1 = collect(range(h, umax - h, length=ngrid - 1))
    H1 = SymTridiagonal(2 / h^2 .+ u1, fill(-1 / h^2, ngrid - 2))
    state1 = eigen(H1, 1:1).vectors[:, 1]
    state1 ./= sqrt(sum(abs2, state1) * h)
    dens1 = abs2.(state1)
    axes[3].plot(u2, dens2, color=BLUE, lw=2.1,
                 label=L"p=2\ \mathrm{(full\ line)}")
    axes[3].plot(u1, dens1, color=RED, lw=2.1,
                 label=L"p=1\ \mathrm{(half\ line)}")
    axes[3].set_xlim(-3.2, 5.0)
    axes[3].set_ylim(0, 0.68)
    axes[3].set_xlabel(L"u", fontsize=15)
    axes[3].set_ylabel(L"|\psi_0(u)|^2", fontsize=15)
    axes[3].legend(loc="upper right", frameon=false, fontsize=11)
    panel!(axes[3], "c"; fontsize=14)

    for ax in axes
        ax.tick_params(which="both", width=1.0, labelsize=13)
    end
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.18, top=0.86, wspace=0.31)
    savefig_clean(fig, "rice_mele_droplet_profiles")
end

function spectral_flow_data(; ncell=48, nalpha=121, nev=22)
    alphas = collect(range(-2, 2, length=nalpha))
    flows = Dict{Float64, Matrix{Float64}}()
    for rho in [0.0, 0.5]
        spectrum = zeros(nalpha, nev)
        for (i, alpha) in enumerate(alphas)
            vals = near_zero(classd_strip_hamiltonian(ncell, alpha, 0.0, rho, 0.0;
                                                      periodic_y=true); nev=nev)
            spectrum[i, :] = vals
        end
        flows[rho] = spectrum
        writedlm(joinpath(DATA_DIR, "classd_spectral_flow_rho_" *
                          replace(string(rho), "." => "p") * ".csv"),
                 hcat(alphas, spectrum), ',')
    end
    return alphas, flows
end

function plot_spectral_flow(alphas, flows)
    fig, axes = subplots(1, 2, figsize=(12.4, 4.7), sharey=true)
    for (ax, rho) in zip(axes, [0.0, 0.5])
        spectrum = flows[rho]
        for n in 1:size(spectrum, 2)
            energies = spectrum[:, n]
            positive = energies .>= 0
            ax.scatter(alphas[positive], energies[positive], s=12,
                       color=BLUE, alpha=0.32, edgecolors="none")
            ax.scatter(alphas[.!positive], energies[.!positive], s=12,
                       color=RED, alpha=0.32, edgecolors="none")
        end
        critical = rho == 0 ? [-1.0, 1.0] : [-1.5, -0.5, 0.5, 1.5]
        for ac in critical
            ax.axvline(ac, color=INK, linestyle="--", linewidth=1.35, alpha=0.9)
        end
        ax.axhline(0, color=INK, linestyle=":", linewidth=1.2)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-0.82, 0.82)
        ax.set_xticks([-2, -1, 0, 1, 2])
        ax.set_xlabel(L"\alpha")
        panel_label = rho == 0 ? L"(a)\ \rho=0" : L"(b)\ \rho=1/2"
        ax.text(0.5, 1.045, panel_label, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=18)
    end
    axes[1].set_ylabel(L"\varepsilon_n(k_x=0)")
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.18, top=0.86, wspace=0.16)
    savefig_clean(fig, "classd_numerical_spectral_flow")
end

function main()
    if "--rice-mele-only" in ARGS
        plot_scaling(scaling_data(; rice_mele_only=true))
        return
    end
    data = scaling_data()
    plot_scaling(data)
    plot_rice_mele_trap_profiles()
    alphas, flows = spectral_flow_data()
    plot_spectral_flow(alphas, flows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
