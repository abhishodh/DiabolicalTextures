using LinearAlgebra
using DelimitedFiles
using LaTeXStrings
using PyPlot

# Rice--Mele spectral flow, neutral gap, staggered density, and entropy.
const ROOT = abspath(get(ENV, "TEXTURE_OUTPUT", joinpath(@__DIR__, "..", "build")))
const DATA = joinpath(ROOT, "numerical_data", "main_rice_mele")
mkpath(DATA)
BLAS.set_num_threads(2)

const BLUE = "#0057b8"
const RED = "#c62828"
const BLACK = "#202020"
const SPECTRUM_L = 200
const FLOW_SAMPLES = 70
const SCALING_SIZES = collect(400:400:7600)
const DENSITY_SIZES = [1000, 2000, 3000, 4000, 5000]

rc("text", usetex=true)
rc("font", family="serif", serif=["Computer Modern"], size=16)
rc("axes", labelsize=18, linewidth=0.9)
rc("xtick", labelsize=16, direction="out")
rc("ytick", labelsize=16, direction="out")
rc("legend", fontsize=14, frameon=false)

function rice_mele(L, alpha, beta=0.0; periodic=true)
    @assert iseven(L)
    x = collect(1:L)
    theta = 2pi * (1 - beta) .* x / L
    potential = (-1.0).^x .* (alpha .+ sin.(theta))
    hopping = (1 .- (-1.0).^x .* cos.(theta)) ./ 2
    h = SymTridiagonal(potential, hopping[1:end-1])
    # For beta=0 the terminal bond vanishes, so this tridiagonal matrix
    # also represents the periodic chain exactly.
    if !periodic || abs(hopping[end]) < 1e-14
        return h
    end
    dense = Matrix(h)
    dense[1,end] = dense[end,1] = hopping[end]
    return Symmetric(dense)
end

function flow_data()
    alpha = collect(range(-2, 2; length=FLOW_SAMPLES))
    beta = collect(range(0, 1; length=FLOW_SAMPLES))
    ea = reduce(vcat, [permutedims(eigvals(rice_mele(SPECTRUM_L, a; periodic=false))) for a in alpha])
    eb = reduce(vcat, [permutedims(eigvals(rice_mele(SPECTRUM_L, 0.0, b; periodic=false))) for b in beta])
    writedlm(joinpath(DATA, "flow_alpha_L200.csv"), hcat(alpha, ea), ',')
    writedlm(joinpath(DATA, "flow_beta_L200.csv"), hcat(beta, eb), ',')
    return alpha, beta, ea, eb
end

function scaling_data()
    gaps = Dict{Int,Float64}()
    entropies = Dict{Int,Float64}()
    densities = Dict{Int,Matrix{Float64}}()
    for L in sort(union(SCALING_SIZES, DENSITY_SIZES))
        F = eigen(rice_mele(L, 1.0))
        nocc = count(<(0), F.values)
        gaps[L] = F.values[nocc+1] - F.values[nocc]
        if L in DENSITY_SIZES
            density = vec(sum(abs2, @view(F.vectors[:,1:nocc]); dims=2))
            sites = collect(1:2:L-1)
            x = (sites .- 3L/4) ./ sqrt(L)
            magnitude = sqrt(L) .* abs.(density[sites] .- density[sites .+ 1])
            densities[L] = hcat(x, magnitude)
            writedlm(joinpath(DATA, "density_L$(L).csv"), densities[L], ',')
        end
        if L in SCALING_SIZES
            # The ground state is pure. The last quarter has the same
            # entropy as the first 3L/4 sites, at a much lower matrix cost.
            V = @view F.vectors[3L÷4+1:L,1:nocc]
            nu = eigvals(Symmetric(V * V'))
            nu = clamp.(nu, 1e-15, 1-1e-15)
            entropies[L] = -sum(nu .* log.(nu) .+ (1 .- nu) .* log1p.(-nu))
        end
        println("L=", L, ": gap=", gaps[L])
        flush(stdout)
        GC.gc()
    end
    table = hcat(SCALING_SIZES, [gaps[L] for L in SCALING_SIZES],
                 [entropies[L] for L in SCALING_SIZES])
    writedlm(joinpath(DATA, "gap_entropy.csv"), table, ',')
    return table, densities
end

function spectrum_panel!(ax, parameter, levels, xlabel, title)
    for n in axes(levels, 2)
        energy = levels[:,n]
        occupied = energy .< 0
        ax.scatter(parameter[occupied], energy[occupied]; s=7.5, color=RED,
                   alpha=0.18, edgecolors="none")
        ax.scatter(parameter[.!occupied], energy[.!occupied]; s=7.5, color=BLUE,
                   alpha=0.18, edgecolors="none")
    end
    # The level indexed L/2+1 tracks the extra occupied state and its expulsion.
    charged = levels[:,SPECTRUM_L÷2+1]
    for (mask, color) in [(charged .<= 0, RED), (charged .>= 0, BLUE)]
        y = copy(charged)
        y[.!mask] .= NaN
        ax.plot(parameter, y; color=color, linewidth=1.5)
        ax.scatter(parameter[mask], charged[mask]; s=8, color=color, edgecolors="none")
    end
    ax.axhline(0; color=BLACK, linewidth=0.9, linestyle="--", zorder=0)
    ax.set_xlabel(xlabel; fontsize=22, labelpad=2)
    ax.set_ylabel(L"\varepsilon_n"; fontsize=22, labelpad=3)
    ax.set_title(title; fontsize=22, pad=7)
    ax.tick_params(labelsize=20, length=3, pad=3)
end

function plot_flow(alpha, beta, ea, eb)
    # Vertically stacked spectral-flow panels (c,d).
    fig = figure(figsize=(360/72, 398.479/72))
    ac = fig.add_axes([0.20, 0.625, 0.76, 0.285])
    ad = fig.add_axes([0.20, 0.145, 0.76, 0.285])
    spectrum_panel!(ac, alpha, ea, L"\alpha", L"(c)\;\beta=0")
    spectrum_panel!(ad, beta, eb, L"\beta", L"(d)\;\alpha=0")
    for a in [-1,1]
        ac.axvline(a; color=BLACK, linewidth=0.9, linestyle="--", zorder=0)
    end
    ac.set_xlim(-2.12,2.12)
    ac.set_ylim(-3.45,3.45)
    ac.set_xticks([-2,-1,0,1,2])
    ac.set_yticks([-2,0,2])
    ad.set_xlim(-0.04,1.04)
    ad.set_ylim(-1.58,1.58)
    ad.set_xticks([0,0.5,1])
    ad.set_xticklabels([L"0",L"1/2",L"1"])
    ad.set_yticks([-1,0,1])
    fig.savefig(joinpath(ROOT,"rice_mele_spectral_panels.pdf"))
    fig.savefig(joinpath(DATA,"spectral_panels.png"); dpi=180)
    close(fig)
end

function plot_scaling(table, densities)
    fig, axs = subplots(1,3; figsize=(7.0,3.05),
                       gridspec_kw=Dict("width_ratios"=>[1.0,1.15,1.0]))
    fig.subplots_adjust(left=0.085,right=0.985,bottom=0.245,top=0.83,wspace=0.57)
    for (ax, letter) in zip(axs, ["a","b","c"])
        ax.set_title("("*letter*")"; fontsize=17, pad=8)
        ax.tick_params(labelsize=14, length=3, pad=3)
    end
    Ls = table[:,1]
    invL = 1 ./ sqrt.(Ls)
    ax = axs[1]
    ax.plot([0,0.052], sqrt(4pi).*[0,0.052]; color=BLACK, linestyle="--", linewidth=1.3)
    ax.plot(invL,table[:,2],"o"; color=BLUE, ms=4.0, label="data")
    ax.set_xlabel(L"1/\sqrt{L}"; fontsize=17, labelpad=3)
    ax.set_ylabel(L"\Delta"; fontsize=18, labelpad=1)
    ax.set_xlim(-0.002,0.053)
    ax.set_ylim(-0.007,0.197)
    ax.set_xticks([0,0.025,0.05])
    ax.set_xticklabels([L"0",L"0.025",L"0.05"])
    ax.set_yticks([0,0.1,0.2])
    ax.set_yticklabels([L"0",L"0.1",L"0.2"])
    ax.text(0.035,0.91,L"\sqrt{4\pi/L}"; transform=ax.transAxes,
            fontsize=15, ha="left",va="top")

    ax = axs[2]
    markers = ["o","s","^","v","D"]
    colors = [BLUE,RED,BLACK,BLUE,RED]
    for (L,marker,color) in zip(DENSITY_SIZES,markers,colors)
        d = densities[L]
        keep = abs.(d[:,1]) .<= 0.55
        ax.plot(d[keep,1],d[keep,2]; marker=marker, ls="none", ms=3.1,
                color=color, markerfacecolor="none", markeredgewidth=0.6,
                label=string(L))
    end
    ax.set_xlabel(L"(x-x^*)/\sqrt{L}"; fontsize=16, labelpad=3)
    ax.set_ylabel(L"\sqrt{L}\,|\langle\hat M(x)\rangle|"; fontsize=16, labelpad=1)
    ax.set_xlim(-0.55,0.55)
    ax.set_ylim(0.40,2.18)
    ax.set_xticks([-0.5,0,0.5])
    ax.set_xticklabels([L"-0.5",L"0",L"0.5"])
    ax.set_yticks([0.5,1,1.5])
    ax.set_yticklabels([L"0.5",L"1",L"1.5"])
    ax.legend(loc="upper center", ncol=2, fontsize=14,
              handlelength=0.6, handletextpad=0.25, columnspacing=0.45,
              labelspacing=0.2, borderaxespad=0.1)

    ax = axs[3]
    S = table[:,3]
    offset = sum(S .- log.(0.75.*Ls)./12)/length(Ls)
    lf = exp.(range(log(minimum(Ls)),log(maximum(Ls));length=150))
    ax.plot(lf, offset .+ log.(0.75.*lf)./12; color=BLACK, ls="--", linewidth=1.3)
    ax.plot(Ls,S,"o"; color=BLUE, ms=4)
    ax.set_xscale("log")
    ax.set_xticks([400,1600,6400])
    ax.set_xticklabels([L"400",L"1600",L"6400"])
    ax.set_xlim(330,9200)
    ax.set_xlabel(L"L"; fontsize=17, labelpad=3)
    ax.set_ylabel(L"S(3L/4)"; fontsize=17, labelpad=2)
    ax.set_yticks([0.75,0.85,0.95,1.05])
    ax.tick_params(axis="x",which="minor",bottom=false)
    fig.savefig(joinpath(ROOT,"Scaling.pdf"))
    fig.savefig(joinpath(DATA,"scaling.png"); dpi=210)
    close(fig)
end

function main()
if "--flow-only" in ARGS
    plot_flow(flow_data()...)
elseif "--plot-only" in ARGS
    aa = readdlm(joinpath(DATA,"flow_alpha_L200.csv"),',',Float64)
    bb = readdlm(joinpath(DATA,"flow_beta_L200.csv"),',',Float64)
    table = readdlm(joinpath(DATA,"gap_entropy.csv"),',',Float64)
    densities = Dict(L=>readdlm(joinpath(DATA,"density_L$(L).csv"),',',Float64) for L in DENSITY_SIZES)
    plot_flow(aa[:,1],bb[:,1],aa[:,2:end],bb[:,2:end])
    plot_scaling(table,densities)
else
    plot_flow(flow_data()...)
    table,densities = scaling_data()
    plot_scaling(table,densities)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
