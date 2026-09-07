# Occupied-band Wilson loops for the two Rice--Mele texture panels.
include("make_rice_mele_texture.jl")
using PyPlot, LaTeXStrings, DelimitedFiles

function texture_panels()
    root = abspath(get(ENV, "TEXTURE_OUTPUT", joinpath(@__DIR__, "..", "build")))
    mkpath(joinpath(root, "numerical_data"))
    rc("text", usetex=true)
    rc("font", family="serif", serif=["Computer Modern"], size=17)
    fig = figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 3], hspace=0.42)
    a, b = fig.add_subplot(gs.__getitem__(0)), fig.add_subplot(gs.__getitem__(1))
    phasecolor(g) = (0.58+0.28cos(g), 0.56+0.22cos(g-2pi/3), 0.56+0.22cos(g+2pi/3))
    rows = zeros(42, 3)
    index = 0
    for (alpha, height) in [(0.0, 0.8), (2.0, 0.25)]
        for s in range(0, 1, length=21)
            gamma = berry_phase(cos(2pi*s), alpha+sin(2pi*s))
            index += 1
            rows[index,:] = [alpha, s, gamma]
            a.quiver(s, height, cos(gamma), sin(gamma), angles="uv",
                     scale_units="inches", scale=5.5, width=0.004,
                     pivot="mid", color=phasecolor(gamma))
        end
        a.text(1.08, height, latexstring("\\alpha=", Int(alpha)), va="center")
    end
    a.set(xlim=(-0.04,1.28), ylim=(-0.1,1.1), yticks=[], xticks=[0,.5,1])
    a.set_xlabel(L"t/T\ \mathrm{or}\ x/L", labelpad=0)
    for side in ["left","right","top"]
        a.spines[side].set_visible(false)
    end
    a.text(-.04,1.03,"(a)",transform=a.transAxes)
    plane = Vector{Vector{Float64}}()
    for J in range(-3,3,length=41), delta in range(-1,1,length=17)
        gamma = berry_phase(delta,J)
        isfinite(gamma) || continue
        push!(plane,[J,delta,gamma])
        b.quiver(J, delta, cos(gamma), sin(gamma), angles="uv",
                 scale_units="inches", scale=7.0, width=.0028,
                 pivot="mid", color=phasecolor(gamma))
    end
    t=range(0,2pi,length=721)
    for (alpha, color) in [(0.,"#0057b8"),(2.,"#c62828")]
        b.plot(alpha .+ sin.(t), cos.(t), color=color,lw=1.5)
    end
    b.plot([0],[0],"ko",ms=3)
    b.set(xlim=(-3.15,3.15),ylim=(-1.13,1.13),xlabel=L"J",ylabel=L"\delta",
          xticks=[-3,-2,-1,0,1,2,3],yticks=[-1,0,1],aspect="equal")
    b.text(-.04,1.05,"(b)",transform=b.transAxes)
    fig.subplots_adjust(left=.1,right=.98,bottom=.1,top=.95)
    fig.savefig(joinpath(root,"texture.pdf"))
    close(fig)
    writedlm(joinpath(root,"numerical_data","berry_rows.csv"),rows,',')
    writedlm(joinpath(root,"numerical_data","berry_plane.csv"),reduce(vcat,permutedims.(plane)),',')
end

if abspath(PROGRAM_FILE) == @__FILE__
    texture_panels()
end
