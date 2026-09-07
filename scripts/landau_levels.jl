"""Rice--Mele Landau-level spectra and signed sublattice amplitudes."""
module RM
include("make_main_rice_mele_figures.jl")
end
using LinearAlgebra, DelimitedFiles, PyPlot, LaTeXStrings

function landau_levels()
    root = RM.ROOT
    mkpath(joinpath(root,"numerical_data"))
    sizes = collect(1000:200:2000)
    fig = figure(figsize=(15.2,4.2))
    gs = fig.add_gridspec(1,4,width_ratios=[1.15,1,1,1],wspace=.43)
    spectral = fig.add_subplot(gs.__getitem__((0,0)))
    modes = [fig.add_subplot(gs.__getitem__((0,j))) for j in 1:3]
    colors=["#0057b8","#c62828","#202020"]
    markers=["o","s","D","^","v","x"]
    lattice_handles = []
    for (ii,L) in enumerate(sizes)
        h=RM.rice_mele(L,1.0)
        # The near-zero mode is slightly shifted at finite L. Use its index
        # to define n=0, then keep both signs of the transverse ladder.
        ev=eigvals(h)
        i0=argmin(abs.(ev))
        levels=eigen(h,i0-25:i0+25)
        spectral.plot(-25:25,sqrt(L).*levels.values,markers[ii],
                      ms=3.8,mfc="none",alpha=.55,mew=1.0,zorder=3,
                      color=colors[mod1(ii,3)],label=latexstring("L=",L))
        writedlm(joinpath(root,"numerical_data","landau_spectrum_L$(L).csv"),
                 hcat(-25:25,levels.values),',')
        if L==2000
            # Positive modes n=0,1,2. Remove the fast Bloch factor (-1)^cell.
            x=collect(1:L); u=(x .-3L/4)./sqrt(L)
            sublatticephase=(-1.).^fld.(x.-1,2)
            for n in 0:2
                psi=levels.vectors[:,26+n].*sublatticephase
                ell=sqrt(L/(2pi)); q=(x .-3L/4)./ell
                f0=exp.(-q.^2/2)./sqrt(sqrt(pi)*ell)
                fs=[f0,sqrt(2).*q.*f0,(2 .*q.^2 .-1)./sqrt(2).*f0]
                # Apply the sigma_z basis transformation in the Landau-level
                # correspondence. In the lattice Bloch gauge both components
                # of a positive-energy mode have the same oscillator sign.
                target=[isodd(i) ? (n==0 ? 0. : fs[n][i]) :
                        (n==0 ? sqrt(2)*f0[i] : fs[n+1][i]) for i in x]
                target ./= norm(target)
                psi .*= sign(dot(psi,target))
                for (parity,color,marker) in [(1,colors[1],"o"),(0,colors[2],"s")]
                    ids=findall(i->mod(i,2)==parity && abs(u[i])<=1.35,x)
                    handle=modes[n+1].plot(u[ids],L^.25 .* psi[ids],
                        linestyle="none",marker=marker,ms=4.5,mfc="none",
                        mec=color,mew=1.1,alpha=.5,zorder=3)[1]
                    n==0 && push!(lattice_handles,handle)
                end
                writedlm(joinpath(root,"numerical_data","landau_mode_n$(n)_L$(L).csv"),
                         hcat(x,psi),',')
            end
        end
    end
    n=collect(-25:25)
    spectral.plot(n,sign.(n).*sqrt.(4pi.*abs.(n)),color="#202020",lw=1.3,
                  zorder=2,label=L"\mathrm{analytic}")
    spectral.set(xlabel=L"n",ylabel=L"\sqrt{L}\,\varepsilon_n")
    spectral.legend(fontsize=9.5,loc="upper left",labelspacing=.22)
    spectral.text(.5,1.025,"(a)",transform=spectral.transAxes,ha="center")
    u=collect(range(-1.35,1.35,length=801))
    q=sqrt(2pi).*u
    f0=2^.25 .* exp.(-q.^2/2)
    f1=sqrt(2).*q.*f0
    f2=(2 .*q.^2 .-1)./sqrt(2).*f0
    fields=[(zero.(u),sqrt(2).*f0),(f0,f1),(f1,f2)]
    analytic_handles = []
    for j in 1:3
        for (field,color) in zip(fields[j],colors[1:2])
            handle=modes[j].plot(u,field,color=color,lw=1.65,zorder=2)[1]
            j==1 && push!(analytic_handles,handle)
        end
        ax=modes[j]
        ax.set(xlim=(-1.35,1.35),ylim=(-1.45,1.8),xlabel=L"(x-x^*)/\sqrt{L}")
        ax.text(.5,1.035,"("*string(Char('a'+j))*")",transform=ax.transAxes,ha="center")
        ax.tick_params(labelsize=12)
        ax.set_title(latexstring("n=",j-1),loc="right",fontsize=12)
    end
    modes[1].set_ylabel(L"L^{1/4}\psi(x)")
    fig.legend([analytic_handles[1],lattice_handles[1],analytic_handles[2],lattice_handles[2]],
               [L"A\ \mathrm{analytic}",L"A\ \mathrm{lattice}",
                L"B\ \mathrm{analytic}",L"B\ \mathrm{lattice}"],
               loc="upper center",bbox_to_anchor=(.635,1.0),ncol=4,
               fontsize=11.5,frameon=false,handlelength=1.5,columnspacing=1.4)
    fig.subplots_adjust(left=.055,right=.99,bottom=.20,top=.82)
    fig.savefig(joinpath(root,"Spectrum collapse and Wavefunctions.pdf"))
    close(fig)
end

if abspath(PROGRAM_FILE)==@__FILE__
    landau_levels()
end
