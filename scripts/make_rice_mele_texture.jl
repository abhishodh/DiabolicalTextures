using LinearAlgebra
using Printf

const N_K = 4001
const N_GRID = 25
const ARROW_SCALE = 1.8

function ghat(delta::Real, J::Real, k::Real)
    gx = (1 + delta)/2 + (1 - delta)*cos(k)/2
    gy = (1 - delta)*sin(k)/2
    gz = -J
    normg = sqrt(gx^2 + gy^2 + gz^2)
    return gx/normg, gy/normg, gz/normg
end

function lower_band_state(n, gauge::Symbol)
    nx, ny, nz = n
    if gauge == :north
        denom = sqrt(2*(1 + nz))
        return (-(nx - im*ny)/denom, (1 + nz)/denom)
    end
    denom = sqrt(2*(1 - nz))
    return (-(1 - nz)/denom, (nx + im*ny)/denom)
end

function preferred_gauge(delta::Real, J::Real)
    north_margin = Inf
    south_margin = Inf
    for k in range(0, 2pi, length=258)[1:end-1]
        nz = ghat(delta, J, k)[3]
        north_margin = min(north_margin, 1 + nz)
        south_margin = min(south_margin, 1 - nz)
    end
    return north_margin >= south_margin ? :north : :south
end

function berry_phase(delta::Real, J::Real)
    hypot(delta, J) > 1e-14 || return NaN
    gauge = preferred_gauge(delta, J)
    ks = range(0, 2pi, length=N_K + 1)
    state = lower_band_state(ghat(delta, J, first(ks)), gauge)
    product = 1.0 + 0.0im
    for k in ks[2:end]
        next_state = lower_band_state(ghat(delta, J, k), gauge)
        overlap = conj(state[1])*next_state[1] + conj(state[2])*next_state[2]
        product *= overlap/abs(overlap)
        state = next_state
    end
    return mod(-angle(product), 2pi)
end

function gamma_color(gamma::Real)
    phi = mod(gamma, 2pi)
    return (
        0.58 + 0.28*cos(phi),
        0.56 + 0.22*cos(phi - 2pi/3),
        0.56 + 0.22*cos(phi + 2pi/3),
    )
end

function write_arrow(io, index, J, delta, gamma, arrow_length)
    direction = (cos(gamma), sin(gamma))
    tail = (J - arrow_length*direction[1]/2,
            delta - arrow_length*direction[2]/2)
    tip = (J + arrow_length*direction[1]/2,
           delta + arrow_length*direction[2]/2)
    perpendicular = (-direction[2], direction[1])
    head_length = 0.56*arrow_length
    head_width = 0.34*arrow_length
    color = gamma_color(gamma)
    @printf(io, "\\definecolor{rmc%d}{rgb}{%.6f,%.6f,%.6f}\n",
            index, color...)
    @printf(io, "\\fill[rmc%d,opacity=.34] (axis cs:%.7f,%.7f) circle[radius=.22pt];\n",
            index, J, delta)
    @printf(io, "\\draw[rmc%d,line width=.58pt,opacity=.95] (axis cs:%.7f,%.7f)--(axis cs:%.7f,%.7f);\n",
            index, tail..., tip...)
    for sign in (-1.0, 1.0)
        wing = (tip[1] - head_length*direction[1] + sign*head_width*perpendicular[1],
                tip[2] - head_length*direction[2] + sign*head_width*perpendicular[2])
        @printf(io, "\\draw[rmc%d,line width=.50pt,opacity=.95] (axis cs:%.7f,%.7f)--(axis cs:%.7f,%.7f);\n",
                index, tip..., wing...)
    end
end

function main()
    values = collect(range(-1.0, 1.0, length=N_GRID))
    step = values[2] - values[1]
    arrow_length = 0.50*ARROW_SCALE*step
    root = abspath(get(ENV, "TEXTURE_OUTPUT", joinpath(@__DIR__, "..", "build")))
    mkpath(root)
    open(joinpath(root, "rice_mele_texture_arrows.tex"), "w") do io
        for value in values
            @printf(io, "\\draw[gray!45,opacity=.12,line width=.25pt] (axis cs:%.7f,-1)--(axis cs:%.7f,1);\n", value, value)
            @printf(io, "\\draw[gray!45,opacity=.12,line width=.25pt] (axis cs:-1,%.7f)--(axis cs:1,%.7f);\n", value, value)
        end
        index = 0
        for J in values, delta in values
            gamma = berry_phase(delta, J)
            isfinite(gamma) || continue
            index += 1
            write_arrow(io, index, J, delta, gamma, arrow_length)
        end
    end
    println("wrote $(N_GRID^2 - 1) Rice--Mele Berry-phase arrows")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
