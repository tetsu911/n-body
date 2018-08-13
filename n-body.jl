using Formatting: printfmt
using BenchmarkTools

mutable struct Star
    name::String
    r::Vector{Float64}
    v::Vector{Float64}
    m::Float64
end

const SOLAR_MASS = 4 * π^2
const DAYS_PER_YEAR = 365.24

function mk_pairs(bodies::Vector{Star})
    local ret = Vector{Tuple{Star, Star}}[]
    @inbounds for (i, x) in enumerate(bodies)
        @inbounds for y in bodies[i+1:end]
            ret = vcat(ret, (x, y))
        end
    end
    return ret
end

function core_advance_v!(dt::Float64, pairs::Vector{Tuple{Star, Star}})
    @inbounds for (s1, s2) in pairs
        local dx = s1.r - s2.r
        local mag = dt * ((dx[1]^2 + dx[2]^2 + dx[3]^2) ^ (-1.5))
        local b1m = s1.m * mag
        local b2m = s2.m * mag
        local s1.v = s1.v - dx * b2m
        local s2.v = s2.v + dx * b1m
    end
end

function core_advance_r!(dt::Float64, bodies::Vector{Star})
    @inbounds for s in bodies
        s.r += dt * s.v
    end
end

function advance!(dt::Float64, n::Int64,
                  bodies::Vector{Star},
                  pairs::Vector{Tuple{Star, Star}})
    @inbounds for i in 1:n
        core_advance_v!(dt, pairs)
        core_advance_r!(dt, bodies)
    end
end

function core_u_energy(pairs::Vector{Tuple{Star, Star}},
                     ene::Float64=0.0)
    @inbounds for (s1, s2) in pairs
        local dx = s1.r - s2.r
        ene -= (s1.m * s2.m) / ((dx[1]^2 + dx[2]^2 + dx[3]^2) ^ 0.5)
    end
    return ene
end

function core_k_energy(bodies::Vector{Star},
                       ene::Float64=0.0)
    @inbounds for s in bodies
        ene += s.m * (s.v[1]^2 + s.v[2]^2 + s.v[3]^2) / 2.
    end
    return ene
end

function report_energy(bodies::Vector{Star},
                       pairs::Vector{Tuple{Star, Star}},
                       ene::Float64=0.0)
    ene = core_u_energy(pairs, ene)
    ene = core_k_energy(bodies, ene)
    printfmt("{:.9f}\n", ene)
    return ene
end

function core_offset_momentum!(bodies::Vector{Star},
                              pr::Vector{Float64}=[0.0, 0.0, 0.0])
    @inbounds for s in bodies
        pr -= s.v * s.m
    end
    return pr
end

function offset_momentum!(ref::Star,
                          bodies::Vector{Star},
                          pr::Vector{Float64}=[0.0, 0.0, 0.0])
    core_offset_momentum!(bodies, pr)
    ref.v = pr / ref.m
end

function main(n::Int, bodies::Vector{Star}, ref::Int64=1)
    pairs::Vector{Tuple{Star, Star}} = mk_pairs(bodies)
    offset_momentum!(bodies[ref], bodies)
    report_energy(bodies, pairs)
    advance!(0.01, n, bodies, pairs)
    report_energy(bodies, pairs)
end

function mk_stars()
    local stars = [
        Star("sun",     [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], SOLAR_MASS),

        Star("jupiter", [4.84143144246472090e+00,
                        -1.16032004402742839e+00,
                        -1.03622044471123109e-01],
                        [1.66007664274403694e-03 * DAYS_PER_YEAR,
                         7.69901118419740425e-03 * DAYS_PER_YEAR,
                        -6.90460016972063023e-05 * DAYS_PER_YEAR],
                         9.54791938424326609e-04 * SOLAR_MASS),

        Star("saturn",  [8.34336671824457987e+00,
                         4.12479856412430479e+00,
                        -4.03523417114321381e-01],
                        [-2.76742510726862411e-03 * DAYS_PER_YEAR,
                          4.99852801234917238e-03 * DAYS_PER_YEAR,
                          2.30417297573763929e-05 * DAYS_PER_YEAR],
                          2.85885980666130812e-04 * SOLAR_MASS),

        Star("uranus",  [1.28943695621391310e+01,
                        -1.51111514016986312e+01,
                        -2.23307578892655734e-01],
                        [2.96460137564761618e-03 * DAYS_PER_YEAR,
                         2.37847173959480950e-03 * DAYS_PER_YEAR,
                        -2.96589568540237556e-05 * DAYS_PER_YEAR],
                         4.36624404335156298e-05 * SOLAR_MASS),

        Star("neptune", [1.53796971148509165e+01,
                        -2.59193146099879641e+01,
                         1.79258772950371181e-01],
                        [2.68067772490389322e-03 * DAYS_PER_YEAR,
                         1.62824170038242295e-03 * DAYS_PER_YEAR,
                        -9.51592254519715870e-05 * DAYS_PER_YEAR],
                         5.15138902046611451e-05 * SOLAR_MASS) ]
    return stars
end
stars = mk_stars()
main(10, stars)
stars = mk_stars()
@benchmark main(1_000_000, stars)

main(1, stars)


#=
-0.169289903
-0.169311556
-0.169406769
-0.169313524
-0.169858047
-0.169895192
-0.170315077
-0.170306906
-0.170367944
-0.170343215
-0.170394214
-0.170397347
-0.170417747
-0.170496314


BenchmarkTools.Trial:
  memory estimate:  6.26 GiB
  allocs estimate:  60000541
  --------------
  minimum time:     3.827 s (8.64% GC)
  median time:      3.831 s (8.71% GC)
  mean time:        3.831 s (8.71% GC)
  maximum time:     3.835 s (8.79% GC)
  --------------
  samples:          2
  evals/sample:     1
=#
