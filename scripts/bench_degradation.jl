# Benchmark + numerical-identity harness for the update_Q! energy-degradation kernel.
#
# Usage:
#   julia --project=. -t 4 scripts/bench_degradation.jl baseline <E_max>
#   julia --project=. -t 4 scripts/bench_degradation.jl compare  <E_max>
#
# "baseline" : runs a full solve, serializes the realistic input flux Ie and the
#              Q array produced by a controlled update_Q! re-sweep (the reference).
# "compare"  : reloads the SAME Ie, re-sweeps update_Q!, and compares Q against the
#              reference -> isolates the kernel change (identical input).
#
# Both modes print the summed wall-clock time of update_Q! over the full energy loop.

using AURORA
using Serialization: serialize, deserialize

const MODE  = length(ARGS) >= 1 ? ARGS[1] : "baseline"
const E_MAX = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 500.0
const TAG   = "emax$(Int(E_MAX))"
const IE_FILE = "/tmp/degr_Ie_$(TAG).bin"
const Q_FILE  = "/tmp/degr_Qref_$(TAG).bin"

const REFROOT = joinpath(@__DIR__, "..", "test", "regression", "reference_results")
const MSIS = joinpath(REFROOT, "msis_20051008-2200_70N-19E.txt")
const IRI  = joinpath(REFROOT, "iri_20051008-2200_70N-19E.txt")

function build_sim()
    altitude_lims = [100, 400]
    θ_lims = 180:-30:0
    B_angle_to_zenith = 13
    model = AuroraModel(altitude_lims, θ_lims, E_MAX, MSIS, IRI, B_angle_to_zenith)
    savedir = mktempdir()
    flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                     beams=1, z_source=1000.0)
    sim = AuroraSimulation(model, flux, savedir;
                           mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                  CFL_number=128, n_loop=1))
    AURORA.initialize!(sim; force_recompute=true, verbose=false)
    return sim
end

# Controlled re-sweep: rebuild Q from a fixed Ie, timing only update_Q!.
function sweep_updateQ!(sim, Ie; nrep=3)
    cache = AURORA.get_cache(sim)
    model = sim.model
    matrices = cache.matrices
    n_E = model.energy_grid.n
    best = Inf
    for _ in 1:nrep
        fill!(matrices.Q, 0.0)
        t = 0.0
        for iE in n_E:-1:1
            B2B = AURORA.update_matrices!(matrices, model, iE, cache.B2B_fragment)
            t += @elapsed AURORA.update_Q!(matrices, Ie, model, cache.t_loop,
                                           B2B, iE, cache.degradation)
        end
        best = min(best, t)
    end
    return best, matrices.Q
end

sim = build_sim()
cache = AURORA.get_cache(sim)
model = sim.model
n_E = model.energy_grid.n
sz = size(cache.Ie)
println("E_max=$(E_MAX)  n_E=$(n_E)  Ie size=$(sz)  n_species=$(length(model.species))  threads=$(Threads.nthreads())")

if MODE == "baseline"
    run!(sim)                                # realistic Ie at every energy
    Ie = copy(cache.Ie)
    serialize(IE_FILE, Ie)
    t, Q = sweep_updateQ!(sim, Ie)
    serialize(Q_FILE, copy(Q))
    println("BASELINE update_Q! total = $(round(t*1e3, digits=2)) ms")
    println("  sum(Q)=$(sum(Q))  sum(abs Q)=$(sum(abs, Q))")
elseif MODE == "compare"
    Ie = deserialize(IE_FILE)::Array{Float64,3}
    @assert size(Ie) == sz "Ie size mismatch; regenerate baseline for this E_max"
    Qref = deserialize(Q_FILE)::Array{Float64,3}
    t, Q = sweep_updateQ!(sim, Ie)
    println("COMPARE  update_Q! total = $(round(t*1e3, digits=2)) ms")
    absdiff = maximum(abs.(Q .- Qref))
    denom = max.(abs.(Q), abs.(Qref), eps())
    reldiff = maximum(abs.(Q .- Qref) ./ denom)
    println("  max |ΔQ|       = $(absdiff)")
    println("  max rel ΔQ     = $(reldiff)")
    println("  sum(Q)=$(sum(Q))  sum(Qref)=$(sum(Qref))")
else
    error("unknown mode $MODE")
end
