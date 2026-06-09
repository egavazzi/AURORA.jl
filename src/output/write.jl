using NCDatasets: NCDataset, defDim, defVar, sync
using JLD2: @save
using Dates: now
import TOML
import LibGit2

# ──────────────────────────────────────────────────────────────────────────────
# TOML config
# ──────────────────────────────────────────────────────────────────────────────

function write_config_toml(sim::AuroraSimulation)
    model = sim.model
    mode = sim.mode

    commit_hash = if isdir(joinpath(pkgdir(AURORA), ".git"))
        LibGit2.head(pkgdir(AURORA))
    else
        "Not available"
    end

    config = Dict{String, Any}(
        "aurora_version"     => string(pkgversion(AURORA)),
        "commit_hash"        => commit_hash,
        "altitude_lims_km"   => [model.altitude_grid.h[1]/1e3, model.altitude_grid.h[end]/1e3],
        "theta_lims"         => collect(Float64, model.pitch_angle_grid.θ_lims),
        "E_max_eV"           => Float64(model.energy_grid.E_max),
        "B_angle_to_zenith"  => Float64(model.B_angle_to_zenith),
    )

    if mode isa TimeDependentMode
        merge!(config, Dict{String,Any}(
            "mode"       => "time_dependent",
            "duration_s" => Float64(mode.duration),
            "dt_s"       => Float64(mode.dt),
        ))
    elseif mode isa SteadyStateMode && !isnothing(mode.duration)
        merge!(config, Dict{String,Any}(
            "mode"       => "steady_state_multi_step",
            "duration_s" => Float64(mode.duration),
            "dt_s"       => Float64(mode.dt),
        ))
    else
        config["mode"] = "steady_state_single_step"
    end

    savefile = joinpath(sim.output.savedir, "config.toml")
    open(savefile, "w") do f
        TOML.print(f, config)
    end
end


# ──────────────────────────────────────────────────────────────────────────────
# Atmosphere NetCDF
# ──────────────────────────────────────────────────────────────────────────────

function write_atmosphere_nc(sim::AuroraSimulation)
    model = sim.model
    ionosphere = model.ionosphere
    savefile = joinpath(sim.output.savedir, "inputs", "atmosphere.nc")
    dl = sim.output.deflatelevel

    NCDataset(savefile, "c") do ds
        defDim(ds, "altitude", length(model.altitude_grid.h))

        alt_v = defVar(ds, "altitude", Float64, ("altitude",); deflatelevel=dl,
                       attrib=["units" => "m", "long_name" => "altitude"])
        alt_v[:] = model.altitude_grid.h

        ne_v = defVar(ds, "ne", Float64, ("altitude",); deflatelevel=dl,
                      attrib=["units" => "m-3", "long_name" => "electron number density"])
        ne_v[:] = ionosphere.ne

        Te_v = defVar(ds, "Te", Float64, ("altitude",); deflatelevel=dl,
                      attrib=["units" => "K", "long_name" => "electron temperature"])
        Te_v[:] = ionosphere.Te

        for sp in model.species
            nsp_v = defVar(ds, "n" * String(sp.name), Float64, ("altitude",); deflatelevel=dl,
                           attrib=["units" => "m-3",
                                   "long_name" => String(sp.name) * " number density"])
            nsp_v[:] = sp.density
        end
    end
end


# ──────────────────────────────────────────────────────────────────────────────
# Physics state JLD2
# ──────────────────────────────────────────────────────────────────────────────

function write_physics_jld2(sim::AuroraSimulation)
    savefile = joinpath(sim.output.savedir, "inputs", "physics_state.jld2")
    model = sim.model
    @save savefile model
end


# ──────────────────────────────────────────────────────────────────────────────
# simulation_data.nc — create, append, write input flux
# ──────────────────────────────────────────────────────────────────────────────

"""
    create_simulation_nc(sim) → NCDataset

Create `<savedir>/simulation_data.nc` with all dimensions and static variables.
Returns the open dataset; the caller is responsible for calling `close(ds)`.

The input boundary flux is written immediately into the file (before any solver loops run) as
a separate `Ie_input` variable on its own fixed `time_input` dimension.  The unlimited `time`
dimension is left empty at creation time and is populated by subsequent
[`append_chunk_nc!`](@ref) calls.
"""
function create_simulation_nc(sim::AuroraSimulation)
    out   = sim.output
    model = sim.model
    dl    = out.deflatelevel

    n_z       = length(model.altitude_grid.h)
    μ_lims    = model.pitch_angle_grid.μ_lims
    n_μ       = length(μ_lims) - 1
    μ_center  = model.pitch_angle_grid.μ_center
    Ω_beam    = beam_weight(model.pitch_angle_grid.θ_lims)
    n_E       = model.energy_grid.n
    E_centers = model.energy_grid.E_centers
    E_edges   = model.energy_grid.E_edges
    ΔE        = model.energy_grid.ΔE
    h         = model.altitude_grid.h

    commit_hash = if isdir(joinpath(pkgdir(AURORA), ".git"))
        LibGit2.head(pkgdir(AURORA))
    else
        "Not available"
    end

    nc_path = joinpath(out.savedir, "simulation_data.nc")
    ds = NCDataset(nc_path, "c",
                   attrib=["aurora_version" => string(pkgversion(AURORA)),
                            "commit_hash"   => commit_hash,
                            "creation_time" => string(now())])

    # Dimensions
    defDim(ds, "altitude",          n_z)
    defDim(ds, "pitch_angle",       n_μ)
    defDim(ds, "energy",            n_E)
    defDim(ds, "energy_bounds",      n_E + 1)
    defDim(ds, "pitch_angle_bounds", n_μ + 1)
    defDim(ds, "time",              Inf)   # unlimited — appended per loop

    # Coordinates
    alt_v = defVar(ds, "altitude", Float64, ("altitude",); deflatelevel=dl,
                   attrib=["units"     => "m",
                            "long_name" => "altitude"])
    alt_v[:] = h

    pa_v = defVar(ds, "pitch_angle", Float64, ("pitch_angle",); deflatelevel=dl,
                  attrib=["units"     => "1",
                           "long_name" => "cosine of pitch angle (beam center)"])
    pa_v[:] = μ_center

    en_v = defVar(ds, "energy", Float64, ("energy",); deflatelevel=dl,
                  attrib=["units"     => "eV",
                           "long_name" => "energy bin center"])
    en_v[:] = E_centers

    defVar(ds, "time", Float64, ("time",); deflatelevel=dl,
           attrib=["units"     => "s",
                    "long_name" => "simulation time"])
    # time values are appended by append_chunk_nc!

    ee_v = defVar(ds, "energy_edges", Float64, ("energy_bounds",); deflatelevel=dl,
                  attrib=["units"     => "eV",
                           "long_name" => "energy bin edges"])
    ee_v[:] = E_edges

    ml_v = defVar(ds, "mu_lims", Float64, ("pitch_angle_bounds",); deflatelevel=dl,
                  attrib=["units"     => "1",
                           "long_name" => "pitch-angle cosine bin boundaries"])
    ml_v[:] = collect(Float64, μ_lims)

    # Bin width and beam solid angle — saved so analysts need not recompute them.
    # `Ie` can be turned into a differential flux by dividing by `dE` and `beam_weight`.
    de_v = defVar(ds, "dE", Float64, ("energy",); deflatelevel=dl,
                  attrib=["units"     => "eV",
                           "long_name" => "energy bin width"])
    de_v[:] = ΔE

    bw_v = defVar(ds, "beam_weight", Float64, ("pitch_angle",); deflatelevel=dl,
                  attrib=["units"     => "sr",
                           "long_name" => "solid-angle beam weight"])
    bw_v[:] = Ω_beam

    # Main output variable — chunked and compressed.
    # NOTE: Ie is the electron *number* flux, already integrated over each energy bin and
    # over each beam's solid angle (not a differential flux). Units are m-2 s-1.
    defVar(ds, "Ie", Float64, ("altitude", "pitch_angle", "time", "energy");
           deflatelevel=dl,
           chunksizes=(n_z, n_μ, 1, n_E),
           attrib=["units"     => "m-2 s-1",
                    "long_name" => "electron number flux"])

    # Input flux — written once at startup before any solver loop
    t_top, Ie_top_3D = input_flux_at_save_cadence(sim)
    n_t_top = length(t_top)

    defDim(ds, "time_input", n_t_top)
    t_top_v = defVar(ds, "time_input", Float64, ("time_input",); deflatelevel=dl,
                     attrib=["units"     => "s",
                              "long_name" => "input flux time"])
    t_top_v[:] = t_top

    Ietop_v = defVar(ds, "Ie_input", Float64, ("pitch_angle", "time_input", "energy");
                     deflatelevel=dl,
                     chunksizes=(n_μ, 1, n_E),
                     attrib=["units"     => "m-2 s-1",
                              "long_name" => "input boundary number flux (precipitation, top of atmosphere)"])
    Ietop_v[:, :, :] = Ie_top_3D
    sync(ds)

    return ds
end


"""
    input_flux_at_save_cadence(sim) → (t_top, Ie_top_3D)

Return the input boundary flux subsampled to the save cadence as
`[n_μ, n_t, n_E]`, together with the corresponding time axis (in seconds).

- **Single-step steady-state**: one slice at t = 0.
- **Multi-step steady-state**: all steps at the uniform time grid.
- **Time-dependent**: subsampled from the internal grid to the save cadence via
  `1 : CFL_factor : end`.
"""
function input_flux_at_save_cadence(sim::AuroraSimulation)
    cache = sim.cache
    time  = sim.time

    if time isa RefinedTimeGrid
        t_top     = collect(Float64, time.t_save)
        Ie_top_3D = cache.Ie_top[:, 1:time.CFL_factor:end, :]
    elseif time isa UniformTimeGrid
        t_top     = collect(Float64, time.t)
        Ie_top_3D = cache.Ie_top
    else  # SingleStepConfig
        t_top     = [0.0]
        Ie_top_3D = cache.Ie_top[:, 1:1, :]
    end

    return t_top, Ie_top_3D
end


"""
    append_chunk_nc!(ds, Ie_chunk, t_chunk, sim)

Reshape `Ie_chunk` from `[n_z·n_μ, n_t, n_E]` to `[n_z, n_μ, n_t, n_E]` and
append it — together with `t_chunk` — to the unlimited `time` dimension of the
open `NCDataset` `ds`.  Calls `sync` after writing so output is crash-safe.
"""
function append_chunk_nc!(ds::NCDataset, Ie_chunk, t_chunk, sim::AuroraSimulation)
    model  = sim.model
    μ_lims = model.pitch_angle_grid.μ_lims
    z      = model.altitude_grid.h
    E_centers = model.energy_grid.E_centers

    t_vec = collect(Float64, t_chunk)
    n_new = length(t_vec)

    # Reshape 3D → 4D
    Ie_4D = restructure_Ie_from_3D_to_4D(Ie_chunk, μ_lims, z, t_vec, E_centers)

    n_existing = length(ds["time"])
    idx = n_existing+1 : n_existing+n_new

    ds["time"][idx]              = t_vec
    ds["Ie"][:, :, idx, :]      = Ie_4D
    sync(ds)
end
