using NCDatasets: NCDataset, defDim, defVar


# ======================================================================================== #
#                                  OPTICAL EMISSIONS                                     #
# ======================================================================================== #

"""
    make_volume_excitation_file(directory_to_process)

Read `simulation_data.nc` and `inputs/atmosphere.nc` from `directory_to_process`,
compute volume-excitation-rates for all tracked optical emissions and ionizations,
and write results to `analysis/volume_excitation.nc`.

Returns a [`VolumeExcitationResult`](@ref).
"""
function make_volume_excitation_file(directory_to_process)
    ## Load coordinates
    coord     = load_coordinates(directory_to_process)
    t         = coord.t
    z         = coord.h_atm
    E_centers = coord.E_centers

    ## Load neutral densities
    atm = read_atmosphere_nc(directory_to_process)
    nN2 = atm.nN2
    nO2 = atm.nO2
    nO  = atm.nO

    ## Load emission cross-sections
    σ_4278    = excitation_4278(E_centers)
    σ_6730    = excitation_6730_N2(E_centers)
    σ_7774_O  = excitation_7774_O(E_centers)
    σ_7774_O2 = excitation_7774_O2(E_centers)
    σ_8446_O  = excitation_8446_O(E_centers)
    σ_8446_O2 = excitation_8446_O2(E_centers)
    σ_O1D     = excitation_O1D(E_centers)
    σ_O1S     = excitation_O1S(E_centers)

    ## Load ionization cross-sections
    σ_N2, σ_O2, σ_O = load_cross_sections(E_centers)
    N2_levels, O2_levels, O_levels = load_excitation_threshold()
    σ_Oi  = σ_O'  * O_levels[:, 2]  # equivalent to sum(σ_for_each_reaction * N_ionizations_per_reaction)
    σ_O2i = σ_O2' * O2_levels[:, 2]
    σ_N2i = σ_N2' * N2_levels[:, 2]

    ## Allocate volume-excitation-rate arrays [n_z, n_t]
    n_z = length(z)
    n_t = length(t)
    Q4278    = zeros(n_z, n_t)
    Q6730    = zeros(n_z, n_t)
    Q7774_O  = zeros(n_z, n_t)
    Q7774_O2 = zeros(n_z, n_t)
    Q8446_O  = zeros(n_z, n_t)
    Q8446_O2 = zeros(n_z, n_t)
    QO1D     = zeros(n_z, n_t)  # quenching is not taken into account
    QO1S     = zeros(n_z, n_t)  # quenching is not taken into account
    QOi      = zeros(n_z, n_t)
    QO2i     = zeros(n_z, n_t)
    QN2i     = zeros(n_z, n_t)

    ## Stream Ie over time chunks; each chunk is summed over pitch-angle beams →
    ## omnidirectional flux [n_z, n_t_chunk, n_E] before integrating over energy.
    foreach_Ie_time_chunk(directory_to_process) do Ie_chunk, t_range
        Ie_omni = dropdims(sum(Ie_chunk, dims=2), dims=2)
        tc = t[t_range]
        Q4278[:, t_range]    .= calculate_volume_excitation(z, tc, Ie_omni, σ_4278,    nN2)
        Q6730[:, t_range]    .= calculate_volume_excitation(z, tc, Ie_omni, σ_6730,    nN2)
        Q7774_O[:, t_range]  .= calculate_volume_excitation(z, tc, Ie_omni, σ_7774_O,  nO)
        Q7774_O2[:, t_range] .= calculate_volume_excitation(z, tc, Ie_omni, σ_7774_O2, nO2)
        Q8446_O[:, t_range]  .= calculate_volume_excitation(z, tc, Ie_omni, σ_8446_O,  nO)
        Q8446_O2[:, t_range] .= calculate_volume_excitation(z, tc, Ie_omni, σ_8446_O2, nO2)
        QO1D[:, t_range]     .= calculate_volume_excitation(z, tc, Ie_omni, σ_O1D,     nO)
        QO1S[:, t_range]     .= calculate_volume_excitation(z, tc, Ie_omni, σ_O1S,     nO)
        QOi[:, t_range]      .= calculate_volume_excitation(z, tc, Ie_omni, σ_Oi,      nO)
        QO2i[:, t_range]     .= calculate_volume_excitation(z, tc, Ie_omni, σ_O2i,     nO2)
        QN2i[:, t_range]     .= calculate_volume_excitation(z, tc, Ie_omni, σ_N2i,     nN2)
    end
    Q7774 = Q7774_O .+ Q7774_O2
    Q8446 = Q8446_O .+ Q8446_O2

    ## Write results to analysis/volume_excitation.nc
    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "volume_excitation.nc")

    NCDataset(savefile, "c") do ds
        n_z = length(z)
        n_t = length(t)
        defDim(ds, "altitude", n_z)
        defDim(ds, "time",     n_t)

        alt_v = defVar(ds, "altitude", Float64, ("altitude",);
                       attrib=["units" => "m", "long_name" => "altitude"])
        alt_v[:] = z
        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = t

        for (name, data, long_name) in (
                ("Q4278",    Q4278,    "Volume excitation rate 4278 Å"),
                ("Q6730",    Q6730,    "Volume excitation rate 6730 Å"),
                ("Q7774",    Q7774,    "Volume excitation rate 7774 Å"),
                ("Q7774_O",  Q7774_O,  "Volume excitation rate 7774 Å (O contribution)"),
                ("Q7774_O2", Q7774_O2, "Volume excitation rate 7774 Å (O2 contribution)"),
                ("Q8446",    Q8446,    "Volume excitation rate 8446 Å"),
                ("Q8446_O",  Q8446_O,  "Volume excitation rate 8446 Å (O contribution)"),
                ("Q8446_O2", Q8446_O2, "Volume excitation rate 8446 Å (O2 contribution)"),
                ("QO1D",     QO1D,     "Volume excitation rate O(1D)"),
                ("QO1S",     QO1S,     "Volume excitation rate O(1S)"),
                ("QOi",      QOi,      "O ionization rate"),
                ("QO2i",     QO2i,     "O2 ionization rate"),
                ("QN2i",     QN2i,     "N2 ionization rate"),
            )
            v = defVar(ds, name, Float64, ("altitude", "time");
                       deflatelevel=4,
                       attrib=["units" => "m-3 s-1", "long_name" => long_name])
            v[:, :] = data
        end
    end

    println("Volume excitation rates saved in $savefile")

    return VolumeExcitationResult(
        Q4278, Q6730, Q7774, Q7774_O, Q7774_O2,
        Q8446, Q8446_O, Q8446_O2,
        QO1D, QO1S, QOi, QO2i, QN2i,
        Vector{Float64}(z), Vector{Float64}(t),
        directory_to_process,
    )
end


"""
    calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)

Calculate the volume-excitation-rate for an excitation of interest, produced by the electron
flux `Ie_ztE_omni` that is summed over the beams (omnidirectional).

The excitation of interest is chosen through the cross-section `σ` given to the function.
Note that the neutral density `n` should match the excitation of interest (e.g. use nN2 when
calculating the volume-excitation-rate of the 4278Å optical emission).

# Calling
`Q = calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)`

# Inputs
- `z`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Ie_ztE_omni`: omnidirectional electron flux (#e⁻/m²/s). 3D array [n\\_z, n\\_t, n\\_E]
- `σ`: excitation cross-section (m⁻²). Vector [n\\_E]
- `n`: density of exciteable atmospheric specie (m⁻³). Vector [n\\_z]
"""
function calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)
    n_z = length(z)
    n_t = length(t)
    Q = zeros(n_z, n_t)
    @views for i_t in eachindex(t)
        Q[:, i_t] .= (Ie_ztE_omni[:, i_t, :] * σ) .* n
    end
    return Q
end


"""
    make_column_excitation_file(directory_to_process)

Read `analysis/volume_excitation.nc` from `directory_to_process`, integrate
volume-excitation-rates in altitude (accounting for finite photon travel time),
and write results to `analysis/column_excitation.nc`.

Returns a [`ColumnExcitationResult`](@ref).
"""
function make_column_excitation_file(directory_to_process)
    ve_file = joinpath(directory_to_process, "analysis", "volume_excitation.nc")

    ## Load volume-excitation-rates
    z, t, Q4278, Q6730, Q7774, Q7774_O, Q7774_O2, Q8446, Q8446_O, Q8446_O2, QO1D, QO1S =
        NCDataset(ve_file, "r") do ds
            z      = Array(ds["altitude"])
            t      = Array(ds["time"])
            Q4278  = Array(ds["Q4278"])
            Q6730  = Array(ds["Q6730"])
            Q7774  = Array(ds["Q7774"])
            Q7774_O  = Array(ds["Q7774_O"])
            Q7774_O2 = Array(ds["Q7774_O2"])
            Q8446  = Array(ds["Q8446"])
            Q8446_O  = Array(ds["Q8446_O"])
            Q8446_O2 = Array(ds["Q8446_O2"])
            QO1D   = Array(ds["QO1D"])
            QO1S   = Array(ds["QO1S"])
            (z, t, Q4278, Q6730, Q7774, Q7774_O, Q7774_O2, Q8446, Q8446_O, Q8446_O2, QO1D, QO1S)
        end

    ## Integrate in altitude (with finite-speed-of-light time shift)
    I_4278    = q2colem(t, z, Q4278)
    I_6730    = q2colem(t, z, Q6730)
    I_7774    = q2colem(t, z, Q7774)
    I_7774_O  = q2colem(t, z, Q7774_O)
    I_7774_O2 = q2colem(t, z, Q7774_O2)
    I_8446    = q2colem(t, z, Q8446)
    I_8446_O  = q2colem(t, z, Q8446_O)
    I_8446_O2 = q2colem(t, z, Q8446_O2)
    I_O1D     = q2colem(t, z, QO1D)
    I_O1S     = q2colem(t, z, QO1S)

    ## Write to analysis/column_excitation.nc
    savefile = joinpath(directory_to_process, "analysis", "column_excitation.nc")
    NCDataset(savefile, "c") do ds
        n_t = length(t)
        defDim(ds, "time", n_t)
        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = t

        for (name, data, long_name) in (
                ("I_4278",    I_4278,    "Column excitation rate 4278 Å"),
                ("I_6730",    I_6730,    "Column excitation rate 6730 Å"),
                ("I_7774",    I_7774,    "Column excitation rate 7774 Å"),
                ("I_7774_O",  I_7774_O,  "Column excitation rate 7774 Å (O contribution)"),
                ("I_7774_O2", I_7774_O2, "Column excitation rate 7774 Å (O2 contribution)"),
                ("I_8446",    I_8446,    "Column excitation rate 8446 Å"),
                ("I_8446_O",  I_8446_O,  "Column excitation rate 8446 Å (O contribution)"),
                ("I_8446_O2", I_8446_O2, "Column excitation rate 8446 Å (O2 contribution)"),
                ("I_O1D",     I_O1D,     "Column excitation rate O(1D)"),
                ("I_O1S",     I_O1S,     "Column excitation rate O(1S)"),
            )
            v = defVar(ds, name, Float64, ("time",);
                       deflatelevel=4,
                       attrib=["units" => "m-2 s-1", "long_name" => long_name])
            v[:] = data
        end
    end

    println("Column excitation rates saved in $savefile")

    return ColumnExcitationResult(
        vec(I_4278), vec(I_6730), vec(I_7774), vec(I_7774_O), vec(I_7774_O2),
        vec(I_8446), vec(I_8446_O), vec(I_8446_O2),
        vec(I_O1D), vec(I_O1S), Vector{Float64}(t),
    )
end


using Integrals: SampledIntegralProblem, TrapezoidalRule, solve
using Interpolations: interpolate, extrapolate, Gridded, Linear
"""
    q2colem(t::Vector, z, Q, A = 1, τ = ones(length(z)))

Integrate the volume-excitation-rate (#exc/m³/s) to column-excitation-rate (#exc/m²/s).

Takes into account the time-delay between light emitted at different altitudes. Photons
emitted at at altitude of 200km will arrive at the detector 100e3/3e8 = 0.333 ms later than
electrons emitted at an altitude of 100km. This is a small time-shift, but it is close to
the time-differences corresponding to the phase-shifts between auroral emissions varying at
~10Hz.

The einstein coefficient `A` and effective lifetime `τ` are optional (equal to one by default).

# Calling
`I = q2colem(t, z, Q, A, τ)`

# Inputs
- `z`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Q`: volume-excitation-rate (#exc/m³/s) of the wavelength of interest. 2D array [n\\_z, n\\_t]
- `A`: einstein coefficient (s⁻¹). Scalar (Float or Int)
- `τ`: effective lifetime (s). Vector [n\\_z].

# Output
- `I`: integrated column-excitation-rate (#exc/m²/s) of the wavelength of interest. Vector [n\\_t]
"""
function q2colem(t::Vector, z, Q, A = 1, τ = ones(length(z)))
    c = 2.99792458e8
    Q = Q .* τ .* A
    # Single time step: no time-shift interpolation needed
    if length(t) == 1
        return solve(SampledIntegralProblem(Q, z; dim=1), TrapezoidalRule()).u
    end
    # Create a 2D interpolator over (altitude, time). For each grid point (z[i], t[j]) we
    # will query it at the shifted time `t[j] - (z[i] - z[1]) / c`, which accounts for the
    # finite travel time of photons from altitude z[i] to the bottom of the column. Values
    # outside the interpolation domain are clamped to 0 (photons not yet arrived).
    #
    # Example (c = 1 for illustration):
    #   Q[z, t] (rows = altitude, cols = time):
    #     0.657  0.065  0.313  0.408  0.812
    #     0.780  0.531  0.546  0.575  0.708
    #     0.573  0.556  0.568  0.888  0.728
    #     0.865  0.142  0.085  0.051  0.657
    #     0.708  0.646  0.960  0.932  0.957
    #   After time-shifting → I[z, t]:
    #     0.657  0.065  0.313  0.408  0.812
    #     0.0    0.780  0.531  0.546  0.575   ← z=2 delayed by 1 unit
    #     0.0    0.0    0.573  0.556  0.568   ← z=3 delayed by 2 units
    #     0.0    0.0    0.0    0.865  0.142
    #     0.0    0.0    0.0    0.0    0.708
    # A single 2D interpolator is used (rather than one 1D interpolator per altitude) for
    # convenience and to match the approach from the legacy Matlab code.
    nodes = (z, t)
    itp = interpolate(nodes, Q, Gridded(Linear()))
    itp = extrapolate(itp, 0.0)  # extrapolated values (before arrival) → 0
    I = [itp(z[i], (t[j] - (z[i] - z[1]) / c)) for i in eachindex(z), j in eachindex(t)]
    # Integrate over altitude to get the column-integrated rate
    problem = SampledIntegralProblem(I, z; dim=1)
    method = TrapezoidalRule()
    I_lambda = solve(problem, method)
    return I_lambda.u
end

"""
    q2colem(t::Real, z, Q, A = 1, τ = ones(length(z)))

Same as above, except time is now a scalar (steady-state results). Simple integration in height.
"""
function q2colem(t::Real, z, Q, A = 1, τ = ones(length(z)))
    Q = Q .* τ .* A
    problem = SampledIntegralProblem(Q, z; dim=1)
    method = TrapezoidalRule()
    I_lambda = solve(problem, method)
    return I_lambda.u
end


# ======================================================================================== #
#                         AuroraSimulation convenience wrappers                          #
# ======================================================================================== #

"""
    make_volume_excitation_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_volume_excitation_file`](@ref) on `sim.output.savedir`.
"""
make_volume_excitation_file(sim::AuroraSimulation) = make_volume_excitation_file(sim.output.savedir)

"""
    make_column_excitation_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_column_excitation_file`](@ref) on `sim.output.savedir`.
"""
make_column_excitation_file(sim::AuroraSimulation) = make_column_excitation_file(sim.output.savedir)
