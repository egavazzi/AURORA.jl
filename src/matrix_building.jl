using SparseArrays: SparseArrays, spdiagm

"""
    TransportMatrices

Container for the matrices used in the electron transport equations.

# Fields
- `A::Vector{Float64}`: Loss term matrix (altitude dimension)
- `B::Array{Float64, 3}`: Scattering term matrix (altitude x angle x angle)
- `D::Array{Float64, 2}`: Diffusion coefficient matrix (energy x angle)
- `Q::Array{Float64, 3}`: Source term array (altitude*angle x time x energy)
- `Ddiffusion::SparseArrays.SparseMatrixCSC{Float64, Int64}`: Diffusion operator (altitude x altitude)

The struct is mutable to allow efficient in-place updates of Q during the energy cascade.
"""
mutable struct TransportMatrices
    A::Vector{Float64}
    B::Array{Float64, 3}
    D::Array{Float64, 2}
    Q::Array{Float64, 3}
    Ddiffusion::SparseArrays.SparseMatrixCSC{Float64, Int64}
end

"""
    TransportMatrices(n_altitude, n_angle, n_time, n_energy)

Construct an empty TransportMatrices container with zeros.

# Arguments
- `n_altitude::Int`: Number of altitude grid points
- `n_angle::Int`: Number of angle bins
- `n_time::Int`: Number of time steps
- `n_energy::Int`: Number of energy bins
"""
function TransportMatrices(n_altitude::Int, n_angle::Int, n_time::Int, n_energy::Int)
    A = zeros(Float64, n_altitude)
    B = zeros(Float64, n_altitude, n_angle, n_angle)
    D = zeros(Float64, n_energy, n_angle)
    Q = zeros(Float64, n_altitude * n_angle, n_time, n_energy)
    Ddiffusion = SparseArrays.spzeros(Float64, n_altitude, n_altitude)
    return TransportMatrices(A, B, D, Q, Ddiffusion)
end

"""
    loss_to_thermal_electrons(E, ne, Te)

Suprathermal electron energy loss function due to electron-electron collisions.

This function calculates the electron energy loss function due to electron-electron
interaction. It uses the analytic form given for the energy-transfer rate from
photoelectrons (or suprathermal electrons) to thermal electrons, given by Swartz
and Nisbet (1971). The expression fits the classical formulation of Itikawa and
Aono (1966) at low energies and gives a smooth transition to fit the quantum
mechanical equation of Schunk and Hays (1971).

# Arguments
- `E::Real`: Energy level [eV]. Scalar value.
- `ne::Vector`: Ambient electron concentration [/m³], length nZ.
- `Te::Vector`: Electron temperature [K], length nZ.

# Returns
- `Le::Vector`: Electron energy loss function [eV/m], length nZ.

# Notes
The paper by Swartz and Nisbet uses electron density in cm⁻³; here the constant
is rescaled to use m⁻³ instead. We calculate the loss function dE/ds(E,ne,Te)
directly and not as in Swartz and Nisbet dE/ds(E,ne,Te)/ne.

The loss is set to zero when the suprathermal electron energy E is below the
thermal electron energy Ee = kB*Te/qₑ.

# References
- Swartz, W. E., J. S. Nisbet, and A. E. S. Green (1971), Analytic expression
  for the energy transfer rate from photoelectrons to thermal electrons,
  J. Geophys. Res., 76(34), 8425-8426, doi: 10.1029/JA076i034p08425.
- Itikawa, Y., and O. Aono (1966), Energy change of a charged particle moving
  in a plasma, Phys. Fluids, 9, 1259-1261.
- Schunk, R. W., and P. B. Hays (1971), Photoelectron energy losses to thermal
  electrons, Planet. Space Sci., 19, 113-117.
"""
function loss_to_thermal_electrons(E::Real, nₑ, Tₑ)
    kB = 1.380662e-23     # Boltzmann constant [J/K]
    qₑ = 1.6021773e-19    # elementary charge [C]

    # Thermal electron energy in eV
    Eₑ = kB / qₑ * Tₑ

    # Calculate energy loss function using Swartz & Nisbet (1971) formula
    Le = 3.0271e-10 * nₑ .^ 0.97 ./ E^0.44 .* ((E .- Eₑ) ./ (E .- 0.53 * Eₑ)) .^ 2.36  / v_of_E(E)

    # Set loss to zero when E < Ee (below thermal energy)
    Le[E .< Eₑ] .= 0

    return Le
end

# Depreciated function, for demo
function beams2beams_demo(phase_fcn, Pmu2mup, BeamWeight_relative)
    B2B = zeros(size(Pmu2mup, 3),size(Pmu2mup, 3));
    for i = size(Pmu2mup, 3):-1:1
        B2B[i, :] = BeamWeight_relative * (@view(Pmu2mup[:, :, i]) * phase_fcn);
    end
    return B2B
end

# The new functions, for faster calculations
function prepare_beams2beams(BeamWeight_relative, Pmu2mup)
    B2B_fragment = zeros(size(BeamWeight_relative, 1), size(Pmu2mup, 2), size(Pmu2mup, 3))
    for i = size(Pmu2mup, 3):-1:1
        B2B_fragment[:, :, i] = BeamWeight_relative * (@view(Pmu2mup[:, :, i]));
    end
    return B2B_fragment
end

function beams2beams(phase_fcn, B2B_fragment)
    B2B = zeros(size(B2B_fragment, 3), size(B2B_fragment, 3));
    for i = size(B2B_fragment, 3):-1:1
        B2B[i, :] = @view(B2B_fragment[:, :, i]) * phase_fcn;
    end
    return B2B
end

## ----------------------------------------------------- ##

function update_A!(A, n_neutrals, σ_neutrals, ne, Te, E, dE, iE)
    fill!(A, 0.0)
    # Loop over the neutral species
    for i1 in 1:length(n_neutrals)
        n = n_neutrals[i1];  # Neutral density
        σ = σ_neutrals[i1];  # Array with collision cross sections

        # add elastic collisions
        A .+= n .* σ[1, iE];

        # add inelastic and ionization collisions
        for i2 in 2:size(σ, 1)      # Loop over the different collisions, because
            A .+= n .* σ[i2, iE];  # they have different cross sections
        end
    end

    # add losses due to electron-electron collisions
    A .+= loss_to_thermal_electrons(E[iE] + dE[iE] / 2, ne, Te) ./ dE[iE];

    return nothing
end

function update_B!(B, n_neutrals, σ_neutrals, E_levels_neutrals, phase_fcn_neutrals, dE, iE, B2B_fragment, finer_θ)
    # Zero out B in place
    fill!(B, 0.0)
    B2B_inelastic_neutrals = Vector{Matrix{Float64}}(undef, length(n_neutrals));
    # Loop over the neutral species
    for i in 1:length(n_neutrals)
        n = n_neutrals[i];                  # Neutral density
        σ = σ_neutrals[i];                  # Array with collision cross sections
        E_levels = E_levels_neutrals[i];    # Array with collision enery levels and number of secondary e-
        phase_fcn = phase_fcn_neutrals[i];   # Tuple with two phase function arrays, the first for elastic collisions
                                                    # and the second for inelastic collisions

        # Convert to 3D the scattering probabilities that are in 1D
        phase_fcn_e = convert_phase_fcn_to_3D(phase_fcn[1][:, iE], finer_θ);
        phase_fcn_i = convert_phase_fcn_to_3D(phase_fcn[2][:, iE], finer_θ);
        B2B_elastic = beams2beams(phase_fcn_e, B2B_fragment);
        B2B_inelastic = beams2beams(phase_fcn_i, B2B_fragment);

        # add scattering from elastic collisions
        for i1 in axes(B2B_elastic, 1)
            for i2 in axes(B2B_elastic, 2)
                B[:, i1, i2] .= @view(B[:, i1, i2]) .+ n .* σ[1, iE] .* B2B_elastic[i1, i2];
            end
        end

        # add scattering from inelastic and ionization collisions
        for i1 in 2:size(σ, 1)
            for i2 in axes(B2B_inelastic, 1)
                for i3 in axes(B2B_inelastic, 2)
                    # The last factor corrects for the case where the energy loss
                    # E_levels[i1, 1] is smaller than the width in energy of the energy bin.
                    # That is, when dE[iE] > E_levels[i1,1], only the fraction
                    # E_levels[i1,1] / dE is lost from the energy bin [E[iE], E[iE] + dE[iE]].
                    B[:, i2, i3] .= @view(B[:, i2, i3]) .+ n .* σ[i1, iE] .* B2B_inelastic[i2, i3] .*
                                                    max(0, 1 - E_levels[i1, 1] ./ dE[iE]);
                end
            end
        end

        # Save the inelastic B2B matrices for the future energy degradations (updates of Q)
        B2B_inelastic_neutrals[i] = copy(B2B_inelastic);
    end
    return B2B_inelastic_neutrals
end

function update_D!(D, E, dE, θ_lims)
    θ_lims_rad = deg2rad.(θ_lims)
    nE = 3
    nθ = 3
    # n_ti = 701
    # n_thi = 401
    for iE in length(E):-1:1
        v = range(v_of_E(E[iE]), v_of_E(E[iE] + dE[iE]), length=nE)
        for iθ in 1:(length(θ_lims_rad) - 1)
            θa = θ_lims_rad[iθ]
            θb = θ_lims_rad[iθ + 1]
            if θ_lims_rad[iθ] == π/2
                θa = θ_lims_rad[iθ] * 0.8 + 0.2 * θ_lims_rad[iθ + 1]
            end
            if θ_lims_rad[iθ + 1] == π/2
                θb = θ_lims_rad[iθ] * 0.2 + 0.8 * θ_lims_rad[iθ + 1]
            end
            θ = range(θa, θb, length=nθ)
            # θ4i = range(minimum(θ), maximum(θ), n_thi)
            v_par = [A * cos(B) for A in v, B in θ]
            t_arrival = 500e3 ./ v_par
            at_a = (maximum(t_arrival) + minimum(t_arrival)) / 2
            dt_a = (maximum(t_arrival) - minimum(t_arrival))
            D_val = (dt_a / 4)^2 / at_a

            D[iE, iθ] = abs(D_val)
        end
    end
    return nothing
end

function update_Ddiffusion!(Ddiffusion, z)
    dzd = z[2:end-1] - z[1:end-2]
    dzu = z[3:end]   - z[2:end-1]

    dsup  = [2 ./ (dzd .* (dzd + dzu)) ; 0]
    dMain = [0 ; -2 ./ (dzd .* dzu) ; 0]
    dsub  = [0 ; 2 ./ (dzu .* (dzd + dzu))]

    # Update the sparse matrix in place by rebuilding it
    # Note: SparseArrays don't support true in-place modification of structure
    D2M = spdiagm( -1 => dsup,
                    0 => dMain,
                    1 => dsub)

    # Copy the new matrix structure into the existing one
    copyto!(Ddiffusion, D2M)

    return nothing
end

"""
    update_matrices!(matrices, n_neutrals, σ_neutrals, ne, Te, E_levels_neutrals,
                     phase_fcn_neutrals, E, dE, iE, B2B_fragment, finer_θ)

Update the A and B matrices in place for a given energy level iE.

# Arguments
- `matrices::TransportMatrices`: Container to update
- `n_neutrals, σ_neutrals, ne, Te, E_levels_neutrals, phase_fcn_neutrals`: Atmosphere and cross section data
- `E, dE, iE`: Energy grid and current energy index
- `B2B_fragment, finer_θ`: Pre-computed beam-to-beam fragments and angle grid

# Returns
- `B2B_inelastic_neutrals`: Array of inelastic beam-to-beam matrices for cascading calculations
"""
function update_matrices!(matrices::TransportMatrices, n_neutrals, σ_neutrals, ne, Te,
                         E_levels_neutrals, phase_fcn_neutrals, E, dE, iE, B2B_fragment, finer_θ)
    # Update A matrix
    update_A!(matrices.A, n_neutrals, σ_neutrals, ne, Te, E, dE, iE)

    # Update B matrix and get B2B_inelastic
    B2B_inelastic_neutrals = update_B!(matrices.B, n_neutrals, σ_neutrals, E_levels_neutrals,
                                                 phase_fcn_neutrals, dE, iE, B2B_fragment, finer_θ)

    return B2B_inelastic_neutrals
end


"""
    initialize_transport_matrices(h_atm, μ_center, t, E, dE, θ_lims)

Create a TransportMatrices container initialized with zeros for A, B, D, Q and Ddiffusion.

# Arguments
- `h_atm`: Altitude grid
- `μ_center`: Cosine of angle centers
- `t`: Time grid
- `E, dE`: Energy grid and bin widths
- `θ_lims`: Angle bin limits

# Returns
- `matrices::TransportMatrices`: Initialized container for the transport matrices
"""
function initialize_transport_matrices(h_atm, μ_center, t, E, dE, θ_lims)
    n_altitude = length(h_atm)
    n_angle = length(μ_center)
    n_time = length(t)
    n_energy = length(E)

    matrices = TransportMatrices(n_altitude, n_angle, n_time, n_energy)

    return matrices
end
