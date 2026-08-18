using LinearAlgebra: mul!
using LoopVectorization: @tturbo

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
    Le = similar(nₑ, Float64)
    return loss_to_thermal_electrons!(Le, E, nₑ, Tₑ)
end

# Calculate energy loss function using Swartz & Nisbet (1971) formula
function loss_to_thermal_electrons!(Le, E::Real, nₑ, Tₑ)
    @assert axes(Le) == axes(nₑ) == axes(Tₑ)

    velocity = v_of_E(E)
    energy_factor = 3.0271e-10 / (E^0.44 * velocity)

    for i in eachindex(Le)
        Eₑ = kB / qₑ * Tₑ[i] # thermal electron energy in eV
        if E < Eₑ
            Le[i] = 0.0
        else
            ratio = (E - Eₑ) / (E - 0.53 * Eₑ)
            Le[i] = energy_factor * nₑ[i]^0.97 * ratio^2.36
        end
    end

    return Le
end

# Beam-to-beam scattering helpers
# - `B2B_kernel` precomputes the contribution of each scattering-grid direction
#   between every pair of pitch-angle beams.
# - `B2B` is computed by applying a species phase function to this kernel. It contains
#   the final scattering weights between pairs of pitch-angle beams.
function beams2beams_kernel(Ω_subbeam_relative, P_scatter)
    B2B_kernel = zeros(size(Ω_subbeam_relative, 1), size(P_scatter, 2), size(P_scatter, 3))
    for i = size(P_scatter, 3):-1:1
        B2B_kernel[:, :, i] = Ω_subbeam_relative * (@view(P_scatter[:, :, i]));
    end
    return B2B_kernel
end

function beams2beams(phase_fcn, B2B_kernel)
    B2B = zeros(size(B2B_kernel, 3), size(B2B_kernel, 3));
    for i = size(B2B_kernel, 3):-1:1
        B2B[i, :] = @view(B2B_kernel[:, :, i]) * phase_fcn;
    end
    return B2B
end

function beams2beams!(B2B, phase_fcn, B2B_kernel)
    for i = size(B2B_kernel, 3):-1:1
        @views mul!(B2B[i, :], B2B_kernel[:, :, i], phase_fcn)
    end
    return B2B
end

function update_A!(matrices::TransportMatrices, model::AuroraModel, iE)
    A = matrices.A
    Le = matrices.Le
    ionosphere = model.ionosphere
    energy_grid = model.energy_grid
    E_centers = energy_grid.E_centers
    ΔE = energy_grid.ΔE

    fill!(A, 0.0)
    # Loop over the neutral species
    for sp in model.species
        n = sp.density
        σ = sp.cross_sections

        # add elastic collisions
        A .+= n .* σ[1, iE];

        # add inelastic and ionization collisions
        for i2 in 2:size(σ, 1)
            A .+= n .* σ[i2, iE];
        end
    end

    # add losses due to electron-electron collisions
    loss_to_thermal_electrons!(Le, E_centers[iE], ionosphere.ne, ionosphere.Te)
    A .+= Le ./ ΔE[iE];

    return nothing
end

function update_B!(matrices::TransportMatrices, model::AuroraModel, iE, B2B_kernel)
    B = matrices.B
    energy_grid = model.energy_grid
    scattering = model.scattering
    ΔE = energy_grid.ΔE
    finer_θ = scattering.θ_scatter

    # Zero out B in place
    fill!(B, 0.0)
    # Pre-allocated per-species inelastic beam-to-beam buffers
    B2B_inelastic_neutrals = matrices.B2B_inelastic_neutrals
    # Loop over the neutral species
    for (i, sp) in enumerate(model.species)
        n = sp.density
        σ = sp.cross_sections
        E_levels = sp.excitation_levels
        phase_fcn = sp.phase_fcn

        # Convert to 3D the scattering probabilities that are in 1D
        convert_phase_fcn_to_3D!(matrices.phase_fcn_e, @view(phase_fcn[1][:, iE]), finer_θ);
        convert_phase_fcn_to_3D!(matrices.phase_fcn_i, @view(phase_fcn[2][:, iE]), finer_θ);
        B2B_elastic = matrices.B2B_elastic
        B2B_inelastic = B2B_inelastic_neutrals[i]
        beams2beams!(B2B_elastic, matrices.phase_fcn_e, B2B_kernel);
        beams2beams!(B2B_inelastic, matrices.phase_fcn_i, B2B_kernel);

        # add scattering from elastic collisions
        n_z = length(n)
        n_angle = size(B2B_elastic, 1)
        σ_elastic = σ[1, iE]

        @tturbo for i2 in 1:n_angle
            for i1 in 1:n_angle
                for iz in 1:n_z
                    B[iz, i1, i2] += n[iz] * σ_elastic * B2B_elastic[i1, i2]
                end
            end
        end

        # add scattering from inelastic and ionization collisions
        n_collisions = size(σ, 1)
        for i_coll in 2:n_collisions
            σ_coll = σ[i_coll, iE]
            E_loss = E_levels[i_coll, 1]
            # The last factor corrects for the case where the energy loss
            # E_levels[i_coll, 1] is smaller than the width in energy of the energy bin.
            # That is, when dE[iE] > E_levels[i_coll,1], only the fraction
            # E_levels[i_coll,1] / dE is lost from the energy bin [E[iE], E[iE] + dE[iE]].
            correction_factor = max(0.0, 1.0 - E_loss / ΔE[iE])

            @tturbo for i3 in 1:n_angle
                for i2 in 1:n_angle
                    for iz in 1:n_z
                        B[iz, i2, i3] += n[iz] * σ_coll * B2B_inelastic[i2, i3] * correction_factor
                    end
                end
            end
        end
    end
    return B2B_inelastic_neutrals
end

"""
    update_matrices!(matrices, model, iE, B2B_kernel)

Update the A and B matrices in place for a given energy level iE.

# Arguments
- `matrices::TransportMatrices`: Container to update
- `model`: `AuroraModel` (grids + atmosphere + physics)
- `iE`: Current energy index
- `B2B_kernel`: Pre-computed beam-to-beam scattering kernel

# Returns
- `B2B_inelastic_neutrals`: Array of inelastic beam-to-beam matrices for cascading calculations
"""
function update_matrices!(matrices::TransportMatrices, model::AuroraModel, iE, B2B_kernel)
    update_A!(matrices, model, iE)
    return update_B!(matrices, model, iE, B2B_kernel)
end


"""
    initialize_transport_matrices(model, t)

Create a `TransportMatrices` container initialized with zeros for A, B and Q.
"""
function initialize_transport_matrices(model::AuroraModel, t)
    n_altitude = length(model.altitude_grid.h)
    n_angle = length(model.pitch_angle_grid.μ_center)
    n_time = length(t)
    n_energy = model.energy_grid.n

    matrices = TransportMatrices(n_altitude, n_angle, n_time, n_energy)

    # Size the scratch buffers that depend on the number of species and the scattering grid.
    n_species = length(model.species)
    n_θ = length(model.scattering.θ_scatter)
    matrices.B2B_inelastic_neutrals = [zeros(Float64, n_angle, n_angle) for _ in 1:n_species]
    matrices.phase_fcn_e = zeros(Float64, n_θ)
    matrices.phase_fcn_i = zeros(Float64, n_θ)

    return matrices
end
