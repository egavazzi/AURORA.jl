✅
❌

#=
STRUCTURE:

Part where parameters are set ❌
Make setup ✅
Loop over n_loop
    Load incoming flux ✅
    Then loop over E
        Make matrix A ✅
        Make matrix B ✅
        Make matrix D ✅
        Call CN scheme ✅
        
        Then it is energy degradation time:
            Calculate e-e collisions ✅
            Calculate inelastic collisions ✅
            Calculate ionizations collisions, with secondary e- and scattering of primary ✅
=#

## This is the control script from where simulations are run
using Aurora
using MAT

## Setting parameters
altitude_max = 400;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 7231;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith 

t = 0:0.001:0.05;           # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile


## Get atmosphere
println("Calling Matlab for the setup...")
h_atm, n_neutrals, ne, Te, E, dE, 
    E_levels_neutrals, σ_neutrals, secondary_e,
    θ_lims, μ_lims, μ_center, μ_scatterings = setup(altitude_max, θ_lims, E_max);


# Initialise variables
Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));


# Load incoming flux
input_file = "/mnt/data/etienne/AURORA/MI_coupling/TIME2PLAY/conversion_1.27e7-1/Ie_incoming.mat"
Ie_top = Ie_top_from_file(input_file, μ_center, t, E, n_loop);


# Create the folder to save the data to
savedir = string(pkgdir(Aurora, "data"), "/",
                        Dates.format(now(), "yyyymmdd-HHMM"))
println("Results will be save in ", savedir)


# Looping over n_loop
for i in 1:n_loop
    println("Starting the calculations for ")

    D = make_D(E, dE, θ_lims);
    # Extract the top flux for the current loop
    Ie_top_local = Ie_top[:, (1 + (i - 1) * (length(t) - 1)) : (length(t) + (i - 1) * (length(t) - 1)), :]

    # Looping over energy
    for iE in length(E):-1:1
        A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);
        
        # Make a finer θ for the scattering calculations
        finer_θ = range(0, π, length=721); 
        # Calculate the phase functions for E[iE] and put them in a Tuple
        phaseN2e, phaseN2i = phase_fcn_N2(finer_θ, E[iE]);
        phaseO2e, phaseO2i = phase_fcn_O2(finer_θ, E[iE]);
        phaseOe, phaseOi = phase_fcn_O(finer_θ, E[iE]);
        phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));

        B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals, 
                                            phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                            μ_scatterings.BeamWeight_relative, finer_θ);

        # Compute the flux of e-
        Ie[:, :, iE] = Crank_Nicolson(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]), 
                                        A, B, D[iE, :], Q[:, :, iE], Ie_top_local[:, :, iE], I0[:, iE]);

        # Update the cascading of e-
        cascading_neutrals = (cascading_N2, cascading_O2, cascading_O)
        update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals, 
                    cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight_discrete, μ_center)
    end

    # Update the starting e- flux profile
    I0 = Ie[:, end, :]

    # Save e- flux
    savefile = string(pkgdir(Aurora, "data"), "/",
                        Dates.format(now())

    file = matopen(filename)
    Ie_top_raw = read(file, "Ie_total")
    close(file

end