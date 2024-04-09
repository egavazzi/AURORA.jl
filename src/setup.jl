using MATLAB

"""
    setup(path_to_AURORA_matlab, top_altitude, θ_lims, E_max)

Load the atmosphere, the energy grid, the collision cross-sections, ... \\
It calls a lot of functions from the original MATLAB code.

# Calling
`h_atm, ne, Te, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, θ_lims, μ_lims, μ_center,
μ_scatterings = setup(altitude_max, θ_lims, E_max);`

# Inputs
- `path_to_AURORA_matlab`: path to where the original Matlab AURORA package is installed
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up. Vector [n_beam]
- `E_max`: upper limit for the energy grid (in eV)

# Outputs
- `h_atm`: altitude (m), vector [nZ]
- `ne`: e- densities (m³), vector [nZ]
- `Te`: e- temperatures (K), vector [nZ]
- `E`: energy grid (eV), vector [nE]
- `dE`: energy bin sizes(eV), vector [nE]
- `n_neutrals`: neutral densities (m³), Tuple of vectors ([nZ], ..., [nZ])
- `E_levels_neutrals`: collisions energy levels and number of secondary e- produced,
    Tuple of matrices ([n`_`levels x 2], ..., [n`_`levels x 2])
- `σ_neutrals`: collision cross-sections (m³), Tuple of matrices ([n`_`levels x nE], ...,
    [n`_`levels x nE])
- `θ_lims`: pitch angle limits of the e- beams (deg), vector [n_beam + 1]
- `μ_lims`: cosine of the pitch angle limits of the e- beams, vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams, vector [n_beam]
- `μ_scatterings`: Tuple with several of the scattering informations, namely
    μ`_`scatterings = `(Pmu2mup, BeamWeight_relative, BeamWeight)`
    + `Pmu2mup`: probabilities for scattering in 3D from beam to beam, matrix [721x721]
    + `BeamWeight_relative`: relative contribution from within each beam, matrix [18 x
        n_beam]
    + `BeamWeight`: solid angle for each stream (ster), vector [n_beam]
"""
function setup(path_to_AURORA_matlab, top_altitude, θ_lims, E_max)
    ## Creating a MATLAB session
    s1 = MSession();
    @mput path_to_AURORA_matlab
    # mat"
    # addpath('/mnt/data/etienne/AURORA','-end')
    # add_AURORA
    # cd /mnt/data/etienne/AURORA/
    # "
    mat"
    addpath(path_to_AURORA_matlab,'-end')
    add_AURORA
    cd(path_to_AURORA_matlab)
    "
    ## Loading atmosphere
    @mput top_altitude
    mat"
    load msis20051008_3.dat

    dz = @(n) 150 + 150/200*(0:(n-1))' +1.2*exp(((0:(n-1))-150)/17)';
    z_atm = 100e3 + cumsum(dz(331)) - dz(1);
    [~, i_zmax] = min(abs(z_atm - double(top_altitude)*1e3));
    h_atm = z_atm(1:i_zmax);

    OPS.atmosphere = interp1(msis20051008_3(:,1),msis20051008_3,h_atm/1e3,'pchip');
    OPS.atmosphere(:,2:4) = OPS.atmosphere(:,2:4)*1e6; % into m^-3
    OPS.atmosphere(:,5) = OPS.atmosphere(:,5)*1e3;     % into kg/m^-3
    nO = OPS.atmosphere(:,2);
    nN2 = OPS.atmosphere(:,3);
    nO2 = OPS.atmosphere(:,4);
    nO(end-2:end)  = 0;
    nO2(end-2:end) = 0;
    nN2(end-2:end) = 0;
    nO(end-5:end-3)  = (erf((1:-1:-1)'/2)+1)/2.*nO(end-5:end-3);
    nO2(end-5:end-3) = (erf((1:-1:-1)'/2)+1)/2.*nO2(end-5:end-3);
    nN2(end-5:end-3) = (erf((1:-1:-1)'/2)+1)/2.*nN2(end-5:end-3);
    "
    h_atm = @mget h_atm;
    nN2 = @mget nN2;
    nO2 = @mget nO2;
    nO = @mget nO;
    n_neutrals = (nN2 = nN2, nO2 = nO2, nO = nO);

    ## Loading n_e(z) and Te as estimated from EISCAT observations
    mat"
    load iri20051008.dat
    Iri20051008 = iri20051008(1:end,:);
    ne = interp1(Iri20051008(:,1),Iri20051008(:,2),h_atm/1e3,'pchip')*10;
    ne(h_atm<82e3) = 1;
    Te = interp1(Iri20051008(:,1),Iri20051008(:,6),h_atm/1e3,'pchip');
    Te(Te==-1) = 350;
    "
    ne = @mget ne;
    Te = @mget Te;

    ## Excitation thresholds
    mat"
    load N2_levels.dat
    load O2_levels.dat
    load O_levels.dat
    "
    N2_levels = @mget N2_levels;
    O2_levels = @mget O2_levels;
    O_levels = @mget O_levels;
    E_levels_neutrals = (N2_levels = N2_levels, O2_levels = O2_levels, O_levels = O_levels);

    ## Energy grid
    mat"
    dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;
    E = cumsum(dEfcn(0:1000,0.15,11.5,0.05,80))+1.9;
    dE = diff(E);dE = [dE,dE(end)];
    "
    E = vec(@mget E);
    dE = vec(@mget dE);

    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    dE = dE[1:iE_max];                      # crop dE accordingly

    ## Collision cross-sections
    mat"
    [XsO,xs_fcnO] = get_all_xs('O',E+dE/2);
    [XsO2,xs_fcnO2] = get_all_xs('O2',E+dE/2);
    [XsN2,xs_fcnN2] = get_all_xs('N2',E+dE/2);
    "
    σ_N2 = @mget XsN2;
    σ_O2 = @mget XsO2;
    σ_O = @mget XsO;
    σ_neutrals = (σ_N2 = σ_N2, σ_O2 = σ_O2, σ_O = σ_O);

    ## X-streams beam-to-beam calculations
    theta_lims2do = reshape(Vector(θ_lims), 1, :);
    @mput theta_lims2do
    mat"
    theta_lims2do = double(theta_lims2do)
    [Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory);
    mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
    c_o_mu = mu_avg(mu_lims);
    "
    θ_lims = vec(@mget theta_lims2do);
    μ_lims = vec(@mget mu_lims);
    μ_center = vec(@mget c_o_mu);

    Pmu2mup = @mget Pmu2mup; # Probability mu to mu prime
    BeamWeight_relative = @mget theta2beamW;

    # This beam weight is calculated in a continuous way
    BeamWeight = 2π .* vec(@mget BeamW);
    # Here we normalize BeamWeight_relative, as it is supposed to be a relative weighting matrix with the relative
    # contribution from withing each beam. It means that when summing up along each beam, we should get 1
    BeamWeight_relative = BeamWeight_relative ./ repeat(sum(BeamWeight_relative, dims=2), 1, size(BeamWeight_relative, 2));



    μ_scatterings = (Pmu2mup = Pmu2mup, BeamWeight_relative = BeamWeight_relative, BeamWeight = BeamWeight);

    ## Closing the MATLAB session
    close(s1)

    return h_atm, ne, Te, E, dE,
        n_neutrals, E_levels_neutrals, σ_neutrals,
        θ_lims, μ_lims, μ_center, μ_scatterings
end


"""
    setup_new(top_altitude, θ_lims, E_max, msis_file, iri_file)

Load the atmosphere, the energy grid, the collision cross-sections, ... \\
It is a partial rework of the `setup` function. Some parts are implemented in pure Julia,
but there are still some calls to the original MATLAB code.

# Calling
`h_atm, ne, Te, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, θ_lims, μ_lims, μ_center,
μ_scatterings = setup_new(top_altitude, θ_lims, E_max, msis_file, iri_file);`

# Inputs
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up. Vector [n_beam]
- `E_max`: upper limit for the energy grid (in eV)
- `msis_file`: path to the msis file to use
- `iri_file`: path to the iri file to use

# Outputs
- `h_atm`: altitude (m), vector [nZ]
- `ne`: e- densities (m³), vector [nZ]
- `Te`: e- temperatures (K), vector [nZ]
- `E`: energy grid (eV), vector [nE]
- `dE`: energy bin sizes(eV), vector [nE]
- `n_neutrals`: neutral densities (m³), Tuple of vectors ([nZ], ..., [nZ])
- `E_levels_neutrals`: collisions energy levels and number of secondary e- produced,
    Tuple of matrices ([n`_`levels x 2], ..., [n`_`levels x 2])
- `σ_neutrals`: collision cross-sections (m³), Tuple of matrices ([n`_`levels x nE], ...,
    [n`_`levels x nE])
- `θ_lims`: pitch angle limits of the e- beams (deg), vector [n_beam + 1]
- `μ_lims`: cosine of the pitch angle limits of the e- beams, vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams, vector [n_beam]
- `μ_scatterings`: Tuple with several of the scattering informations, namely
    μ`_`scatterings = `(Pmu2mup, BeamWeight_relative, BeamWeight)`
    + `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x
    n`_`direction]
    + `BeamWeight_relative`: relative contribution from within each beam. Matrix [n`_`beam x
    n`_`direction]
    + `BeamWeight`: solid angle for each stream (ster). Vector [n_beam]
    + `theta1`: scattering angles used in the calculations. Vector [n_direction]
"""
function setup_new(top_altitude, θ_lims, E_max, msis_file, iri_file)
    h_atm = make_altitude_grid(top_altitude)
    E, dE = make_energy_grid(E_max)
    μ_lims, μ_center, μ_scatterings = make_scattering_matrices(θ_lims)
    n_neutrals = load_neutral_densities(msis_file, h_atm)
    ne, Te = load_electron_properties(iri_file, h_atm)
    E_levels_neutrals = load_excitation_threshold()

    ## Collision cross-sections
    print("Calling Matlab for the cross-sections...")
    # Translating all of this to Julia is a big task. So for now we still use the
    # Matlab code.
    # TODO: translate that part
    # Creating a MATLAB session
    path_to_AURORA_matlab = pkgdir(AURORA, "MATLAB_dependencies")
    s1 = MSession();
    @mput path_to_AURORA_matlab
    mat"
    addpath(genpath(path_to_AURORA_matlab))
    cd(path_to_AURORA_matlab)
    "
    @mput E
    @mput dE
    mat"
    E = E';
    dE = dE';
    [XsO,xs_fcnO] = get_all_xs('O',E+dE/2);
    [XsO2,xs_fcnO2] = get_all_xs('O2',E+dE/2);
    [XsN2,xs_fcnN2] = get_all_xs('N2',E+dE/2);
    "
    σ_N2 = @mget XsN2;
    σ_O2 = @mget XsO2;
    σ_O = @mget XsO;
    σ_neutrals = (σ_N2 = σ_N2, σ_O2 = σ_O2, σ_O = σ_O);
    # Closing the MATLAB session
    close(s1)
    println(" done.")

    return h_atm, ne, Te, E, dE,
    n_neutrals, E_levels_neutrals, σ_neutrals,
    μ_lims, μ_center, μ_scatterings
end


"""
    make_altitude_grid(top_altitude)

Create an altitude grid based on the `top_altitude` given as input.

# Calling
`h_atm = make_altitude_grid(top_altitude)`

# Inputs
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation

# Outputs
- `h_atm`: altitude (m), vector [nZ]
"""
function make_altitude_grid(top_altitude)
    Δz(n) = 150 .+
            150 / 200 * (0:(n - 1)) .+
            1.2 * exp.(Complex.(((0:(n - 1)) .- 150) / 22) .^ .9)
    h_atm = 100e3 .+ cumsum(real.(Δz(450))) .- real.(Δz(1))
    i_zmax = findmin(abs.(h_atm .- top_altitude * 1e3))[2]
    h_atm = h_atm[1:i_zmax]
    return h_atm
end

"""
    make_energy_grid(E_max)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`E, dE = make_energy_grid(E_max)`

# Inputs
- `E_max`: upper limit for the energy grid (in eV)

# Outputs
- `E`: energy grid (eV), vector [nE]
- `dE`: energy bin sizes(eV), vector [nE]
"""
function make_energy_grid(E_max)
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:2000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    dE = diff(E); dE = [dE; dE[end]]
    return E, dE
end

"""
    make_scattering_matrices(θ_lims)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`μ_lims, μ_center, μ_scatterings = make_scattering_matrices(θ_lims)`

# Inputs
- `θ_lims`: pitch angle limits of the e- beams (deg), vector [n_beam + 1]

# Outputs
- `μ_lims`: cosine of the pitch angle limits of the e- beams, vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams, vector [n_beam]
- `μ_scatterings`: Tuple with several of the scattering informations, namely
    μ`_`scatterings = `(Pmu2mup, BeamWeight_relative, BeamWeight)`
    + `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x
    n`_`direction]
    + `BeamWeight_relative`: relative contribution from within each beam. Matrix [n`_`beam x
    n`_`direction]
    + `BeamWeight`: solid angle for each stream (ster). Vector [n_beam]
    + `theta1`: scattering angles used in the calculations. Vector [n_direction]
"""
function make_scattering_matrices(θ_lims)
    μ_lims = cosd.(θ_lims);
    μ_center = mu_avg(θ_lims);
    BeamWeight = beam_weight(θ_lims); # this beam weight is calculated in a continuous way
    Pmu2mup, _, BeamWeight_relative, θ₁ = load_scattering_matrices(θ_lims, 720)
    μ_scatterings = (Pmu2mup = Pmu2mup, BeamWeight_relative = BeamWeight_relative, BeamWeight = BeamWeight, theta1 = θ₁);

    return μ_lims, μ_center, μ_scatterings
end


# This function is here just for testing. It is not supposed to be used in production.
function load_old_scattering_matrices(path_to_AURORA_matlab, θ_lims)
    ## Creating a MATLAB session
    s1 = MSession();
    @mput path_to_AURORA_matlab
    mat"
    addpath(path_to_AURORA_matlab,'-end')
    add_AURORA
    cd(path_to_AURORA_matlab)
    "

    theta_lims2do = reshape(Vector(θ_lims), 1, :);
    @mput theta_lims2do
    mat"
    theta_lims2do = double(theta_lims2do)
    [Pmu2mup,theta2beamW,BeamW,mu_lims] = e_scattering_result_finder(theta_lims2do,AURORA_root_directory);
    mu_scatterings = {Pmu2mup,theta2beamW,BeamW};
    c_o_mu = mu_avg(mu_lims);
    "
    θ_lims = vec(@mget theta_lims2do);
    μ_lims = vec(@mget mu_lims);
    μ_center = vec(@mget c_o_mu);

    Pmu2mup = @mget Pmu2mup; # Probability mu to mu prime
    BeamWeight_relative = @mget theta2beamW;

    # This beam weight is calculated in a continuous way
    BeamWeight = 2π .* vec(@mget BeamW);
    # Here we normalize BeamWeight_relative, as it is supposed to be a relative weighting matrix with the relative
    # contribution from withing each beam. It means that when summing up along each beam, we should get 1
    BeamWeight_relative = BeamWeight_relative ./ repeat(sum(BeamWeight_relative, dims=2), 1, size(BeamWeight_relative, 2));

    μ_scatterings = (Pmu2mup = Pmu2mup, BeamWeight_relative = BeamWeight_relative, BeamWeight = BeamWeight);

    ## Closing the MATLAB session
    close(s1)

    return μ_lims, μ_center, μ_scatterings
end


using DelimitedFiles
using PythonCall
using SpecialFunctions
using Term
function load_neutral_densities(msis_file, h_atm)
    data_msis = readdlm(msis_file, skipstart=14)
    z_msis = data_msis[:, 6]
    # import interpolate function from python
    pyinterpolate = pyimport("scipy.interpolate")
    # create the interpolator
    msis_interpolator = pyinterpolate.PchipInterpolator(z_msis, data_msis);
    # interpolate the msis data over our h_atm grid
    msis_interpolated = msis_interpolator(h_atm / 1e3)
    # the data needs to be converted from a Python array back to a Julia array
    msis_interpolated = pyconvert(Array, msis_interpolated)

    nO = msis_interpolated[:, 9] * 1e6 # from cm⁻³ to m⁻³
	nN2 = msis_interpolated[:, 10] * 1e6 # from cm⁻³ to m⁻³
	nO2 = msis_interpolated[:, 11] * 1e6 # from cm⁻³ to m⁻³

	nO[end-2:end] .= 0
	nN2[end-2:end] .= 0
	nO2[end-2:end] .= 0
    erf_factor = (erf.((1:-1:-1) / 2) .+ 1) / 2
	nO[end-5:end-3] .= erf_factor .* nO[end-5:end-3]
	nN2[end-5:end-3] .= erf_factor .* nN2[end-5:end-3]
	nO2[end-5:end-3] .= erf_factor .* nO2[end-5:end-3]

    n_neutrals = (nN2 = nN2, nO2 = nO2, nO = nO)
    return n_neutrals
end

using PythonCall
function load_electron_properties(iri_file, h_atm)
    data_iri = readdlm(iri_file, skipstart=41)
    z_iri = data_iri[:, 1]
    # import interpolate function from python
    pyinterpolate = pyimport("scipy.interpolate")
    # create the interpolator
    iri_interpolator = pyinterpolate.PchipInterpolator(z_iri, data_iri);
    # interpolate the iri data over our h_atm grid
    iri_interpolated = iri_interpolator(h_atm / 1e3)
    # the data needs to be converted from a Python array back to a Julia array
    iri_interpolated = pyconvert(Array, iri_interpolated)

    ne = iri_interpolated[:, 2] * 1e6 # from cm⁻³ to m⁻³
	Te = iri_interpolated[:, 6]
	Te[Te .== -1] .= 350

    return ne, Te
end

function load_excitation_threshold()
    file_N2_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "N2_levels.dat")
	file_O2_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "O2_levels.dat")
	file_O_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "O_levels.dat")

	N2_levels = readdlm(file_N2_levels, comments=true, comment_char='%')
	O2_levels = readdlm(file_O2_levels, comments=true, comment_char='%')
	O_levels = readdlm(file_O_levels, comments=true, comment_char='%')

	E_levels_neutrals = (N2_levels = N2_levels, O2_levels = O2_levels, O_levels = O_levels)
    return E_levels_neutrals
end
