using MATLAB

"""
This is the setup function, loading the atmosphere, the energy grid, the collision cross-sections, ... \\
It calls a lot of functions from the original MATLAB code. 

# Inputs
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation
- `θ_lims`: range of angles for the limits of our electron beams, i.e 180:-10:0
- `E_max`: upper limit for the energy grid (in eV)
"""
function setup(top_altitude, θ_lims, E_max)
    ## Creating a MATLAB session
    s1 = MSession();
    mat"
    addpath('/mnt/data/etienne/AURORA','-end')
    add_AURORA
    cd /mnt/data/etienne/AURORA/
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
    [XsO,xs_fcnO] =   get_all_xs('O',E+dE/2);
    [XsO2,xs_fcnO2] = get_all_xs('O2',E+dE/2);
    [XsN2,xs_fcnN2] = get_all_xs('N2',E+dE/2);
    "
    σ_N2 = @mget XsN2;
    σ_O2 = @mget XsO2;
    σ_O = @mget XsO;
    σ_neutrals = (σ_N2 = σ_N2, σ_O2 = σ_O2, σ_O = σ_O);

    ## Pre-calculations of cascading electron-spectra for ionizations
    mat"
    S2ndO = O_e_2nd_dist(E,E(end),O_levels(end,1),'c',AURORA_root_directory);
    S2ndO2 = O2_e_2nd_dist(E,E(end),O2_levels(end,1),'c',AURORA_root_directory);
    S2ndN2 = N2_e_2nd_dist(E,E(end),N2_levels(end,1),'c',AURORA_root_directory);
    "    
    secondary_e_N2 = vec(@mget S2ndN2);
    secondary_e_O2 = vec(@mget S2ndO2);
    secondary_e_O = vec(@mget S2ndO);
    secondary_e = (secondary_e_N2 = secondary_e_N2, secondary_e_O2 = secondary_e_O2, secondary_e_O = secondary_e_O);

    ## X-streams beam-to-beam calculations
    theta_lims2do = reshape(Vector(θ_lims), 1, :);
    @mput theta_lims2do
    mat"
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
    BeamWeight_continuous = 2π .* vec(@mget BeamW);
    # Whereas this one is calculated in a discrete way. This is to ensure the conservation of the number of e-
    # when calculating the scatterings, as these calculations are discretized in pitch angle (around 721
    # different angles) which leads to a slightly different normalization factor
    BeamWeight_discrete = sum(BeamWeight_relative, dims=2); 
    BeamWeight_discrete[μ_center .< 0] = 2π .* BeamWeight_discrete[μ_center .< 0] ./ sum(BeamWeight_discrete[μ_center .< 0]);
    BeamWeight_discrete[μ_center .> 0] = 2π .* BeamWeight_discrete[μ_center .> 0] ./ sum(BeamWeight_discrete[μ_center .> 0]);

    # Here we normalize BeamWeight_relative, as it is supposed to be a relative weighting matrix with the relative
    # contribution from withing each beam. It means that when summing up along each beam, we should get 1
    BeamWeight_relative = BeamWeight_relative ./ repeat(sum(BeamWeight_relative, dims=2), 1, size(BeamWeight_relative, 2));



    μ_scatterings = (Pmu2mup = Pmu2mup, BeamWeight_relative = BeamWeight_relative, BeamWeight_continuous = BeamWeight_continuous, BeamWeight_discrete = BeamWeight_discrete);

    ## Closing the MATLAB session
    close(s1)

    return h_atm, n_neutrals, ne, Te, E, dE, 
        E_levels_neutrals, σ_neutrals, secondary_e,
        θ_lims, μ_lims, μ_center, μ_scatterings
end