"""
This is the setup function. Loading the atmosphere, the energy grid, collision cross-sections, ... \\
It calls a lot of functions from the original MATLAB code.

# Inputs
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation
- `θ_lims`: range of angles for the limits of our electron beams, i.e 180:-10:0
"""
function setup(top_altitude, θ_lims)
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
    nO = @mget nO;
    nN2 = @mget nN2;
    nO2 = @mget nO2;

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

    ## Energy grid
    mat"
    dEfcn = @(X,DEinitial,DEfinal,C,X0) DEinitial+(1+tanh(C*(X-X0)))/2*DEfinal;
    E = cumsum(dEfcn(0:1000,0.15,11.5,0.05,80))+1.9;
    dE = diff(E);dE = [dE,dE(end)];
    "
    E = vec(@mget E);
    dE = vec(@mget dE);

    ## Collision cross-sections
    mat"
    [XsO,xs_fcnO] =   get_all_xs('O',E+dE/2);
    [XsO2,xs_fcnO2] = get_all_xs('O2',E+dE/2);
    [XsN2,xs_fcnN2] = get_all_xs('N2',E+dE/2);
    "
    σ_N2 = @mget XsN2;
    σ_O2 = @mget XsO2;
    σ_O = @mget XsO;

    ## Pre-calculations of cascadning electron-spectra for ionizations
    mat"
    S2ndO = O_e_2nd_dist(E,E(end),O_levels(end,1),'c',AURORA_root_directory);
    S2ndO2 = O2_e_2nd_dist(E,E(end),O2_levels(end,1),'c',AURORA_root_directory);
    S2ndN2 = N2_e_2nd_dist(E,E(end),N2_levels(end,1),'c',AURORA_root_directory);
    "
    secondary_e_O = vec(@mget S2ndO);
    secondary_e_O2 = vec(@mget S2ndO2);
    secondary_e_N2 = vec(@mget S2ndN2);

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

    Pmu2mup = @mget Pmu2mup; # Great name
    θ_to_BeamWeight = @mget theta2beamW;
    BeamWeight = vec(@mget BeamW);
    μ_scatterings = (Pmu2mup, θ_to_BeamWeight, BeamWeight);

    ## Closing the MATLAB session
    close(s1)

    return h_atm, nN2, nO2, nO, ne, Te, E, dE, 
        N2_levels, O2_levels, O_levels, 
        σ_N2, σ_O2, σ_O,
        secondary_e_N2, secondary_e_O2, secondary_e_O,
        θ_lims, μ_lims, μ_center, μ_scatterings
end