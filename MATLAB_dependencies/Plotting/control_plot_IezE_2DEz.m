results_dir = "/Users/ost051/Documents/PhD/AURORA-20220613(1)/EtiennesVersion/AURORA.jl/data/ionprod/10.0-10.23352047097258eV"
addapth("/Users/ost051/Documents/PhD/AURORA-20220613(1)/AURORA/Plotting")
cd(results_dir)
dDir = dir;
if dDir(1).isdir
    Ie_matfiles = dir('IeFlickering-*.mat');
    Qzt_file = dir('Qzt_all_L.mat');
    atm_file = dir('neutral_atm.mat');
    if ~isempty(atm_file)
      load(atm_file.name,'nN2','nO2','nO','h_atm')
    end
end

BeamSA = 2*pi * diff(mu_lims)
spp = [2 * ones(18, 1), 9*ones(18, 1), [1:9,18:-1:10]']
plot_IezE_2DEz(h_atm,E,Ie_ZTE,dE,BeamSA,cx_lims,spp, theta_strs)