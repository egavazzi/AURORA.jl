results_dir = "/mnt/data/oliver/AURORA.jl/data/ionprod/9942.457296427418-10174.634027476064eV"
addpath("/home/bgu001/matlab/AURORA/Plotting")
cd(results_dir)
dDir = dir;
if dDir(1).isdir
    Ie_matfiles = dir('IeFlickering-*.mat');
    Qzt_file = dir('Qzt_all_L.mat');
    atm_file = dir('neutral_atm.mat');
    if ~isempty(atm_file)
      load(atm_file.name,'nN2','nO2','nO','h_atm')
    end
    [t,h_atm,E,mu_lims,Ie_ZTE] = Ie_ztE_loader({'.'});
end
dE = diff(E);
dE(end +1) = dE(end);
BeamSA = 2*pi * diff(mu_lims)
spp = [2 * ones(18, 1), 9*ones(18, 1), [1:9,18:-1:10]']
cx_lims = []
theta_lims_str = string(acos(mu_lims)/pi*180)
theta_strs = theta_lims_str(1:end-1) +'-'+ theta_lims_str(2:end)
theta_strs = cellstr(theta_strs)
%
plot_IezE_2DEz(h_atm,E,Ie_ZTE,dE,BeamSA,cx_lims,spp, theta_strs)