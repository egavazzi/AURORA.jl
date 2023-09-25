function make_all_Ie_at_z(altitude_to_extract, results_dir, RunDirs)
%% Electron-flux-at-top-of-ionosphere-extraction
% This script extracts the electron-fluxes at a given altitude for all
% pitch-angle-streams.
%
% Calling:
%   make_all_Ie_at_z(altitude_to_extract);
% Input:
%  altitude_to_extract - (km)

%% Root result-directories
if ~exist('results_dir','var') || isempty(results_dir)
	disp('Please enter a "results_dir"')
	return
end

%% Run-directories:
if ~exist('RunDirs','var') || isempty(RunDirs)
	disp('Please enter one or several "RunDirs"')
	return
end

for i2 = 1:numel(RunDirs)
	cd(results_dir)
	cd(RunDirs{i2})
	dDir = dir;
	if dDir(1).isdir
		CD = pwd;
		try
		[t,h_atm,E,mu_lims,Ie_ZTE,mu_scatterings] = Ie_ztE_loader({'.'});
		dE = diff(E);
		dE = dE([1:end,end]);
		BeamW = mu_scatterings{3};
		nZ = numel(h_atm);

		% Find the iz of the altitude_to_extract
		[~, iZ] = min(abs(h_atm - altitude_to_extract * 1e3));
		% and extract
		Ie_raw = Ie_ZTE(iZ:nZ:end,:,:);
		altitude_extracted = h_atm(iZ) / 1e3
		% then save
		name_of_the_file = ['Ie_', num2str(round(h_atm(iZ) / 1e3)), 'km.mat']
		save(name_of_the_file, 'Ie_raw', 'altitude_to_extract', 'altitude_extracted', 'E','dE','BeamW','t')
		clear Ie_raw
		fprintf(':::Processed OK: %s\n',CD)
		catch
			fprintf('Everythings not right in directory: %s\n',CD)
		end
	end
	cd(results_dir)
	cd(RunDirs{i2}) 
end

end