%% Electron-flux-at-top-of-ionosphere-extraction
% This script extracts the electron-fluxes at the top of the
% ionosphere for all pitch-angle-streams.


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

%% Loop away!
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
      nZ = numel(h_atm);
      Ie_top_raw = Ie_ZTE(nZ:nZ:end,:,:);
      Ie_top = Ie_top_raw;
      BeamW = mu_scatterings{3};
      for i1 = size(Ie_top,1):-1:1,
        % to electrons per steradian from electrons  [m^-2s^-1]
        Ie_top(i1,:,:) = Ie_top(i1,:,:)/BeamW(i1)/2/pi;
      end
      for i3 = size(Ie_top,3):-1:1,
        % to electrons per steradian per eV per m^2 per s 
        Ie_top(:,:,i3) = Ie_top(:,:,i3)/dE(i3);
      end
      save('Ie_top.mat','Ie_top','Ie_top_raw','E','dE','BeamW','t')
      clear Ie_top Ie_top_raw
      fprintf(':::Processed OK: %s\n',CD)
    catch
      fprintf('Everythings not right in directory: %s\n',CD)
    end
  end
end

cd(results_dir)
cd(RunDirs{i2})


%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
