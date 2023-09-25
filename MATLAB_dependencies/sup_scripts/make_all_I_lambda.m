%% Column-emission-calculations
% This script integrates the volume-emission(excitation)-rates
% taking photon time-of-flight into account. 

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
    try
      % Loading all volume-emission(excitation)-rates
      load Qzt_all_L.mat Q4278 Q6730 Q7774 Q8446 QO1D QO1S h_atm t Q7774_O Q8446_O Q7774_O2 Q8446_O2
      if numel(t) ~= size(Q4278,2)
        t = [t,t(end)+t(2:end)];
      end
      % Integrating along photon-trajectories/in altitude
      I_4278 = q2colem(t,h_atm,Q4278);
      I_6730 = q2colem(t,h_atm,Q6730);
      I_7774 = q2colem(t,h_atm,Q7774);
      I_8446 = q2colem(t,h_atm,Q8446);
      I_O1D  = q2colem(t,h_atm,QO1D);
      I_O1S  = q2colem(t,h_atm,QO1S);
      
      I_7774_O = q2colem(t,h_atm,Q7774_O);
      I_8446_O = q2colem(t,h_atm,Q8446_O);
      I_7774_O2 = q2colem(t,h_atm,Q7774_O2);
      I_8446_O2 = q2colem(t,h_atm,Q8446_O2);
      % and save
      save I_lambda_of_t.mat t I_4278 I_6730 I_7774 I_8446 I_O1D I_O1S I_7774_O I_7774_O2 I_8446_O I_8446_O2
      fprintf('### Thing worked out OK in directory: %s\n',pwd)
    catch
      fprintf('Everythings not right in directory: %s\n',pwd)
    end
  end
end


%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
