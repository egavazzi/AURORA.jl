%% Make all animations-script
% This script attempts to produce 4 animations per
% results-directory, 3 with plane-cuts through the z-t-E
% electron-flux array for each pitch-angle-stream and one with
% pitch-angle-spectra at four altitudes.

%% Intensity-limits. 
%  To make it possible or easy to have consistent intensity-limits
%  across multiple result-files there is currently no automatic
%  intensity-limits set in this script. This means that an array
%  for the colour-limits, cxmax has to be generated manually. 

%% Root result-directories
if ~exist('results_dir','var') || isempty(results_dir)
  disp('Please enter a "results_dir"')
  return
end

%% Run-directories
if ~exist('RunDirs','var') || isempty(RunDirs)
  disp('Please enter one or several "RunDirs"')
  return
end
if ~exist('movies2make')
  movies2make = [1 1 1 1 1];
end

%% Check cxmax
if ~exist('cxmax')
  disp('Please enter a value for "cxmax"')
  return
end

%% This sub-plot layout might need to be modified to adapt to layouts
% more suitable to other pitch-angle-stream configurations:
spp = [2*ones(8,1),4*ones(8,1),[1:4,(8:-1:5)]'];

%% These figure-sizes might need to be modified to adapt to other
%pitch-angle-stream configurations too.
figPwide = [692, 556, 997, 418];   % trying this one for 10-beam calculations
figPsquare = [1021, 421, 668, 553];% trying this one for 9-beam
                                   % calculations and for
                                   % pitch-angle animation
figPHigh = [1214, 326, 475, 648];
fig_sz = figPwide;


figure
whitebg([1 1 1])


%% Loop away!
for i2 = 1:numel(RunDirs)
  cd(results_dir)
  cd(RunDirs{i2})
  dDir = dir;
  
  % Variable for diagnostics
  iFaulties = 1;
  %# cxmax = 12*ones(numel(dDir),2);
  iRF = 1;
    if dDir(1).isdir
    try
      % loading the electron-transport results
      [t,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader({dDir(iRF).name});
      dE = diff(E);dE(end+1) = dE(end);
      if abs(sum(2*pi*mu_scatterings{3}) - 4*pi) < 0.01
        BeamW = 2*pi*mu_scatterings{3}; % works with old matlab results
      elseif abs(sum(mu_scatterings{3}) - 4*pi) < 0.01
        BeamW = mu_scatterings{3}; % works with newer matlab and julia results
      end
      if abs(sum(BeamW) - 4*pi) > 0.01 % because sum of mu_scattering should be = 4pi
        disp(" ")
        disp(" ")
        disp(" ")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp("WARNING : The sum of BeamW is not equal to 4pi")
        disp(" ")
        disp(" ")
        disp(" ")
        break
      end
      
      if ~exist("theta_lims_2_plot")
        theta_lims_2_plot = [180 160 130 110 90 70 50 20 0];
      end
      % make a string for the plot titles, based on the given array
      % theta_lims_2_plot
      mu_lims_2_plot = cosd(theta_lims_2_plot);
      for i1 = numel(mu_lims_2_plot)-1:-1:1
        if mu_lims_2_plot(i1) < 0
          theta_str{i1} = sprintf('%3.1f - %3.1f DOWN',...
                        180-180/pi*acos(mu_lims_2_plot(i1)),...
                        180-180/pi*acos(mu_lims_2_plot(i1+1)));
        else
          theta_str{i1} = sprintf('%3.1f - %3.1f UP',...
                        180/pi*acos(mu_lims_2_plot(i1+1)),...
                        180/pi*acos(mu_lims_2_plot(i1)));
        end
      end
      % [~, ~, BeamW_plot, ~] = e_scattering_result_finder(theta_lims_2_plot,AURORA_root_directory);
      for iMu = (numel(mu_lims_2_plot)-1):-1:1
        BeamW_plot(iMu) = 2 * pi * abs(integral(@(pa) sin(pa),acos(mu_lims_2_plot(iMu)),acos(mu_lims_2_plot(iMu+1))));
      end
      if abs(sum(BeamW_plot) - 4*pi) > 0.01 % because sum of mu_scattering should be = 4pi
        disp(" ")
        disp(" ")
        disp(" ")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp("WARNING : The sum of BeamW_plot is not equal to 4pi")
        disp(" ")
        disp(" ")
        disp(" ")
        break
      end
      % Sum the fluxes from the streams that are contained in between the
      % theta_lims_2_plot values.
      % Exemple: if we want to plot the total flux of electrons with
      % pitch-angles between 20° and 40°, it will sum the fluxes from the
      % streams between 20°-30° and 30°-40°, if these streams exist.
      IeZTE_2_plot = zeros(numel(h_atm)*(numel(theta_lims_2_plot)-1),size(IeZTE,2),size(IeZTE,3));
      try
        for i1 = 1:numel(mu_lims_2_plot)-1
          % find the index of the streams to sum
          index1 = find(abs(mu_lims - mu_lims_2_plot(i1)) < 0.001);
          index2 = find(abs(mu_lims - mu_lims_2_plot(i1+1)) < 0.001);
          if isempty(index1)
            disp(['Error : the pitch-angle to plot ', num2str(theta_lims_2_plot(i1)),...
                    '° does not match any of the stream limits used in the simulation.'])
          elseif isempty(index2)
            disp(['Error : the pitch-angle to plot ', num2str(theta_lims_2_plot(i1+1)),...
                    '° does not match any of the stream limits used in the simulation.'])
          end
          % and sum them
          for i3 = index1:(index2-1)
            hindex = (1:numel(h_atm))+(i1-1)*numel(h_atm);
            hindex_2_sum = (1:numel(h_atm))+(i3-1)*numel(h_atm);
            IeZTE_2_plot(hindex,:,:) =  IeZTE_2_plot(hindex,:,:) + IeZTE(hindex_2_sum,:,:);
          end
        end
      catch
%         disp(['Error : the pitch-angle to plot ', num2str(theta_lims_2_plot(i1)),...
%             '° does not match any of the stream limits used in the simulation'])
      end  
      
      try
        % Producing the first animation with subplots for energy
        % (x) - altitude (y) with time-variation 
        if movies2make(1)
          filename = fullfile(dDir(iRF).name,'IeztE_3DEzoft.avi');
          fprintf(['Making animation: ',filename, '\n'])
          colormap(jet)
          set(gcf,'position',fig_sz)
%             set(gcf,'WindowState','maximized');
          animate_IeztE_3DEzoft(t,h_atm,E(1:size(IeZTE,3)),...
                                IeZTE_2_plot,...
                                dE(1:size(IeZTE,3)),BeamW_plot,...
                                [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with doing: ',filename])
        Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
      
      try
        % Producing the second animation with subplots for energy
        % (x) - altitude (y) with energy-variation 
        if movies2make(2)
          filename = fullfile(dDir(iRF).name,'IeztE_3DtzofE.avi');
          fprintf(['Making animation: ',filename, '\n'])
          colormap(jet)
          set(gcf,'position',fig_sz)
          animate_IeztE_3DtzofE(t,h_atm,E(1:size(IeZTE,3)),...
                                IeZTE_2_plot,...
                                dE(1:size(IeZTE,3)),BeamW_plot,...
                                [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with doing: ',filename])
        Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
      
      try
        % Producing the third animation with subplots for time
        % (x) - Energy (y) with altitude-variation 
        if movies2make(3)
          filename = fullfile(dDir(iRF).name,'IeztE_3DtEofz.avi');
          fprintf(['Making animation: ',filename, '\n'])
          set(gcf,'position',fig_sz)
          colormap(jet)
          animate_IeztE_3DtEofz(t,h_atm,E(1:size(IeZTE,3)),...
                                IeZTE_2_plot,...
                                dE(1:size(IeZTE,3)),BeamW_plot,...
                                [-5 0]+max(cxmax(min(end,iRF),:)),spp, theta_str,filename);
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with doing: ',filename])
        Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
      
      try
        % Producing the fourth animation with pitch-angle
        % distribution at four heights
        if movies2make(4) 
          filename = fullfile(dDir(iRF).name,'IeztE_pitchangledist.avi');
          fprintf(['Making animation: ',filename, '\n'])
          [dh,i115] = min(abs(h_atm/1e3-225));
          [dh,i175] = min(abs(h_atm/1e3-325));
          [dh,i300] = min(abs(h_atm/1e3-450));
          [dh,i600] = min(abs(h_atm/1e3-600));
          iZ = [i115,i175,i300,i600];
          set(gcf,'position',figPsquare)
          animate_IeztE_pitchangleflux(t,h_atm/1e3,...
                                        E,dE,...
                                        IeZTE_2_plot,...
                                        BeamW_plot,mu_lims_2_plot,...
                                        iZ,...
                                        filename,...
                                        [-3.5 0]+max(cxmax(min(end,iRF),:)));
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with doing: ',filename])
        Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
      %%
      try
        % Producing the animation with fluxes as frunction of pitch-angle
        % and height at highest energy
        if movies2make(5)
          filename = fullfile(dDir(iRF).name,'IeztE_mu_z_at_E.avi');
          set(gcf,'position',figPHigh)
          fprintf(['Making animation: ',filename, '\n'])
          animate_IeztE_3DzmuatEoft(t,h_atm,E,...
                                    IeZTE_2_plot,...
                                    dE,BeamW_plot,...
                                    size(IeZTE,3),...
                                    [6.5 12.5]-4,...
                                    ... {'0','10','30','60','80','90','100','120','150','170','180'},...
                                    theta_str, ...
                                    filename);
        end
      catch
        ERR_MSG{iFaulties} = lasterr;
        disp(['something wrong with doing: ',filename])
        Faulties{iFaulties} = filename;
        iFaulties = iFaulties + 1;
      end
    catch
      ERR_MSG{iFaulties} = lasterr;
      disp(['something wrong with processing directory: ',dDir(iRF).name])
      Faulties{iFaulties} = dDir(iRF).name;
      iFaulties = iFaulties + 1;
    end
  end
end


%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
