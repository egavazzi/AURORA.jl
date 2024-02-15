%% Volume and Column intensity-plots
% This script produces plots of altitude-time-variation of volume
% emission-rates and time-variation of normalized
% column-emission(excitation) intensity plots.

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

%% Optional setting of legend location:
if ~exist('leg_location')
%   leg_location = 'northwest';
%   leg_location = 'northeast';
  leg_location = 'eastoutside';
end
% For possible options see documentation for legend.

%% Hard-coded output-directory
FigDir = sprintf('Figures-%s',datestr(now,'yyyymmdd'))
mkdir(fullfile(results_dir,FigDir))

%% Loop away!
for i2 = 1:numel(RunDirs)
  cd(results_dir)
  cd(RunDirs{i2})
  dDir = dir;
	if dDir(1).isdir	
		try
			load Qzt_all_L.mat Q4278 Q6730 Q7774 Q8446 QO1D QO1S h_atm t
		if numel(t) ~= size(Q4278,2)
			t = [t,t(end)+t(2:end)];
		end
		clf
		colormap('default')
		subplot(3,2,1)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,Q4278),
		title('4278')
		shading flat
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('photons/m^3/s')
		subplot(3,2,2)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,Q6730),
		title('6730')
		shading flat
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('photons/m^3/s')
		subplot(3,2,3)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,Q7774),
		title('7774')
		ylabel('height (km)')
		shading flat
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('photons/m^3/s')
		subplot(3,2,4)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,Q8446),
		title('8446')
		shading flat
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('#exc/m^3/s')
		subplot(3,2,5)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,QO1D),
		title('O1D')
		shading flat
		xlabel('time (s)')
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('#exc/m^3/s')
		subplot(3,2,6)
		pcolor(t(1:size(Q4278,2)),h_atm/1e3,QO1S),
		title('O1S')
		shading flat
		xlabel('time (s)')
		ax = axis;
		%         axis([ax(1:3) 300])
		colorbar_labeled('photons/m^3/s')
		orient tall
		        % print('-depsc2','-painters',fullfile(results_dir,FigDir,strcat(RunDirs{i2},'-Qtz-01.eps')))
		        print('-dpng','-painters',fullfile(results_dir,FigDir,strcat(RunDirs{i2},'-Qtz-01.png')))
		load I_lambda_of_t.mat t I_4278 I_6730 I_7774 I_7774_O I_7774_O2 I_8446 I_8446_O I_8446_O2 I_O1D I_O1S
		figure
		%         clf
		%         ph = plot(t(1:numel(I_4278)),...
		%                    100*[I_4278/max(I_4278);
		%                        I_6730/max(I_6730);
		%                        I_7774/max(I_7774);
		%                        I_8446/max(I_8446);
		%                        I_O1D/max(I_O1D);
		%                        I_O1S/max(I_O1S);
		%                        I_7774_O/max(I_7774);
		%                        I_7774_O2/max(I_7774);
		%                        I_8446_O/max(I_8446);
		%                        I_8446_O2/max(I_8446)]);
		ph = semilogy(t(1:numel(I_4278)),...
										[I_4278;
										I_6730;
										I_7774;
										I_8446;
										I_O1D;
										I_O1S;
										]);
		xlabel('time (s)','fontsize',15)
		%         ylabel('(%)','fontsize',15)
		ylabel('#exc/m^2∕s','fontsize',15)
		%         title('Intensity Modulation','fontsize',15)
		title('Intensity','fontsize',15)
		ax = axis;
    axis([ax(1:2) 1e4 ax(4)*1e1])
% 		axis([ax(1:3) 105])
		set(ph,'linewidth',2)
		set(ph(1),'color','b')
		set(ph(2),'color','r')
		set(ph(3),'color',[0.5 0 0])
		set(ph(4),'color','k')
		set(ph(5),'color',[1 0.2 0],'linestyle','--','linewidth',1)
		set(ph(6),'color','g','linestyle','--','linewidth',1)
    grid on
		%         set(ph(7),'color',[0.5 0 0],'linestyle','-.','linewidth',1)
		%         set(ph(8),'color',[0.5 0 0],'linestyle','--','linewidth',1)
		% 
		%         set(ph(9),'color','k','linestyle','-.','linewidth',1)
		%         set(ph(10),'color','k','linestyle','--','linewidth',1)
		%         
		% legend(ph,...
		% 		'I_{4278}',...
		% 		'I_{6730}',...
		% 		'I_{7774}',...
		% 		'I_{8446}',...
		% 		'I_{O(^1D)}',...
		% 		'I_{O(^1S)}',...
		% 		'I_{7774}(O)',...
		% 		'I_{7774}(O_2)',...
		% 		'I_{8446}(O)',...
		% 		'I_{8446}(O_2)',...
		% 		'location',leg_location)
		legend(ph,...
				'I_{4278}',...
				'I_{6730}',...
				'I_{7774}',...
				'I_{8446}',...
				'I_{O(^1D)}',...
				'I_{O(^1S)}',...
				'location',leg_location)
		orient portrait
		% Figures saved as both .png and .eps-files
		        % print('-depsc2','-painters',fullfile(results_dir,FigDir,strcat(RunDirs{i2},'-It-01.eps')))
		        print('-dpng','-painters',fullfile(results_dir,FigDir,strcat(RunDirs{i2},'-It-01.png')))
		catch
			fprintf('Everythings not right in directory: %s\n',pwd)
		end
	end
end


%   Copyright � 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
