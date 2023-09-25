% This script is here to control and call all the make_all_* scripts.
%
% - RunDirs: 
% The simulation folder to analyse/plot.
%       exemple: RunDirs = {'Alfven_536s'}
% Several folders can be analysed (if they have the same root).
%       exemple: RunDirs = {'Alfven_536s', 'Alfven_325s'}
%
% - results_dir: 
% Root of the RunDirs folder(s).
%       exemple: results_dir = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/"



results_dir = ""
RunDirs = {''} 
%%
make_all_Ie_top     % precipitation-spectra extraction at top
make_all_Ie_top_MI  % precipitation-spectra extraction at top and top-1
%%
make_all_Q_lambda   % excitation and ionization-rates calculations
make_all_I_lambda   % column-emission calculations
%%
make_all_IQ_lambda_plots    % volume/column-emission plots
%%
cxmax = 10                % maximum value for the colorbar(s)
movies2make = [0 0 0 1 0] % choose animations to make
theta_lims_2_plot = [180 165 135 105 90 75 45 15 0] % as the name says

make_all_animations       % electron-flux animation production
