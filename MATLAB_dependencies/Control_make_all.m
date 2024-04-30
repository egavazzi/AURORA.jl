% This script is here to control and call all the make_all_* scripts.
%
% - RunDirs: 
%   The simulation folder to analyse/plot.
%       exemple: RunDirs = {'Alfven_536s'}
%   Several folders can be analysed (if they have the same root).
%       exemple: RunDirs = {'Alfven_536s', 'Alfven_325s'}
%
% - results_dir: 
%   Root of the RunDirs folder(s).
%       exemple: results_dir = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/"


results_dir = "/Users/ost051/Documents/PhD/AURORA-20220613(1)/EtiennesVersion/AURORA.jl/data/ionprod"
cd(results_dir)
list = dir("*eV")
RunDirs = {list.name} 
%%
%make_all_Ie_top     % precipitation-spectra extraction at top
%make_all_Ie_top_MI  % precipitation-spectra extraction at top and top-1
%%
make_all_Q_lambda   % excitation and ionization-rates calculations
%make_all_I_lambda   % column-emission calculations
%%
% This has been translated to Julia
% make_all_IQ_lambda_plots    % volume/column-emission plots
%%
%cxmax = 9                 % maximum value for the colorbar(s)
%movies2make = [1 0 0 0 0] % choose animations to make
%theta_lims_2_plot = [180 160 130 110 90 70 50 20 0] % as the name says

%make_all_animations       % electron-flux animation production
