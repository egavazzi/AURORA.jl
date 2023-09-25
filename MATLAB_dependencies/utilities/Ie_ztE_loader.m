function [t,h_atm,E,mu_lims,IeZTE,mu_scatterings] = Ie_ztE_loader(data_paths,ops)
% Ie_ztE_loader reads flickering-aurora data from the
% time-dependent electron transport runs for flickering
% aurora. This run was done in sequences of short time-intervalls,
% typically 0.1 s, with initial fluxes from the end of the
% preceding time-period. To ease further analysis we need to strip
% off the fluxes at the first time-step from data from the second
% file onwards. This function allows loading of data from two
% separate runs, stored in separate directories for addition or
% subtraction of the fluxes - this provided that the altitude,
% pithc-angle, time and energy-grids all match up
% Input:
%  data_path - path, absolute or relative path to directory with
%              output from the time-dependent electron transport
%              runs for flickering aurora, cell array [1 x 1] or
%              [1 x 2]
%  ops      - options struct with fields: operator, display, this
%             is an optional input, if not given and the data_paths
%             operator, [ {'+'} | '-' ]
%             display [{'quiet'},'verbose']
% Output:
%  t        - time [1 x n_t] double aray with time in suitable
%             units for dispolaying time-variation of densities
%  h_atm    - altitude (m), double array [n_z x 1]
%  E        - energy grid (eV), double array [1 x n_E]
%  mu_lims  - cosine-of-pitch-angle-limits for the number of
%             streams, double array [1 x (n_streams+1)]
%  IeZTE    - electron number flux (#/m^2/s), double array 
%             [(n_z*n_streams) x n_t x n_e] with number of
%             electrons per energy-bin, sorted along the first
%             dimension as stream-above-stream from most paralell
%             to B downward to most parallel to B upward.
%  mu_scatterings - cell-array with the stream-to-stream arrays,
%             the pitch-angle-size/(2*pi)

%   Copyright � 2018-2019 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later



if nargin < 1 || isempty(data_paths)
  data_paths{1} = pwd;
end
if nargin < 2 || isempty(ops)
  ops.operator = '+';
  ops.display  = 'quiet';
end


if numel(data_paths) == 1 % Simplest case with data from one directory
  
  %  E               1x94                     752  double
  %  Ie_ztE       4833x31x94            112666896  double
  %  h_atm         537x1                     4296  double
  %  mu_lims         1x10                      80  double
  %  t_run           1x31                     248  double
  if ~strcmp(ops.display,'quiet')
    res_files = dir(fullfile(data_paths{1},'IeFl*.mat'));
    disp(['loading: ',fullfile(data_paths{1},res_files(1).name)])
  else
    res_files = dir(fullfile(data_paths{1},'IeFl*.mat'));
  end
  load(fullfile(data_paths{1},res_files(1).name))
  [mu_scatterings, E, t_run, mu_lims] = convert_from_julia(mu_scatterings, E, t_run, mu_lims);
  IeZTE = Ie_ztE;
  t = t_run;
  for i1 = 2:length(res_files),
    if ~strcmp(ops.display,'quiet')
      disp(['loading: ',fullfile(data_paths{1},res_files(i1).name)])
    end
    load(fullfile(data_paths{1},res_files(i1).name),...
        'E','Ie_ztE','h_atm','mu_lims','t_run')
    [mu_scatterings, E, t_run, mu_lims] = convert_from_julia(mu_scatterings, E, t_run, mu_lims);
    IeZTE = cat(2,IeZTE,Ie_ztE(:,2:end,:));
    t = [t,t_run(2:end)];
  end
  
elseif numel(data_paths) == 2 && ops.operator == '+'% Case with data from two directories and ops:'+'
  res_files1 = dir(fullfile(data_paths{1},'IeFl*.mat'));
  res_files2 = dir(fullfile(data_paths{2},'IeFl*.mat'));
  
  disp(['loading: ',fullfile(data_paths{1},res_files1(1).name)])
  load(fullfile(data_paths{1},res_files1(1).name))
  Ie1 = Ie_ztE;
  disp(['loading: ',fullfile(data_paths{2},res_files2(1).name)])
  load(fullfile(data_paths{2},res_files2(1).name))
  switch ops.operator
   case '+' % Case with data from two directories and ops:'+'
    IeZTE = Ie1 + Ie_ztE;
   case '-' % Case with data from two directories and ops:'-'
    IeZTE = Ie1 - Ie_ztE;
   otherwise
  end
  t = t_run;
  for i1 = 2:length(res_files1),
    disp(['loading: ',fullfile(data_paths{1},res_files1(i1).name)])
    load(fullfile(data_paths{1},res_files1(i1).name))
    Ie1 = Ie_ztE;
    disp(['loading: ',fullfile(data_paths{2},res_files2(i1).name)])
    load(fullfile(data_paths{2},res_files2(i1).name))
    IeZTE = cat(2,IeZTE,Ie1(:,2:end,:)+Ie_ztE(:,2:end,:));
    t = [t,t_run(2:end)];
  end
  
elseif numel(data_paths) == 2 && ops.operator == '-'% Case with data from two directories and ops:'-'
  res_files1 = dir(fullfile(data_paths{1},'IeFl*.mat'));
  res_files2 = dir(fullfile(data_paths{2},'IeFl*.mat'));
  
  disp(['loading: ',fullfile(data_paths{1},res_files1(1).name)])
  load(fullfile(data_paths{1},res_files1(1).name))
  Ie1 = Ie_ztE;
  disp(['loading: ',fullfile(data_paths{2},res_files2(1).name)])
  load(fullfile(data_paths{2},res_files2(1).name))
  switch ops.operator
   case '+' % Case with data from two directories and ops:'+'
    IeZTE = Ie1 + Ie_ztE;
   case '-' % Case with data from two directories and ops:'-'
    IeZTE = Ie1 - Ie_ztE;
   otherwise
  end
  t = t_run;
  for i1 = 2:length(res_files1),
    disp(['loading: ',fullfile(data_paths{1},res_files1(i1).name)])
    load(fullfile(data_paths{1},res_files1(i1).name))
    Ie1 = Ie_ztE;
    disp(['loading: ',fullfile(data_paths{2},res_files2(i1).name)])
    load(fullfile(data_paths{2},res_files2(i1).name))
    IeZTE = cat(2,IeZTE,Ie1(:,2:end,:)-Ie_ztE(:,2:end,:));
    t = [t,t_run(2:end)];
  end
  
else
  disp('I''m afraid I cannot do that.')
end
end




function [mu_scatterings, E, t_run, mu_lims] = convert_from_julia(mu_scatterings, E, t_run, mu_lims)
  % if we load data produced by the Julia version of Aurora, the type/class
  % of mu_scatterings is going to be different. To remediate to this problem,
  % we convert mu_scatterings in the following block
  if isa(mu_scatterings, 'struct')
    mu_scatterings_temp{1} = mu_scatterings.Pmu2mup;
    mu_scatterings_temp{2} = mu_scatterings.BeamWeight_relative;
    mu_scatterings_temp{3} = mu_scatterings.BeamWeight_discrete;

    mu_scatterings = mu_scatterings_temp;
  end
  % similar thing, dimensions of E are inversed between Julia and Matlab
  if size(E, 1) > size(E, 2)
    E = E';
  end
  if size(t_run, 1) > size(t_run, 2)
    t_run = t_run';
  end
  if size(mu_scatterings{3}, 1) > size(mu_scatterings{3}, 2)
    mu_scatterings{3} = mu_scatterings{3}';
  end
  if size(mu_lims, 1) > size(mu_lims, 2)
    mu_lims = mu_lims';
  end
end