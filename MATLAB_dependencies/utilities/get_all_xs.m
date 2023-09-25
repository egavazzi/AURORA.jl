function [Xs,xs_fcn] = get_all_xs(name,E)
% GET_ALL_XS - collect all electron-neutral impact cross sections
% All cross sections for species 'NAME' at energy E (in eV)
% NAME should be the name of a species 'O', 'N2' or 'O2',
% for which to calculate the electron cross section. E - electron
% energy in electron volts. The cross sections are typicaly
% obtained from the compilations done by Itikawa et al.
% 
% Calling:
%  [Xs,xs_fcn] = get_all_xs(name,E);
% Input:
%  name - name of species,{'O','O2','N2'}
%  E    - Electron energy vector [1 x n] (eV)
% Output:
%  Xs - electron neutral collision cross sections (m^2). The
%       collision cross sections are sorted: [elastic;
%       inelastic;ionization], the excitation energies and the name
%       of the excited states are to be found in files
%       [O/N2/O2]_level.[name/dat].
%  xs_fcn - function handle to function giving the cross sections.
%

%  Copyright © Bjorn Gustavsson 20030403, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

fp = fopen([name,'_levels.name'],'r');

if fp > 0 
  
  % state_names = fscanf(fp,'%c');
  state_names = textscan(fp,'%s');
  fclose(fp);
  
  fcn_i = 1;
  for i1 = 1:numel(state_names{1})
    if state_names{1}{i1}(1) ~= '%'
      fcn{fcn_i} = str2func(['e_',name,state_names{1}{i1}]);
      fcn_i = fcn_i + 1;
    end
  end
%   [T,state_names] = strtok(state_names);
%   fcn{fcn_i} = str2func(['e_',name,T]);
%   
%   while ~isempty(state_names)
%     
%     fcn_i = fcn_i + 1;
%     [T,state_names] = strtok(state_names);
%     if ~isempty(T)
%       fcn{fcn_i} = str2func(['e_',name,T]);
%     end
%     
%   end
  
end

for fcn_i = length(fcn):-1:1,
  
  Xs(fcn_i,:) = feval(fcn{fcn_i},E);
  
end

if nargout==2
  
  xs_fcn = fcn;
  
end
