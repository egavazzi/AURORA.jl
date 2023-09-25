function Xs = N2_e_2nd_dist(Es,E0in,Eionization,s_or_c,AURORA_root_directory)
% N2_e_2nd_dist - 2nd:ary electron energy spectra ionisation of N2.
%
% Calling:
%  Xs = N2_e_2nd_dist(Es,E0in)
% Input:
%  Es   - energy of secondar electrons (eV), double array [1 x nEs]
%  E0in - energy of primary electrons (eV), double scalar
%  Eionization - ionization threshold (eV), double scalar so that
%         we can properly calculate the energy-distribution of
%         degrading electrons due to ionization to excited N2+ and
%         dissociative ionizations
%  s_or_c - spectra of secondary electrons or cascading electrons
% Output:
%  Xs - spectra of secondary electrons, double array [1 x nEs]
% 
% NOTE: Unnormalized! So that the user have to calculate the total
% ionization by normalizing Xs by its sum and scale with the
% ionizing cross-section, the N2 density and the electron flux.
% 
% Equation from Ittikawa 1986 J. Phys. Chem. Ref. Data

%  Copyright � Bjorn Gustavsson 20190408, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


persistent Q E4Q Eionizations
if nargin < 5
  disp('Error with the N2_e_2nd_dist function, lacks the AURORA root directory argument')
end

% trying to find a cascading spectra file with matching energy grid Es 
if isempty(Q) || ~all(E4Q(1:min(numel(Es),numel(E4Q)))==Es(1:min(numel(Es),numel(E4Q)))) || size(Es,2) > size(E4Q,2)
  try
    e_2nd_s_files = dir(fullfile(AURORA_root_directory,'E_cascadings','N2'));
    for i1 = 1:numel(e_2nd_s_files),
      if ~ e_2nd_s_files(i1).isdir
        try
          load(fullfile(AURORA_root_directory,...
                        'E_cascadings','N2',...
                        e_2nd_s_files(i1).name),...
               'E4Q')
%           if isequal(E4Q,Es)
            if all(E4Q(1:min(numel(Es),numel(E4Q)))==Es)
            % then we have found a match, so load
            fprintf('Loading cascading-matrices from file: %s\n',e_2nd_s_files(i1).name)
            load(fullfile(AURORA_root_directory,...
                         'E_cascadings','N2',...
                         e_2nd_s_files(i1).name),...
                'E4Q',...
                'Q',...
                'Eionizations')
            foundem = 1;
            % And break the loop already
            break
          end
        catch
        end
      end
    end
  catch
  fprintf('Could not find file with matching energy grid\n')
  fprintf('Starting to calculate the requested cascading-matrices')
  end
end

%% Precalculation of degrading spectra
E_hat = 11.4;
if isempty(Q) || ~all(E4Q(1:min(numel(Es),numel(E4Q)))==Es(1:min(numel(Es),numel(E4Q)))) || size(Es,2) > size(E4Q,2)
  % ~all(E4Q(1:numel(Es))==Es)
  Eionizations = [15.581
                  16.73
                  18.75
                  24
                  42];
  Q = zeros(numel(Es),numel(Es),numel(Eionizations));
  E4Q = Es;
  dE = diff(Es);dE = dE([1:end,end]);
  disp('Pre-calculating all energy-degradations for e - N2-ionizing collisions.')
  disp('This will take a bit of time.')
  disp(['Starting at: ',datestr(now,'HH:MM:SS')])
  for i3 = numel(Eionizations):-1:1,
    iLim = find(Es>Eionizations(i3),1,'first');
    for i2 = iLim:numel(Es),
      iHalf = find(Es<(Es(i2)-Eionizations(i3))/2,1,'last');
      for i1 = iHalf:(i2-1),
        Eprime = [Es(i2),Es(i2)+dE(i2)];
        Edeg = [Es(i1),Es(i1)+dE(i1)];
        fun =  @(E1,E2) 1./((E2-Eionizations(i3)-E1).^2+E_hat^2);
        % Lower boundary of primary electrons, This correctly
        % accounts for the limit when Edeg(i1+1) is larger than
        % Eprimary(i2+1)-Eionizations(i3), cutting a corner off the
        % [Ei1 - Ei1+dEi1] x [Ei2 + dEi2] square we're integrating over.
        Emin = @(E) max(Eprime(1),E+Eionizations(i3));
        % Correct upper Edeg integration limit, this cuts the
        % region into a triangle when needed
        Edegmax = min(Es(i1)+dE(i1),Es(i2)+dE(i2)-Eionizations(i3));
        % This if-condition just makes sure that we have the
        % physically meaningful limits, i.e. we integrate in
        % increasing Edeg - which is not otherwise guaranteed.
        if Edegmax > Es(i1)
          % Edges(i1) = {Emin,Edegmax};
          Q(i2,i1,i3) = integral2(fun,Es(i1),Edegmax,Emin,Eprime(end));
        end
      end
    end
    disp(['Done with level ',num2str(i3),' at: ',datestr(now,'HH:MM:SS')])
  end

  % Save the results for future use
  save_filename = sprintf('CascadingSpecN2ionization_%s.mat',...
  datestr(now,'yyyymmdd-HHMMSS'));
  save(fullfile(AURORA_root_directory,...
  'E_cascadings','N2',...
  save_filename),...
  'Q',...
  'E4Q',...
  'Eionizations')
end

% $$$         funH = @(E1,E2) heaviside(E2-Eionizations(i3)-E1)./((E2-Eionizations(i3)-E1).^2+E_hat^2);
% Q2(i2,i1,i3) = integral2(funH,Es(i1),Edegmax,Emin,Eprime(end));

if strcmp(s_or_c,'s')
  %% Calculating the spectra of secondary electrons
  Xs = 1./(E_hat^2+Es.^2).*(Es < (E0in - Eionization)/2 );
  
else
  %% or extracting the degrading primary electron
  [DE,iLevel] = min(abs(Eionization-Eionizations));
  [DE2,iPrimary] = min(abs(E4Q-E0in));
  Xs = Q(iPrimary,1:numel(Es),iLevel);
  
end
