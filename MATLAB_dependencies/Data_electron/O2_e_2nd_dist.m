function Xs = O2_e_2nd_dist(Es,E0in,Eionization,s_or_c,AURORA_root_directory)
% O2_e_2nd_dist - 2nd:ary electron energy spectra ionisation of O2.
%
% Calling:
%  Xs = O2_e_2nd_dist(Es,E0in)
% Input:
%  Es   - energy of secondar electrons (eV), double array [1 x nEs]
%  E0in - energy of primary electrons (eV), double scalar
%  Eionization - ionization threshold (eV), double scalar so that
%         we can properly calculate the energy-distribution of
%         degrading electrons due to ionization to excited O2+ and
%         dissociative ionizations
%  s_or_c - spectra of secondary electrons or cascading electrons
% Output:
%  Xs - spectra of secondary electrons, double array [1 x nEs]
% 
% NOTE: Unnormalized! So that the user have to calculate the total
% ionization by normalizing Xs by its sum and scale with the
% ionizing cross-section, the N2 density and the electron flux.
% 
% Equation from Ittikawa et al. 1989
% 
% See also: N2_e_2nd_dist, O_e_2nd_dist

persistent Q E4Q Eionizations

if nargin < 5
  disp('Error with the O2_e_2nd_dist function, lacks the AURORA root directory argument')
end

if isempty(Q)  || ~all(E4Q(1:min(numel(Es),numel(E4Q)))==Es(1:min(numel(Es),numel(E4Q)))) || size(Es,2) > size(E4Q,2)
  try
    e_2nd_s_files = dir(fullfile(AURORA_root_directory,'E_cascadings','O2'));
    for i1 = 1:numel(e_2nd_s_files),
      if ~ e_2nd_s_files(i1).isdir
        try
          load(fullfile(AURORA_root_directory,...
                        'E_cascadings','O2',...
                        e_2nd_s_files(i1).name),...
               'E4Q')
%           if isequal(E4Q,Es)
            if all(E4Q(1:min(numel(Es),numel(E4Q)))==Es)
            % then we have found a match, so load
            fprintf('Loading cascading-matrices from file: %s\n',e_2nd_s_files(i1).name)
            load(fullfile(AURORA_root_directory,...
                         'E_cascadings','O2',...
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
  fprintf('Starting to calculate the requested cascading-matrices\n')
  end
end

E_hat = 15.2;
if isempty(Q) || ~all(E4Q(1:min(numel(Es),numel(E4Q)))==Es(1:min(numel(Es),numel(E4Q)))) || size(Es,2) > size(E4Q,2)
  Eionizations = [12.072
                  16.1
                  16.9
                  18.2
                  18.9
                  32.51];
  Q = zeros(numel(Es),numel(Es),numel(Eionization));
  E4Q = Es;
  dE = diff(Es);dE = dE([1:end,end]);
  disp('Pre-calculating all energy-degradations for e - O2-ionizing collisions.')
  disp('This will take a bit of time.')
  disp(['Starting at: ',datestr(now,'HH:MM:SS')])
  for i3 = numel(Eionizations):-1:1,
    iLim = find(Es>Eionizations(i3),1,'first');
    for i2 = iLim:numel(Es),
      iHalf = find(Es<(Es(i2)-Eionizations(i3))/2,1,'last');
      for i1 = iHalf:(i2-1),
        Eprime = [Es(i2),Es(i2)+dE(i2)];
        Edeg = [Es(i1),Es(i1)+dE(i1)];
        %funH = @(E1,E2) heaviside(E2-Eionizations(i3)-E1)./((E2-Eionizations(i3)-E1).^2+E_hat^2);
        fun = @(E1,E2) 1./((E2-Eionizations(i3)-E1).^2+E_hat^2);
        % fun = @(x,y) interp2(E1,E2,f,x,y);
        Emin = @(E) max(Eprime(1),E+Eionizations(i3));           
        % Emin = @(E) max(Eprime(1),E+Eworks(i3));
        Edegmax = min(Es(i1)+dE(i1),Es(i2)+dE(i2)-Eionizations(i3));
        if Edegmax > Es(i1)
          % Edges(i1) = {Emin,Edegmax};
          Q(i2,i1,i3) = integral2(fun,Es(i1),Edegmax,Emin,Eprime(end));
         % Q2(i2,i1,i3) = integral2(funH,Es(i1),Edegmax,Emin,Eprime(end));
         end
      end
    end
    disp(['Done with level ',num2str(i3),' at: ',datestr(now,'HH:MM:SS')])
  end
    
  % Save the results for future use
  save_filename = sprintf('CascadingSpecO2ionization_%s.mat',...
  datestr(now,'yyyymmdd-HHMMSS'));
  save(fullfile(AURORA_root_directory,...
  'E_cascadings','O2',...
  save_filename),...
  'Q',...
  'E4Q',...
  'Eionizations')
end
% save Q_both_O2.mat Q Q2
if isempty(Eionization)
  Eionization = 12.072;
end

if strcmp(s_or_c,'s') % Secondary spectra
  
  Xs = 1./( E_hat^2 + Es.^2 ).*( Es < ( E0in - Eionization )/2 );
  
else % Cascading spectra
  
  [DE,iLevel] = min(abs(Eionization-Eionizations));
  [DE2,iPrimary] = min(abs(E4Q-E0in));
  Xs = Q(iPrimary,1:numel(Es),iLevel);
  
end

%% Below here there be trash and garbage, only comments/BG 20190408
%% Equation and parameters from M.H. Rees [1989]
%Ebar = 17.4; %(eV) Shape parameter
%E_I = 12.2; % (eV) Ionization energy
%
%Xs = (1+Es/Ebar).^(-2.1)./tanh((E0in-E_I)/(2*Ebar))/Ebar;
%Xs(Es>E0in/2) = 0;
%Xs(~isfinite(Xs)) = 0;

% $$$ 
% $$$ 
% $$$ for i2 = iLim:numel(E),
% $$$   iHalf = find(E<(E(i2)-12.75)/2,1,'last');
% $$$   for i1 = iHalf:(i2+1),     
% $$$     Eprime = linspace(E(i2),E(i2)+dE(i2),10); 
% $$$     Efine = linspace(E(i1),E(i1)+dE(i1),10);
% $$$     [E1,E2] = meshgrid(Efine,Eprime);
% $$$     f = @(E1,E2) heaviside(E2-12.75-E1)./((E2-12.75-E1).^2+15^2);
% $$$     F = f(E1,E2); 
% $$$     Q(i1,i2) = trapz(E2(:,1),trapz(E1(1,:),F,2));
% $$$   end
% $$$   if rem(i2,10) == 0
% $$$     disp(sprintf('%d: %s',i2,datestr(now,'HH:MM:SS')))
% $$$   end
% $$$ end
% $$$ if isempty(Q)
% $$$   load CascadingSpecO2ionization.mat Q Eall Eworks
% $$$   Eworks = [12.072
% $$$             16.1
% $$$             16.9
% $$$             18.2
% $$$             18.9
% $$$             32.51];
% $$$ end
