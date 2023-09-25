function b_s_r = back_scattering_42stream(species,Ei)
% BACK_SCATTERING_42STREAM - Back scatter ratios for electron 2-stream transport
%  The back-scattering ratios used are from "The auroral
%  6300~{{\AA}} emisson: Observation and modelling ,S. C. Solomon,
%  P. B. Hays, V. J. Abreu, J. Geophys Res, Vol 93, No A9,
%  9867-9882, 1988,
%
% Calling:
%   b_s_r = back_scattering_42stream(species,Ei);
% Input:
%   species - Atmospheric species, ['N2' | 'O2' | 'O']
%   Ei      - input energy (eV) [1 x N]
%   B_S_R   - back-scatter ratio [2 x N] with B_S_R(1,:) being the
%             elastic and B_S_R(2,:) the inelastic back-scattering
%             ratios.

% Copyright B. Gustavsson 20080512
% Contact bjorn@irf.se

load back-scatter-ratios-solomon.mat E_Oe E_O2e E_N2e ...
                                     S_Oe S_O2e S_N2e ...
                                     E_Oi E_O2i E_N2i ...
                                     S_Oi S_O2i S_N2i


switch species
 case 'N2'
  bscre = exp(interp1(E_N2e,log(1./S_N2e),Ei,'pchip','extrap'));
  bscre(Ei>E_N2e(end)) = 1/S_N2e(end).*(E_N2e(end)./Ei(Ei>E_N2e(end))).^(1/1);
  bscri = exp(interp1(E_N2i,log(1./S_N2i),Ei,'pchip','extrap'));
  bscri(Ei>E_N2i(end)) = 1/S_N2i(end).*(E_N2i(end)./Ei(Ei>E_N2i(end)));
  b_s_r = [bscre;bscri];
 case 'O2'
  bscre = exp(interp1(E_O2e,log(1./S_O2e),Ei,'pchip','extrap'));
  bscre(Ei>E_O2e(end)) = 1/S_O2e(end).*(E_O2e(end)./Ei(Ei>E_O2e(end))).^(1/1);
  bscri = exp(interp1(E_O2i,log(1./S_O2i),Ei,'pchip','extrap'));
  bscri(Ei>E_O2i(end)) = 1/S_O2i(end).*(E_O2i(end)./Ei(Ei>E_O2i(end)));
  b_s_r = [bscre;bscri];
 case 'O'
  bscre = exp(interp1(E_Oe,log(1./S_Oe),Ei,'pchip','extrap'));
  bscre(Ei>E_Oe(end)) = 1/S_Oe(end).*(E_Oe(end)./Ei(Ei>E_Oe(end))).^(1/1);
  bscri = exp(interp1(E_Oi,log(1./S_Oi),Ei,'pchip','extrap'));
  bscri(Ei>E_Oi(end)) = 1/S_Oi(end).*(E_Oi(end)./Ei(Ei>E_Oi(end)));
  b_s_r = [bscre;bscri];
 otherwise
  warning(['Unknown species: ',species])
  b_s_r = zeros([2,length(Ei)]);
end
