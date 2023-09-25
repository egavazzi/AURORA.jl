function Q = exc_tz_of_Ie_ztE(t,h_atm,E,Ie_tzE,n,Xs)
% EXC_TZ_OF_IE_ZTE - volume excitation rate altitude-time variation
% from altitude-time-energy varying electron fluxes, height-varying
% neutral densities and excitation cross-section.
%   
% Calling:
%   Q = exc_tz_of_Ie_ztE(t,h_atm,E,Ie_tzE,n,Xs)
% Input:
%  t      - time (s), double array [1 x n_t]
%  h_atm  - height (m), double array [n_z x 1]
%  E      - energy (eV), double array [1 x n_E]
%  Ie_tzE - electron particle flux (/m^2/s/Delta_E), double array
%           [n_z*n_beams x n_t x n_E]
%  n      - density of exciteable atmospheric constituent (/m^3),
%           dobule array [n_z x 1]
%  Xs     - excitation cross-section (/m^2), double array [1 x n_E]
% Output:
%  Q      - volume excitation rate (/m^3/s), time-dependent
%           variation of volume excitation-rate, double array
%           [n_z x n_t]
% Example:
%  TBD
% This function calculates the excitation rates from
% electron-fluxes produced by the time-dependent multistream
% electron transport-code. This function automatically determines
% the number of streams used from the size-ratio
% size(Ie_tzE,1)/size(h,1) - since the multi-stream code stacks the
% fluxes in all streams in the first dimension.

% Copyright © B. Gustavsson 20190125
%  This is free software, licensed under GNU GPL version 2 or later

n_beams = size(Ie_tzE,1)/numel(h_atm);
for i_t = numel(t):-1:1,
  
  I_current = squeeze(Ie_tzE(:,i_t,:));
  
  n_x_sigma = n*Xs(end,:);
  Ie = 0;
  for iB = 1:n_beams,
    Ie = Ie + I_current((1:numel(h_atm))+(iB-1)*numel(h_atm),:);% + ...
  end
  Q(:,i_t) = sum(n_x_sigma.*Ie,2);
  
end
