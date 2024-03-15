function [J_down,J_up,IeEdown,IeEup] = Ie2currents(t,h_atm,E,Ie_ztE,c_o_mu,idx_z)
% Ie2currents - calculate current from multistream electron fluxes
% 
% Calling:
%  [J_down,J_up,IeEdown,IeEup] = Ie2currents(t,h_atm,E,Ie_ztE,c_o_mu,idx_z)
% Input:
%  t,     - time-scale (s), double array [1 x n_t]
%  h_atm  - altitudes (km), double array [n_z x 1]
%  E      - energies (eV), double array [1 x n_E]
%  Ie_ztE - electron number-flux, double array [n_z*n_mu x n_t x n_E]
%  c_o_mu - centre-of-pitch-angle-stream, double array [1 x n_mu]
%           or [n_mu x 1]
%  idx_z  - index in altitude-dimension to plot, integer array [n_Z] in
%           range of <1, 2,... ,numel(h_atm) >
% Output:
%  J_down - downward current (A) as a function of time at selected
%           altitudes, double array [n_Z x n_t]
%  J_up   - upward current (A) as a function of time at selected
%           altitudes, double array [n_Z x n_t]
%  IeEdown - downward energy-flux (eV/m^2/s) as a function of time at selected
%            altitudes, double array [n_Z x n_t]
%  IeEup   - upward energy-flux (eV/m^2/s) as a function of time at selected
%            altitudes, double array [n_Z x n_t]


%   Copyright © 2019 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

nZ = numel(h_atm);
clims = [];
J_up = zeros(numel(idx_z),numel(t));
J_down = zeros(numel(idx_z),numel(t));
IeEup = zeros(numel(idx_z),numel(t));
IeEdown = zeros(numel(idx_z),numel(t));

q_e = 1.602176620898e-19;     % elementary charge [C]

if length(t) > 1
  for i_z = 1:numel(idx_z),
    for i_mu = 1:numel(c_o_mu),
      % keyboard
      if c_o_mu(i_mu) > 0
        J_up(i_z,:)  = J_up(i_z,:)  + abs(c_o_mu(i_mu))*sum(squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)),2)';
        IeEup(i_z,:) = IeEup(i_z,:) + abs(c_o_mu(i_mu))*sum(repmat(E,numel(t),1).*...
                                                            squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)),2)';
      else
        J_down(i_z,:)  = J_down(i_z,:) + abs(c_o_mu(i_mu))*sum(squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)),2)';
  
        IeEdown(i_z,:) = IeEdown(i_z,:) + abs(c_o_mu(i_mu))*sum(repmat(E,numel(t),1).*...
                                                                squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)),2)';
      end
    end
  end
elseif length(t) == 1
  for i_z = 1:numel(idx_z),
      for i_mu = 1:numel(c_o_mu),
        % keyboard
        if c_o_mu(i_mu) > 0
          J_up(i_z,:)  = J_up(i_z,:)  + abs(c_o_mu(i_mu))*sum(squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))';
          IeEup(i_z,:) = IeEup(i_z,:) + abs(c_o_mu(i_mu))*sum(repmat(E,numel(t),1)' .*...
                                                              squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))';
        else
          % disp(size(J_down))
          % disp(size(Ie_ztE))
          % disp(size((squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))))
          % disp(size(sum(repmat(E,numel(t),1).*squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))'));
          % disp(size(repmat(E, numel(t), 1)))
          J_down(i_z,:)  = J_down(i_z,:) + abs(c_o_mu(i_mu))*sum(squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))';    
          IeEdown(i_z,:) = IeEdown(i_z,:) + abs(c_o_mu(i_mu))*sum(repmat(E,numel(t),1)' .*...
                                                                  squeeze(Ie_ztE(idx_z(i_z)+(i_mu-1)*numel(h_atm),:,:)))';
        end
      end
  end
end






J_up   = q_e * J_up;
J_down = q_e * J_down;
