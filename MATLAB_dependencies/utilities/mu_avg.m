function c_o_mu = mu_avg(mu_lims)
% MU_AVG - pitch-angle averages between limits of cosine-of-pitch-angles
%   for isotropically distributed fluxes within each pitch-angle
%   limits, i.e. flux weighted with sin(theta).
%
% Calling:
%   c_o_mu = mu_avg(mu_lims),
% Input:
%   mu_lims - pitch-angle limits for beams, double array [1 x n_beams+1]
% Output:
%   c_o_mu - pitch-angle averages, or centre-of-mu, double array [1 x n_beams+1]
% 
% Example:
%   mu_mean = mu_avg([-1 -cos(pi/4) 0 cos(pi/4) 1]);
%   disp(mu_mean)
%   disp(acos(mu_mean)*180/pi)
%   mu_mean = mu_avg([-1 0 1]);
%   disp(mu_mean)
%   disp(acos(mu_mean)*180/pi)

% Copyright © 20180430 Björn Gustavsson, bjorn.gustavsson@uit.no
% This is free software GPL copyright license 2.0 or later applies

theta_lims = acos(mu_lims);

c_o_mu = zeros(1,(numel(mu_lims)-1));
for i_th = 1:(numel(mu_lims)-1),
  c_o_mu(i_th) = integral(@(theta) cos(theta).*sin(theta),...
                          theta_lims(i_th),...
                          theta_lims(i_th+1)) / ...
      integral(@(theta) sin(theta),...
               theta_lims(i_th),...
               theta_lims(i_th+1));
  
end
