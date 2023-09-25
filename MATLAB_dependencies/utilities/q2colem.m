function I_lambda = q2colem(t,h,Q,A,tau)
% Q2COLEM - Integrate volume excitation-rates to column emission
%  Q2COLEM integrates altitude-time varying volume excitation-rates
%  to column intensities, taking the speed-of-light into
%  account.
% 
% Calling:
%  I_lambda = q2colem(t,h,Q,[A],[tau])
% Input:
%  t   - time (s), double array [1 x n_t]
%  h   - distance along magnetic field (m), double array [n_z x 1]
%  Q   - Volume excitation-rate (/m^3/s), double array [n_z x n_t]
%  A   - Einstein coefficiend (1/s), double scalar
%  tau - effective lifetime (s), double array [n_z x 1]
% 
% This function takes into account the time-delay between light
% emitted at different altitudes. For auroral photons emitted
% downwards from 100 and 200 km altitude the photons emitted from
% 200 km will arrive at the detector 100e3/3e8 = 0.333 ms later at
% the detector. This is a small time-shift but with the phase-shifts
% between auroral emissions varying at ~10 Hz this time-difference
% are getting close to the corresponding time-shifts.


%  Copyright © Bjorn Gustavsson 20190117, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


c0	= 2.99792458e8;		    % Speed of light [m/s]


if nargin < 4 || isempty(A)
  A = 1;
end
if nargin > 4 & ~isempty(tau)
  Q = Q.*repmat(tau,1,size(Q,2));
end

[t2,h2] = meshgrid(t,h);
I = interp2(t,h,Q*A,t2-(h2-h(1))/c0,h2,'linear',0);
I_lambda = trapz(h,I);
