function ny_ee = ny_e_e(E,Te,ne,Ti)
% NY_E_E - equivalent electron-electron collision cross section.frequency
%   
% 
% Calling:
%   sigma_ee = e_e(E,Te,ne,Ti)
% Input:
%   E - electron energies [1xN] (eV)
%   Te - electron temperature (K)
%   ne - electron concentration (m^-3)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]
q_e     = 1.6021773e-19;            % elementary charge [C]
m_e     = 9.10939e-31;              % electron rest mass [kg]
kB	= 1.380662e-23;		    % Boltzmann constant [J/K]

%v = v_of_E(E);


% sigma_ee = 4*pi*b_0^2*(1+ln((b_0/l_D)^2))^(1/2)

% b_0 = q_e^2/(4*pi*Ep_0)/(m_e*v^2/2) = q_e^2/(4*pi*Ep_0)/(E)
if nargin == 3
  
  l_D = debye_length(ne,Te);
  
else
  
  l_D = debye_length(ne,Te,Ti);
  
end
v = v_of_E(E);

%b_0 = (q_e)./(E)/(4*pi*Ep_0);

a = (m_e./Te/2/kB).^(1/2);

if numel(l_D) == 1
  ny_ee = ne*q_e^4*log(4*pi*ne*l_D)./(2*pi*Ep_0^2*m_e*v.^3).*(erf(a*v)-(erf(a*v)-2/pi^.5*(a*v).*exp(-a^2*v^2))./(a*v).^2/2);
else
  a = repmat(a,size(E));
  v = repmat(v,size(Te));
  keyboard
  ny_ee = repmat(ne*q_e^4.*log(4*pi*ne.*l_D),size(E))./(2*pi*Ep_0^2*m_e*v.^3).*(erf(a.*v)-(erf(a.*v)-2/pi^.5*(a.*v).*exp(-a.^2.*v.^2))./(a.*v).^2/2);
end

