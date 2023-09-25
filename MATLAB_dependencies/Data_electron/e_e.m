function sigma_ee = e_e(E,Te,ne)
% E_E - equivalent electron-electron collision cross section.
%   Calculates the Coulomb cross section (M. Mitchner and
%   C. H. Kruger, Jr, Partially ionized gases, 1973, pp 54-62)
% 
% Calling:
%   sigma_ee = e_e(E,Te,ne)
% Input:
%   E - electron energies [1xN] (eV)
%   Te - electron temperature (K)
%   ne - electron concentration (m^-3)


Ep_0    = 1e-9/36/pi;               % Permittivity [As/Vm]
q_e     = 1.6021773e-19;            % elementary charge [C]
m_e     = 9.10939e-31;              % electron rest mass [kg]

%v = v_of_E(E);


% sigma_ee = 4*pi*b_0^2*(1+ln((b_0/l_D)^2))^(1/2)

% b_0 = q_e^2/(4*pi*Ep_0)/(m_e*v^2/2) = q_e^2/(4*pi*Ep_0)/(E)
l_D = debye_length(ne,Te);

b_0 = (q_e)./(E)/(4*pi*Ep_0);


if numel(l_D) == 1
  sigma_ee = 4*pi*b_0.^2.*(1+log((b_0/l_D).^2)).^(1/2);
else
  sigma_ee = 4*pi*repmat(b_0,size(l_D)).^2.*(1+log((repmat(b_0,size(l_D))./repmat(l_D,size(b_0))).^2)).^(1/2);
end
