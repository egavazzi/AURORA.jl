function [phfcnE,phfcnI] = phase_fcn_O2(E,theta)
% PHASE_FCN_N2 - 
%   


for iE = numel(E):-1:1,
  for i_theta = numel(theta):-1:1,
    DCS(i_theta,iE) = DCSO2(theta(i_theta),E(iE));
  end
end
i3 = find(E>=4.15,1,'first');
for iE = 1:(i3-1),
  DCS(:,iE) = DCS(:,i3);
end

phfcnE = DCS./repmat(sum(DCS,1),size(theta));

for iE = numel(E):-1:1,
  
  if E(iE)<=150
    Etmp = E(iE)*exp(0.025*(E(iE)-5.0));
  else
    Etmp = 50*E(iE);
  end
  for i_theta = numel(theta):-1:1,
    DCS(i_theta,iE) = DCSO2(theta(i_theta),E(iE));
  end
end

phfcnI = DCS./repmat(sum(DCS,1),size(theta));
