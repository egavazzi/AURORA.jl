function [phfcnE,phfcnI] = phase_fcn_N2ei2I(E,theta)
% PHASE_FCN_N2 - 
%   

persistent E0 ph2nd phI

for iE = numel(E):-1:1,
  for i_theta = numel(theta):-1:1,
    DCS(i_theta,iE) = DCSN2(theta(i_theta),E(iE));
  end
end
i3 = find(E>=4.796,1,'first');
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
    DCS(i_theta,iE) = DCSN2(theta(i_theta),E(iE));
  end
end

phfcnI = DCS./repmat(sum(DCS,1),size(theta));

% $$$ for iE0 = 1:7,
% $$$   D = load(N2files(iE0).name);
% $$$   phiD = D(1,2:end-1)*pi/180;
% $$$   E_D = D(2:end,1);
% $$$   E_peakN2(iE0) = max(E_D);
% $$$   sigmaD = D(2:end,2:end-1);
% $$$   sigma_totD = D(2:end,end);
% $$$   sigma_2nde_N2(iE0,1:numel(sigma_totD)) = sigma_totD;
% $$$   sigma2EphiC = ddsigma_2nde(phiD,E_D,sigmaD,sigma_totD,721);
% $$$   ph_2nde_N(iE0,1:numel(E_D),:) = sigma2EphiC./repmat(sum(sigma2EphiC,2),[1,size(sigma2EphiC,2)]);
% $$$   scrollsubplot(3,2,1+2*(iE0-1))
% $$$   surf(phiD*180/pi,(E_D),log10(sigmaD)),view(0,90),shading flat     
% $$$   set(gca,'yscale','log')
% $$$   cx = caxis;
% $$$   axis([0 180,4 10^2.6])
% $$$   scrollsubplot(3,2,2*iE0)      
% $$$   surf(phi721*180/pi,(E_D),log10(medfilt2(sigma2EphiC,[3 3]))),shading flat     
% $$$   set(gca,'yscale','log')
% $$$   caxis(cx)
% $$$   view(0,90)
% $$$   axis([0 180,4 10^2.6]) 
% $$$ end
