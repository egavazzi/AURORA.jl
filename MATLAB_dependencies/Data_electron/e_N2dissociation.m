function Xs = e_N2dissociation(Ep)
% Xs = e_N2dissociation(E)
% 
% e_N2dissociation - dissociation cross section (m^2)
% E electron energy (eV)

E = [20   23.6 27.8 35   45  55   65   80   96  110  125 148  170 ...
     195  245  295 2985.1 49183 55306];
s = [ 0.87 1.13 1.39 1.54 1.7 1.87 1.87 2.04 2.07 1.96 1.9 1.87 1.78 ...
      1.74 1.57 1.48 5.0198e-1/3  4.5953e-2/3 4.47e-2/3]*1e-20;%-16;


Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

% warning('Arbitrary correction of diss-cross-section')
CF = 0.3+0.7./(1+exp((Ep-120)/20));
Xs = Xs.*CF;

inans = find(Ep<20.6);
Xs(inans) = 0;

