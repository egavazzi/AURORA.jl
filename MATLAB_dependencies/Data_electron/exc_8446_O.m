function sigma = exc_8446_O(E)
% Xs = exc_8446(E)
% 
% exc_8446 - emission cross section (m^2) for 8446 AA
% from atomic oxygen. E electron energy (eV)
%
% Electron-impact excitation of OI 3p 3P from ground state
% in the excitation of O
% From Itikawa and Ichimura, J Phys Chem Ref Data 1990
%
%	Given by the generous Bjorn, March 2009
% parent       O
% products     0    (emission)
% threshold    10.99   (Itikawa) Confirmed/BG-20191016
% units        eV
% valid from   0.0   Remove this eventually - or replace with other comments ...
% valid to     0.0
% points       34

EnX = [ 1.1000000e+01   2.5781459e-01
        1.2059649e+01   3.6743271e-01
        1.3221375e+01   4.8399335e-01
        1.4495013e+01   5.9391952e-01
        1.5891342e+01   6.8434932e-01
        1.7422182e+01   7.4632073e-01
        1.9100490e+01   7.7643659e-01
        2.0940473e+01   7.6493066e-01
        2.2957704e+01   6.7097864e-01
        2.5169259e+01   5.5059069e-01
        2.7593856e+01   4.5203505e-01
        3.0252019e+01   3.9695472e-01
        3.3166248e+01   3.6965228e-01
        3.6361209e+01   3.5009898e-01
        3.9863946e+01   3.3411540e-01
        4.3704108e+01   3.1832706e-01
        4.7914199e+01   2.9997511e-01
        5.2529855e+01   2.7743037e-01
        5.7590145e+01   2.5248291e-01
        6.3137901e+01   2.3229768e-01
        6.9220083e+01   2.1658731e-01
        7.5888171e+01   2.0184430e-01
        8.3198607e+01   1.8801995e-01
        9.1213271e+01   1.7506665e-01
        1.0000000e+02   1.6293805e-01
        2.1544347e+02   8.8234000e-02
        4.6415888e+02   4.6805251e-02
        1.0000000e+03   2.4440707e-02
        2.1544347e+03   1.2604857e-02
        4.6415888e+03   6.4357220e-03
        1.0000000e+04   3.2587610e-03
        2.1544347e+04   1.6386314e-03
        4.6415888e+04   8.1909190e-04
        1.0000000e+05   4.0734512e-04];
        
sigma = exp(interp1(log(EnX(:,1)),log(EnX(:,2)),log(E),'pchip'));
sigma(E<10.99) = 0;
sigma = sigma*1e-21;
