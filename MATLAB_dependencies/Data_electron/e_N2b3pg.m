function Xs = e_N2b3pg(Ep)
% Xs = e_N2b3pg(E)
% 
% e_N2b3pg - electron excitation cross section (m^2)
% E electron energy (eV)

E = [7  8 9   10 12.5 15   17   20 30 50];
s = [.2 2 4.5 25 32   22.5 18.5 13  8  3]*1e-22;%-18;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

for ie=1:length(Ep),
  if Ep(ie) > 7.35 & Ep(ie) < 15.52,
    cross_section(ie)=(1-7.35/Ep(ie))*exp(-1139.542+1583.892*log(Ep(ie))-844.6222*log(Ep(ie))^2+198.2095*log(Ep(ie))^3-17.29356*log(Ep(ie))^4);
  elseif Ep(ie) >= 15.52 & Ep(ie) < 43.1,
    cross_section(ie)=exp(-0.9660707-18.27795*log(Ep(ie))-5.686039*log(Ep(ie))^2+4.259477*log(Ep(ie))^3-0.5800993*log(Ep(ie))^4);
  elseif Ep(ie) >= 43.1 & Ep(ie) < 100.0,
    cross_section(ie)=exp(-1087.040+1098.993*log(Ep(ie))-429.7112*log(Ep(ie))^2+74.32909*log(Ep(ie))^3-4.807258*log(Ep(ie))^4);    
  elseif Ep(ie) >=100,
    cross_section(ie)=6.148782e-13/Ep(ie)^3;
  else
    cross_section(ie)=0;
  end
end

Xs(Ep>95) = cross_section(Ep>95)/1e4;

inans = find(Ep<7.3532);
Xs(inans) = 0;
