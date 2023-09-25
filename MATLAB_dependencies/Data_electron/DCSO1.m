function [DCS]=DCSO1(alfa,Energy)

CoeffA1=[-12.61 7.773 -3.922 0.825 -0.05546 0.0 0.0];
CoeffA2=[-20.30 16.37 -6.604 1.218 -0.08263 0.0 0.0];
CoeffA3=[-16.22 4.452  -0.3577  0.0  0.0  0.0  0.0];
CoeffEt1=[-4.008 1.017 -0.06066 0.0 0.0 0.0 0.0];
CoeffB1=[11.44 -1.986 0.0];
CoeffB3=[-19.04 5.805 -0.5186];
CoeffEt2=[3.291 -0.7007 0.0];

Elog(1)=1;
Elog(2)=log(Energy);
  
  for i=3:7
    Elog(i)=Elog(i-1)*Elog(2);
  end
  
  if Energy <= 500
    B(1)=exp(sum(Elog(1:5).*CoeffA1(1:5)));
  else
    B(1)=1.-exp(sum(Elog(1:2).*CoeffB1(1:2)));
  end
  
 B(2)=exp(sum(Elog(1:5).*CoeffA2(1:5)));
  
  if Energy <= 100
    B(3)=exp(sum(Elog(1:5).*CoeffA3(1:5)));
  else
    B(3)=exp(sum(Elog(1:3).*CoeffB3(1:3)));
  end
  
  if Energy <= 500
    Etta=exp(sum(Elog(1:3).*CoeffEt1(1:3))); 
  else
	Etta=1.3198-exp(sum(Elog(1:2).*CoeffEt2(1:2)));
  end
  
 T=Energy/510879.;
 Etta=Etta*6.8e-5/T/(T+2);
  
  if Energy <= 500
    x=4.4549-0.003114*Energy-0.4663*Elog(2);
	  x=exp(x);
	  x2=x.*x;
    FunAlfa=alfa./(x2-8100.).*(x2-8100+x*180.*(1.-alfa/3.14159265));
  else
    FunAlfa=alfa;
  end
  x=(1.-cos(FunAlfa)+2.*Etta);
  DCS=B(1)./x.^2+B(2)./x.^3+B(3).*x.^4;


