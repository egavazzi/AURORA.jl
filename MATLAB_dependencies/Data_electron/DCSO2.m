function [DCS]=DCSO2(alfa,Energy)


CoeffA1=[-59.543 79.8  -49.7886 16.3775 -2.90625 0.26045 -0.9109e-2];
CoeffA2=[-20.3557 33.1387 -30.4416 14.6647 -3.7744 0.489333 -0.024934];
CoeffA3=[-15.861 4.7481 -0.432 0.0 0.0 0.0 0.0];
CoeffEt1=[0.02369 0.5232 -0.585 0.240122 -0.041798 0.2622e-2 0.0];
CoeffB1=[3.381 -0.4622 0.0];
CoeffB2=[-121.6 41.0 -3.447];
CoeffB3=[-86.19 33.04 -3.259]; 
CoeffEt2=[2.308 -0.3645 0.0];
Elog(1)=1;
Elog(2)=log(Energy);
  
  for i=3:7
   Elog(i)=Elog(i-1)*Elog(2);
  end
  
  if Energy <= 500
    B(1)=exp(sum(Elog(:).*CoeffA1(:)));
  else
    B(1)=2.-exp(sum(Elog(1:2).*CoeffB1(1:2)));
  end
  
  if Energy <= 250
    B(2)=exp(sum(Elog(:).*CoeffA2(:)));
  else
    B(2)=exp(sum(Elog(1:3).*CoeffB2(1:3)));
  end
  
  if Energy <= 100
    B(3)=exp(sum(Elog(1:3).*CoeffA3(1:3)));
  else
    B(3)=exp(sum(Elog(1:3).*CoeffB3(1:3)));
  end
  
  if Energy <= 500
	Etta=(sum(Elog(1:6).*CoeffEt1(1:6)));
  else
	Etta=1.3141-exp(sum(Elog(1:2).*CoeffEt2(1:2)));
  end
T=Energy/510879;
Etta=Etta*6.8e-5/T/(T+2);
x=(1.-cos(alfa)+2.*Etta);
DCS=B(1)./x.^2+B(2)./x+B(3).*x.^4;


