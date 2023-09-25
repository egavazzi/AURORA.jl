function [DCS]=DCSN2(Ang,Energy)


CoeffA1=[-7.55013 6.4651 -4.61116 1.51264 -0.24804 0.020829 -0.7241e-3];
CoeffA2=[-121.79 161.31 -88.9543 25.0846 -3.79789  0.294225 -0.91716e-2];
CoeffA3=[-24.4482 11.6556 -0.109564 -1.50976 0.479943 -0.0597 0.263263e-2];
CoeffEt1=[-9.06218 11.1088 -6.1393 1.43937 -0.149179 0.570119e-2 0.0];
CoeffB1=[3.481694 -0.47699  0.0];
CoeffB2=[-42.1667 13.79589 -1.123838];
CoeffB3=[-35.49183 12.7418 -1.247576]; 
CoeffEt2=[1.590443 -0.207755 0.0];


Elog(1,1)=1.0;
Elog(1,2)=log(Energy);
  
  for i=3:7
   Elog(1,i)=Elog(1,i-1)*Elog(1,2);
  end
    
  if Energy <= 500.0
    B(1)=exp(sum(Elog.*CoeffA1));
  else
    B(1)=2.0-exp(sum(Elog(1:3).*CoeffB1));
  end
  if Energy <= 250.0
    B(2)=exp(sum(Elog.*CoeffA2));
  else
    B(2)=exp(sum(Elog(1:3).*CoeffB2));
  end
  if Energy <= 100.0
    B(3)=exp(sum(Elog.*CoeffA3));
  else
    B(3)=exp(sum(Elog(1:3).*CoeffB3));
  end
  if Energy <= 500.0
    Etta=exp(sum(Elog(1:6).*CoeffEt1(1:6)));
  else
    Etta=1.6038-exp(sum(Elog(1:2).*CoeffEt2(1:2)));
  end
  

T=Energy/510879.0;
Etta=Etta*6.221e-5/T/(T+2);
x=(1.-cos(Ang)+2*Etta);
DCS=B(1)./x.^2+B(2)./x+B(3).*x.^4;


