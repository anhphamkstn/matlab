function Eno=opamp(R1,R2,Rs,f,set);
% Eno=opamp(R1,R2,Rs,f,set)
% Nois by different op-amps.
%
%                 En [V/sqrt(Hz)fcv[Hz]In [A/sqrt(Hz)] fci[Hz] 
% 1    LT1028  	 0.8e-9       3.5    1.0e-12         250     
% 2    LF411   	25.0e-9      60      1.0e-14         500     
% 3    OP27/37 	 3.0e-9       2.7    0.4e-12         140     
% 4    OP97    	18.0e-9       2      8.0e-15         300     
% 5    muA741  	10.0e-9      20      0.5e-12         150     
% 6    TLE12027 	 3.0e-9       8      0.7e-12          20     

Av = 1 + R2/R1;
T=300;
k=1.38e-23;
Et1=sqrt(4*k*T*R1);
Et2=sqrt(4*k*T*R2);
Ets=sqrt(4*k*T*Rs);

ar=[ 			0.8e-9   3.5 1.0e-12 250		%	LT1028
            25.0e-9  60   1.0e-14 500		%	LF411
             3.0e-9   2.7 0.4e-12 140		%	OP27/37
            18.0e-9   2   8.0e-15 300		%	OP97
            10.0e-9  20   0.5e-12 150		%  muA741
            3.0e-9   8   0.7e-12  20		%	TLE12027
         ];
ar1=sum(set.*ar(:,1));
ar2=sum(set.*ar(:,2));
ar3=sum(set.*ar(:,3));
ar4=sum(set.*ar(:,4));

En=sqrt(ar1^2 *(1 + ar2/f));
In=sqrt(ar3^2 *(1 + ar4/f));
Eni=sqrt(Ets^2 + En^2 + (Et2^2)*(R1/(R1+R2))^2 + (Et1^2)*(R2/(R1+R2))^2 +....
   (In*R1*R2/(R1 + R2))^2 + (In^2)*Rs^2);
Eno = Av*Eni;

