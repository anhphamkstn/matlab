function F=F3DCoils(x,B,f,Ra,Rb,d_isol,rho);

rb1 = x(1);		
NV1 = x(2);		
NV2 = x(3);
NV3 = x(4);
y(1:16)=x(5:20);

% d_cu         : diameter of copper wire. (set 2)
% d            : diameter of copper wire with isolation. (set 2)
% NH           : number of windings placed next to each other on a layer on the coil. (integer)
% NV1, NV2, NV3: number of layers of windings on top of each other for each coil. (integers or continuous)
% rb1          : inside radius smallest coil, coil 1. (continuous)
% d_isol       : isolation thickness between the coils. (constant)
% rho          : copper material constant. (constant)
% B,f          : strength and frequency of magnetic fields. (constants)
% Ra,Rb   		: configuration of op-amp. (constants)

% set 1 op-amp (implemented in noise.m)
% 		    op-amp
% y(1)    LT1028
% y(2)    LF411
% y(3)    OP27/37
% y(4)    OP97
% y(5)    muA741
% x(6)    TLE12027

% set 2 wire
%			d_cu(mm)
% y(7) 	0.05 
% y(8)	0.10
% y(9)	0.15
% y(10)	0.20
% y(11)	0.25
% y(12)	0.30
% y(13)	0.35
% y(14)	0.40
% y(15)	0.45
% y(16)0.50
% d = 1.25*d_cu

d_cu = 0.05e-3+1e-3*sum(y(7:16).*[0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45]);
d = 1.25 * d_cu;
NH=0.02/d;

% B = 1e-6; f=0.1;
% Ra = 100; Rb = 3*1e6;
% d_isol = 0.001;
% rho = 1.724*1e-8;

a = 3*d^2*pi*NH/4;
a1 = (2*rb1/d + (1 - sqrt(3)))/sqrt(3);
A1 = a*(a1^2*NV1 + a1*NV1*(NV1+1) + sum([1:NV1].^2));
l1 = pi*NH*d*sqrt(3)*(a1*NV1 + NV1*(NV1+1)/2);
routside1 = rb1 + (1 - sqrt(3))*d + NV1*sqrt(3)*d/2;

rb2 = routside1 + d_isol;
a2 = (2*rb2/d + (1 - sqrt(3)))/sqrt(3);
A2 = a*(a2^2*NV2 + a2*NV2*(NV2+1) + sum([1:NV2].^2));
l2 = pi*NH*d*sqrt(3)*(a2*NV2 + NV2*(NV2+1)/2);
routside2 = rb2 + (1 - sqrt(3))*d + NV2*sqrt(3)*d/2;

rb3 = routside2 + d_isol;
a3 = (2*rb3/d + (1 - sqrt(3)))/sqrt(3);
A3 = a*(a3^2*NV3 + a3*NV3*(NV3+1) + sum([1:NV3].^2));
l3 = pi*NH*d*sqrt(3)*(a3*NV3 + NV3*(NV3+1)/2);
routside3 = rb3 + (1 - sqrt(3))*d + NV3*sqrt(3)*d/2;

Ri = 4*rho*[l1,l2,l3]/(pi*d_cu^2);
Vi = 2*pi*f*B*[A1,A2,A3];

Vo =(1+Rb/Ra)*Vi;
Vnoise(1)=noise(Ra,Rb,Ri(1),f,y(1:6)');
Vnoise(2)=noise(Ra,Rb,Ri(2),f,y(1:6)');
Vnoise(3)=noise(Ra,Rb,Ri(3),f,y(1:6)');

Ratio = Vo./Vnoise;

F = -min(Ratio);

