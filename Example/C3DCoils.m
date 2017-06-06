function [C,Ceq]=C3DCoils(x,B,f,Ra,Rb,d_isol,rho);

rb1 = x(1);		
NV1 = x(2);		
NV2 = x(3);
NV3 = x(4);
y(1:16)=x(5:20);

d_cu = 0.05e-3+1e-3*sum(y(7:16).*[0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45]);
d = 1.25 * d_cu;
NH=0.02/d;

a = 3*d^2*pi*NH/4;
a1 = (2*rb1/d + (1 - sqrt(3)))/sqrt(3);
A1 = a*(a1^2*NV1 + a1*NV1*(NV1+1) + sum([1:NV1].^2));
routside1 = rb1 + (1 - sqrt(3))*d + NV1*sqrt(3)*d/2;

rb2 = routside1 + d_isol;
a2 = (2*rb2/d + (1 - sqrt(3)))/sqrt(3);
A2 = a*(a2^2*NV2 + a2*NV2*(NV2+1) + sum([1:NV2].^2));
routside2 = rb2 + (1 - sqrt(3))*d + NV2*sqrt(3)*d/2;

rb3 = routside2 + d_isol;
a3 = (2*rb3/d + (1 - sqrt(3)))/sqrt(3);
A3 = a*(a3^2*NV3 + a3*NV3*(NV3+1) + sum([1:NV3].^2));
routside3 = rb3 + (1 - sqrt(3))*d + NV3*sqrt(3)*d/2 + d_isol;

C=[
   30-A1         % < 0, minimal effective area
   A1-32         % < 0, maximal effective area
   30-A2         % < 0, minimal effective area
   A2-32         % < 0, maximal effective area
   30-A3         % < 0, minimal effective area
   A3-32         % < 0, maximal effective area
   routside3-0.08 % < 0, maximal outside radius
];

% Ceq=[
%    30 - A1
%    30 - A2
%    30 - A3
% ];

Ceq=[];