function [data] = sod_analytic(time,XL)
%time  = 0.2;
%Initial conditions
gam  = 1.4;
rho4 = 1;
u4   = 0;
P4   = 1;
T4   = P4/rho4;
a4   = sqrt(gam*P4/rho4);

rho1 = 0.125;
u1   = 0;
P1   = 0.1;
T1   = P4/rho1;
a1   = sqrt(gam*P1/rho1);

%------------ Shock-tube calculation ------------
P4P1 = P4/P1;
func = @(Ms) (P4P1-(2*gam*Ms(1)^2-(gam-1))/(gam+1) ...
    *(1-(gam-1)/(gam+1)*a1/a4*(Ms(1)-1/Ms(1)))^(-2*gam/(gam-1)))^2;
[Ms,eps] = fminsearch(func,1);
val  = (2*gam*Ms^2-(gam-1))/(gam+1) ...
   *(1-(gam-1)/(gam+1)*a1/a4*(Ms-1/Ms))^(-2*gam/(gam-1));
rho2 = ((gam+1)*Ms^2)/((gam-1)*Ms^2+2)*rho1;
U2   = 2/(gam+1)*(Ms-1/Ms)*a1;
P2   = (2*gam*Ms^2-(gam-1))/(gam+1)*P1;
T2   = P2/rho2;
a2   = sqrt(gam*P2/rho2);
P3   = P2;
rho3 = (P3/P4)^(1/gam)*rho4;
U3   = U2;
a3   = sqrt(gam*P3/rho3);
T3   = P3/rho3;

Us   = Ms*a1;
Utail= U3-a3;
Uhead=-a4;
x1   = Us*time;
x2   = U2*time;
x3   = Utail*time;
x4   = Uhead*time;


nx   = 10000;
dx  = XL/nx;
xxx = -XL/2+dx/2:dx:XL/2-dx/2;
Pre = zeros(nx,1);
Tem = zeros(nx,1);
Den = zeros(nx,1);
Vel = zeros(nx,1);

for i = 1:1:nx
   if xxx(i)>x1
      Pre(i) = P1;
      Tem(i) = P1/rho1;
      Den(i) = rho1;
      Vel(i) = 0.0;
   elseif xxx(i)>x2
      Pre(i) = P2;
      Tem(i) = P2/rho2;
      Den(i) = rho2;
      Vel(i) = U2;
   elseif xxx(i)>x3
      Pre(i) = P3;
      Tem(i) = P3/rho3;
      Den(i) = rho3;
      Vel(i) = U3;
   elseif xxx(i)>x4
      aa = (x3-xxx(i))/(x3-x4);
      Tem(i) = T3+(T4-T3)*aa;
      Pre(i) = P4*(Tem(i)/T4)^(gam/(gam-1));
      Den(i) = (Pre(i)/P3)^(1/gam)*rho3;
      Vel(i) = U3*(1-aa);
   else
      Pre(i) = P4;
      Tem(i) = P4/rho4;
      Den(i) = rho4;
      Vel(i) = 0.0;
   end
end

data.x = xxx;
data.rho = Den;
data.u = Vel;
data.P = Pre;
data.T = Tem;


