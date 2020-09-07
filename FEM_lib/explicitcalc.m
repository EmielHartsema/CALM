function [Ux,Uy] = explicitcalc(myCFD)
%EXPLICITCALC Summary of this function goes here
%   Detailed explanation goes here
Ax = diag(myCFD.Residual.Ax);
Ay = diag(myCFD.Residual.Ay);

Hx = myCFD.Residual.Hx;
Hy = myCFD.Residual.Hy;
[gradpx,gradpy] = calc_gradp(myCFD);

% calculate velocity fields
Ux = (Ax\Hx)-(Ax\gradpx);
Uy = (Ay\Hy)-(Ay\gradpy);
