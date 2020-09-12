function [Ux,Uy] = explicitcalc(myCFD)
%EXPLICITCALC Summary of this function goes here
%   Detailed explanation goes here
Ax = diag(myCFD.Residual.Ax);
Ay = diag(myCFD.Residual.Ay);

Hx = myCFD.Residual.Hx;
Hy = myCFD.Residual.Hy;
%[gradpx,gradpy] = calc_gradp(myCFD);
[~,gradpx] = ux_BuildMatricesandVectors(myCFD);
[~,gradpy] = uy_BuildMatricesandVectors(myCFD);

% calculate velocity fields
Ux = (Ax\Hx)+(Ax\gradpx);
Uy = (Ay\Hy)+(Ay\gradpy);

% add under relaxation
alpha = myCFD.sim_settings.under_relax_fac;
ux_old = myCFD.Solution.Ux;
uy_old = myCFD.Solution.Uy;

Ux = alpha*Ux+(1-alpha)*ux_old;
Uy = alpha*Uy+(1-alpha)*uy_old;

end
