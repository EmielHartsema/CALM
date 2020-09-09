%{
-----------------------------------------------------
CALM - Computational Applicatoin for Liquid Mechanics
Author: Emiel Hartsema
Version: 0.0.1
Date: 16-4-2020
-----------------------------------------------------
%}

% add path to libraries
addpath('readmsh');
addpath('FEM_lib')

% Read settings
proj_folder = uigetdir;
%proj_folder = 'C:\Users\Emiel\Documents\CALM\Testcase1_Cilinder_flow';
myCFD.sim_settings = jsondecode(fileread(strcat(proj_folder,'\simulation_settings.json')));
myCFD.transport_prop = jsondecode(fileread(strcat(proj_folder,'\transport_prop.json')));
myCFD.boundaries.Ux = jsondecode(fileread(strcat(proj_folder,'\boundaries\Ux.json')));
myCFD.boundaries.Uy = jsondecode(fileread(strcat(proj_folder,'\boundaries\Uy.json')));
myCFD.boundaries.p = jsondecode(fileread(strcat(proj_folder,'\boundaries\P.json')));
% Load mesh
myCFD.Mesh = importmsh(strcat(proj_folder,'\mesh.msh'));
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements')

%initialize solution vectors
myCFD.Solution.Ux = zeros(size(myCFD.Mesh.Nodes,2),1);
myCFD.Solution.Uy = zeros(size(myCFD.Mesh.Nodes,2),1);
myCFD.Solution.p = zeros(size(myCFD.Mesh.Nodes,2),1);

%create figure to show residuals
residual_fig = figure('Name','Residuals');
for outerloop = 1:myCFD.sim_settings.num_iter
%% Momentum predictor
% assemble matrix for ux
[M,f] = ux_BuildMatricesandVectors(myCFD);
%figure;spy(M);

myCFD.Residual.rx(outerloop) = norm(M*myCFD.Solution.Ux-f,2);
% Solve system
myCFD.Solution.Ux = solvematrix_momentum_x(myCFD,M,f);
%calculate residuals
myCFD.Residual.Ax = diag(M);
A = diag(diag(M));
myCFD.Residual.Hx = A*myCFD.Solution.Ux-M*myCFD.Solution.Ux;


% assemble matrix for uy
[M,f] = uy_BuildMatricesandVectors(myCFD);
%figure;spy(M);
myCFD.Residual.ry(outerloop) = norm(M*myCFD.Solution.Uy-f,2);
% Solve system
myCFD.Solution.Uy = solvematrix_momentum_y(myCFD,M,f);

%calculate residuals
myCFD.Residual.Ay = diag(M);
A = diag(diag(M));
myCFD.Residual.Hy = A*myCFD.Solution.Uy - M*myCFD.Solution.Uy;

%% Pressure corrector
%assemble matrix
[M,f] = pres_BuildMatricesandVectors(myCFD);
%figure;spy(M);
% solve system
myCFD.Solution.p = solvematrix_pressure(myCFD,M,f);

%% explicit calculation of momentum
[Ux,Uy] = explicitcalc(myCFD);
myCFD.Solution.Ux = Ux;
myCFD.Solution.Uy = Uy;
disp_residual(residual_fig,myCFD)
end%end of outer loop
% Post Processing
myCFD.Solution.U = sqrt(myCFD.Solution.Ux.^2+myCFD.Solution.Uy.^2);

% Display solution
figure
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.Ux)
title('Ux')
figure
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.Uy)
title('Uy')
figure
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.p)
title('p')
