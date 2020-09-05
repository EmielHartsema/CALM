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

%create waitbar
bar = waitbar(0,'computing solution...','Name','calculating solution...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(bar,'canceling',0);
for outerloop = 1:myCFD.sim_settings.num_iter
%% Momentum predictor
% assemble matrix for ux
[M,f] = ux_BuildMatricesandVectors(myCFD);
%figure;spy(M);

% Solve system
myCFD.Solution.Ux = solvematrix_momentum_x(myCFD,M,f);
%pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.Ux)
%calculate residuals
myCFD.Residual.Ax = diag(M);
A = diag(diag(M));
myCFD.Residual.Hx = (M-A)*myCFD.Solution.Ux;

% assemble matrix for uy
[M,f] = uy_BuildMatricesandVectors(myCFD);
%figure;spy(M);

% Solve system
myCFD.Solution.Uy = solvematrix_momentum_y(myCFD,M,f);

%calculate residuals
myCFD.Residual.Ay = diag(M);
A = diag(diag(M));
myCFD.Residual.Hy = (M-A)*myCFD.Solution.Uy;

%% Pressure corrector
%assemble matrix
[M,f] = pres_BuildMatricesandVectors(myCFD);
%figure;spy(M);
% solve system
myCFD.solution.p = solvematrix_pressure(myCFD,M,f);

%% explicit calculation of momentum
[Ux,Uy] = explicitcalc(myCFD);
myCFD.Solution.Ux = Ux;
myCFD.Solution.Uy = Uy;

%waitbar controls
waitbar(outerloop/myCFD.sim_settings.num_iter,bar,'computing solution...');
% Check for clicked Cancel button
if getappdata(bar,'canceling')
    break
end
end%end of outer loop
delete(bar)
% Post Processing
myCFD.Solution.U = sqrt(myCFD.Solution.Ux.^2+myCFD.Solution.Uy.^2);

% Display solution
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.Ux)
