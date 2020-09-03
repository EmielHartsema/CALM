%{
-----------------------------------------------------
CALM - Computational Applicatoin for Liquid Mechanics
Author: Emiel Hartsema
Version: 0.0.1
Date: 16-4-2020
-----------------------------------------------------
%}

%Simulation parameters
n_it = 1;
% add path to libraries
addpath('readmsh');
addpath('FEM_lib')

% Read settings
%proj_folder = uigetdir;
proj_folder = 'C:\Users\Emiel\Documents\CALM\Testcase1_Cilinder_flow';
myCFD.transport_prop = jsondecode(fileread(strcat(proj_folder,'\transport_prop.json')));
myCFD.boundaries.Ux = jsondecode(fileread(strcat(proj_folder,'\boundaries\U.json')));
myCFD.boundaries.Uy = myCFD.boundaries.Ux;
myCFD.boundaries.p = jsondecode(fileread(strcat(proj_folder,'\boundaries\P.json')));
% Load mesh
myCFD.Mesh = importmsh(strcat(proj_folder,'\mesh.msh'));
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements')

%initialize solution vectors
myCFD.Solution.Ux = zeros(size(myCFD.Mesh.Nodes,2),1);
myCFD.Solution.Uy = zeros(size(myCFD.Mesh.Nodes,2),1);
myCFD.Solution.p = zeros(size(myCFD.Mesh.Nodes,2),1);

%create waitbar
bar = waitbar(0,'1','Name','calculating solution...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(bar,'canceling',0);
for outerloop = 1:n_it
%% Momentum predictor
% assemble matrix
[M,fx,fy] = BuildMatricesandVectors(myCFD);
%figure;spy(M);
% Solve system
myCFD.Solution.Ux = solvematrix_momentum(myCFD,M,fx);
myCFD.Solution.Uy = solvematrix_momentum(myCFD,M,fy);

%% Pressure corrector
myCFD.Residual.A = diag(M);
A = diag(diag(M));
myCFD.Residual.Hx = (M-A)*myCFD.Solution.Ux;
myCFD.Residual.Hy = (M-A)*myCFD.Solution.Uy;

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
waitbar(outerloop/n_it,bar,'computing solution');
% Check for clicked Cancel button
if getappdata(bar,'canceling')
    break
end
end%end of outer loop
delete(bar)
% Post Processing
pdeplot(myCFD.Mesh.Nodes,myCFD.Mesh.Elements','XYData',myCFD.Solution.Ux)


% Display solution

