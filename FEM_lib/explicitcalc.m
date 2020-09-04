function [Ux,Uy] = explicitcalc(myCFD)
%EXPLICITCALC Summary of this function goes here
%   Detailed explanation goes here
Ax = diag(myCFD.Residual.Ax);
Ay = diag(myCFD.Residual.Ay);

Hx = myCFD.Residual.Hx;
Hy = myCFD.Residual.Hy;
cp = myCFD.Solution.p;
% calc grad P
gradpx = zeros(size(myCFD.Mesh.Nodes,2),1);
gradpy = zeros(size(myCFD.Mesh.Nodes,2),1);
for el_index = 1:size(myCFD.Mesh.Elements,1)
    topology = myCFD.Mesh.topology.element;
    elmat = myCFD.Mesh.Elements;
    x = myCFD.Mesh.Nodes(1,:);
    y = myCFD.Mesh.Nodes(2,:);

    xc = zeros(1,topology);
    yc = zeros(1,topology);
    for index1 = 1:topology
        xc(index1) = x(elmat(el_index,index1));
        yc(index1) = y(elmat(el_index,index1));
    end
    B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);
    
    beta  = B_mat(2,1:3);
    gamma = B_mat(3,1:3);
    
    for i=1:topology
        global_index = elmat(el_index,i);
        gradpx(global_index) = gradpx(global_index) + beta(i)*cp(global_index)/topology;
        gradpy(global_index) = gradpy(global_index) + gamma(i)*cp(global_index)/topology;
    end
end

Ux = Ax\Hx-Ax\gradpx;
Uy = Ay\Hy-Ay\gradpy;
