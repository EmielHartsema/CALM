function [felem] = uy_GenerateElementVector(el_index,myCFD)
%GENERATEELEMENTVECTOR Summary of this function goes here
%   Detailed explanation goes here
% Module for element mass matrix for reactive term
%
% Output: felem  ====== vector of two components
%
% felem(1), felem(2) to be computed in this routine.
% extract data from myCFD
topology = myCFD.Mesh.topology.element;
elmat = myCFD.Mesh.Elements;
x = myCFD.Mesh.Nodes(1,:);
y = myCFD.Mesh.Nodes(2,:);

cp = myCFD.Solution.p;
xc = zeros(1,topology);
yc = zeros(1,topology);

for index1 = 1:topology
	xc(index1) = x(elmat(el_index,index1));
	yc(index1) = y(elmat(el_index,index1));
end

Delta = det([1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)]);

B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

%alpha = B_mat(1,1:3);
%beta  = B_mat(2,1:3);
gamma = B_mat(3,1:3);

sum_cpgamma = 0;
for index1 = 1:topology
	global_index = elmat(el_index,index1);
    sum_cpgamma = sum_cpgamma + cp(global_index)*gamma(index1);
end

felem(1:topology) = -sum_cpgamma*abs(Delta)/6;

end

