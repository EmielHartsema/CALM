function [Selem] = pres_GenerateElementMatrix(el_index,myCFD)

% Module for element mass matrix for reactive term
%
% Output: Selem  ====== 2 x 2 matrix
%
% Selem(1,1), Selem(1,2), Selem(1,3), Selem(2,1), Selem(2.2), Selem(2,3),
%
% Selem(3,1), Selem(3,2), Selem(3,3) to be computed in this routine.
%
% elmat(i,1), elmat(i,2), elmat(i,3) give the index numbers of the vertices of element i
%
% x(elmat(i,j)), y(elmat(i,j)) give the coordinates of the vertices
%
% i = index number of element, imported from AssemblyStepStiffnessMatrix.m
%
% Selem(index1,index2) = (grad phi(elmat(i,index1)),grad phi(i,index2))
%
%extract data from Mesh
topology = myCFD.Mesh.topology.element;
elmat = myCFD.Mesh.Elements;
x = myCFD.Mesh.Nodes(1,:);
y = myCFD.Mesh.Nodes(2,:);
Ax = myCFD.Residual.Ax;
Ay = myCFD.Residual.Ay;


xc = zeros(1,topology);
yc = zeros(1,topology);
for index1 = 1:topology
	xc(index1) = x(elmat(el_index,index1));
	yc(index1) = y(elmat(el_index,index1));
end

Delta = det([1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)]);

B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

%alpha = B_mat(1,1:3);
beta  = B_mat(2,1:3);
gamma = B_mat(3,1:3);

Selem = zeros(3);
for index1 = 1:topology
	for index2 = 1:topology
        sum_ax = 0;
        sum_ay = 0;
        for i=1:3
            globalindex = elmat(el_index,i);
            sum_ax = sum_ax + 1/Ax(globalindex);
            sum_ay = sum_ay + 1/Ay(globalindex);
        end
        Selem(index1,index2) = (beta(index1)*beta(index2)*sum_ax+...
                                gamma(index1)*gamma(index2)*sum_ay)*...
                                abs(Delta)/6;
	end
end

end

