function [areafunc_elem_x,areafunc_elem_y] = GenerateAreafuncElem(el_index,myCFD)
%UNTITLED2 Summary of this function goes here
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

Delta = det([1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)]);

B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

%alpha = B_mat(1,1:3);
beta  = B_mat(2,1:3);
gamma = B_mat(3,1:3);

areafunc_elem_x = zeros(3);
areafunc_elem_y = zeros(3);
for index1 = 1:topology
	for index2 = 1:topology
        areafunc_elem_x(index1,index2) = beta(index2)*abs(Delta)/6;
        areafunc_elem_y(index1,index2) = gamma(index2)*abs(Delta)/6;
	end
end
end

