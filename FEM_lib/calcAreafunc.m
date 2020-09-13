function [areafunc_x,areafunc_y] = calcAreafunc(myCFD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n        = size(myCFD.Mesh.Nodes,2);
areafunc_x = sparse(n,n); % stiffness matrix
areafunc_y = sparse(n,n); % stiffness matrix

topology = myCFD.Mesh.topology.element;
elmat = myCFD.Mesh.Elements;

for el_index = 1:length(elmat(:,1)) % for all internal elements
	[areafunc_elem_x,areafunc_elem_y] = GenerateAreafuncElem(el_index,myCFD); % Selem	
    for ind1 = 1:topology
        for ind2 = 1:topology
            areafunc_x(elmat(el_index,ind1),elmat(el_index,ind2))	= areafunc_x(elmat(el_index,ind1),elmat(el_index,ind2)) + areafunc_elem_x(ind1,ind2);
            areafunc_y(elmat(el_index,ind1),elmat(el_index,ind2))	= areafunc_y(elmat(el_index,ind1),elmat(el_index,ind2)) + areafunc_elem_y(ind1,ind2);
        end
    end
end

end

