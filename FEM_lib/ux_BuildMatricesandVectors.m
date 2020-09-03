function [S,fx,fy] = BuildMatricesandVectors(myCFD)
%
% This routine constructs the large matrices and vector.
%
% The element matrices and vectors are also dealt with.
%
%
% First the internal element contributions
%
% First Initialisation of large discretisation matrix, right-hand side vector
n = size(myCFD.Mesh.Nodes,2);
S 		= sparse(n,n); % stiffness matrix
fx		= zeros(n,1); % right-hand side vector
fy		= zeros(n,1); % right-hand side vector

%
% Treatment of the internal (triangular) elements
%
topology = myCFD.Mesh.topology.element;
topologybnd = myCFD.Mesh.topology.boundary;

elmat = myCFD.Mesh.Elements;
elmatbnd = myCFD.Mesh.BndElements;

for el_index = 1:length(elmat(:,1)) % for all internal elements
	Selem = GenerateElementMatrix(el_index,myCFD); % Selem	
    for ind1 = 1:topology
        for ind2 = 1:topology
            S(elmat(el_index,ind1),elmat(el_index,ind2))	= S(elmat(el_index,ind1),elmat(el_index,ind2)) + Selem(ind1,ind2);
        end
    end
    %source term for x momentum
    felem = GenerateElementVector(el_index,myCFD,1); % felem
    for ind1 = 1:topology
        fx(elmat(el_index,ind1)) = fx(elmat(el_index,ind1)) + felem(ind1);
    end
    
    %source term for y momentum
    felem = GenerateElementVector(el_index,myCFD,2); % felem
    for ind1 = 1:topology
        fy(elmat(el_index,ind1)) = fy(elmat(el_index,ind1)) + felem(ind1);
    end
end

% Next the boundary contributions

for el_index = 1:length(elmatbnd(:,1)) % for all boundary elements extension of mass matrix M and element vector f
	BMelem = GenerateBoundaryElementMatrix(el_index,myCFD); % BMelem
    for ind1 = 1:topologybnd
        for ind2 = 1:topologybnd
            S(elmatbnd(el_index,ind1),elmatbnd(el_index,ind2)) = S(elmatbnd(el_index,ind1),elmatbnd(el_index,ind2)) + BMelem(ind1,ind2);
        end
    end
    
    % boundary source term x direction
	bfelem = GenerateBoundaryElementVector(el_index,myCFD,1); % bfelem   
    for ind1 = 1:topologybnd
        fx(elmatbnd(el_index,ind1)) = fx(elmatbnd(el_index,ind1)) + bfelem(ind1);
    end
    
    % boundary source term y direction
    bfelem = GenerateBoundaryElementVector(el_index,myCFD,2); % bfelem   
    for ind1 = 1:topologybnd
        fy(elmatbnd(el_index,ind1)) = fy(elmatbnd(el_index,ind1)) + bfelem(ind1);
    end
end

end
