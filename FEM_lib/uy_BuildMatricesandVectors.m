function [S,f] = uy_BuildMatricesandVectors(myCFD)
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
f		= zeros(n,1); % right-hand side vector

%
% Treatment of the internal (triangular) elements
%
topology = myCFD.Mesh.topology.element;
topologybnd = myCFD.Mesh.topology.boundary;

elmat = myCFD.Mesh.Elements;
elmatbnd = myCFD.Mesh.BndElements;

for el_index = 1:length(elmat(:,1)) % for all internal elements
	Selem = uy_GenerateElementMatrix(el_index,myCFD); % Selem	
    for ind1 = 1:topology
        for ind2 = 1:topology
            S(elmat(el_index,ind1),elmat(el_index,ind2))	= S(elmat(el_index,ind1),elmat(el_index,ind2)) + Selem(ind1,ind2);
        end
    end
    %source term for x momentum
    felem = uy_GenerateElementVector(el_index,myCFD); % felem
    for ind1 = 1:topology
        f(elmat(el_index,ind1)) = f(elmat(el_index,ind1)) + felem(ind1);
    end
end

% Next the boundary contributions

for el_index = 1:length(elmatbnd(:,1)) % for all boundary elements extension of mass matrix M and element vector f
	BMelem = uy_GenerateBoundaryElementMatrix(el_index,myCFD); % BMelem
    for ind1 = 1:topologybnd
        for ind2 = 1:topologybnd
            S(elmatbnd(el_index,ind1),elmatbnd(el_index,ind2)) = S(elmatbnd(el_index,ind1),elmatbnd(el_index,ind2)) + BMelem(ind1,ind2);
        end
    end
    
    % boundary source term x direction
	bfelem = uy_GenerateBoundaryElementVector(el_index,myCFD); % bfelem   
    for ind1 = 1:topologybnd
        f(elmatbnd(el_index,ind1)) = f(elmatbnd(el_index,ind1)) + bfelem(ind1);
    end
end

%% modify for fixed points
%identify fixed value boundaries
fixedvalues = false(size(myCFD.Mesh.PhysicalTag,2),1);
for i=1:size(myCFD.Mesh.PhysicalNames,2)
    boundarytag = myCFD.Mesh.PhysicalNames(i);
    if strcmp(myCFD.boundaries.Uy.(boundarytag).type,"Fixed value")
        fixedvalues = fixedvalues | ...
                      strcmp(myCFD.Mesh.PhysicalTag,myCFD.Mesh.PhysicalNames(i))';
    end
end

%init c with fixed values
c = zeros(length(fixedvalues),1);
index_fv = find(fixedvalues);
for i=1:nnz(fixedvalues)
    tag = myCFD.Mesh.PhysicalTag(index_fv(i));
    c(index_fv(i)) = myCFD.boundaries.Uy.(tag).value;
end

%extract collums that belong to fixed values and add coefficnents in the rows.
Sinhom = S(:,fixedvalues);
finhom = Sinhom*c(fixedvalues);

% subtract from rhsvector
f = f-finhom;


S(fixedvalues,fixedvalues) = speye(sum(fixedvalues));
f(fixedvalues) = c(fixedvalues);


%add under relaxation
alpha = myCFD.sim_settings.under_relax_fac;
U_old = myCFD.Solution.Uy;
S = S + (1-alpha)/alpha*diag(diag(S));
f = f + (1-alpha)/alpha*(diag(diag(S))*U_old);
end

