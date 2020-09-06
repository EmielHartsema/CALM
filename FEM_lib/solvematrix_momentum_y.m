function [c] = solvematrix_momentum_y(myCFD,S,f)
%SOLVEMATRIX Summary of this function goes here
%   Detailed explanation goes here

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

% remove rows and collums that do not need to be evaluated
S2 = S(~fixedvalues,~fixedvalues);
f2 = f(~fixedvalues);

% solve the system for the non-fixed values
c(~fixedvalues) = S2\f2;
end

