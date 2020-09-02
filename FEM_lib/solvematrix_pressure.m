function [c] = solvematrix_pressure(myCFD,S,f)
%SOLVEMATRIX Summary of this function goes here
%   Detailed explanation goes here

%identify fixed value boundaries
fixedvalues = strcmp(myCFD.Mesh.PhysicalTag,"outflow");

% remove rows and collums that do not need to be evaluated
S2 = S(~fixedvalues,~fixedvalues);
f2 = f(~fixedvalues);

% solve the system for the non-fixed values
c = zeros(length(fixedvalues),1);
c(~fixedvalues) = S2\f2;

% add fixed value to solution vector
index_fv = find(fixedvalues);
for i=1:nnz(fixedvalues)
    tag = myCFD.Mesh.PhysicalTag(index_fv(i));
    c(index_fv(i)) = myCFD.boundaries.p.(tag).value;
end
end

