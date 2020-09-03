function [bfelem] = ux_GenerateBoundaryElementVector(~,myCFD)
%GENERATEBOUNDARYELEMENTVECTOR Summary of this function goes here
%   Detailed explanation goes here

topologybnd = myCFD.Mesh.topology.boundary;

%for index1 = 1:topologybnd
%	xc(index1) = x(elmatbnd(i,index1));
%	yc(index1) = y(elmatbnd(i,index1));
%end;

%lek = sqrt((xc(2)-xc(1))^2+(yc(2)-yc(1))^2);

%for index1 = 1:topologybnd
%		bfelem(index1) = 0;
%end
bfelem(1:topologybnd) = zeros(1,topologybnd);
end

