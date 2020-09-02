function [bfelem] = pres_GenerateBoundaryElementVector(el_index,myCFD)
%GENERATEBOUNDARYELEMENTVECTOR Summary of this function goes here
%   Detailed explanation goes here

topologybnd = myCFD.Mesh.topology.boundary;
elmatbnd = myCFD.Mesh.BndElements;
bndtag = myCFD.Mesh.PhysicalTag;

%for index1 = 1:topologybnd
%	xc(index1) = x(elmatbnd(i,index1));
%	yc(index1) = y(elmatbnd(i,index1));
%end;

%lek = sqrt((xc(2)-xc(1))^2+(yc(2)-yc(1))^2);

%for index1 = 1:topologybnd
%		bfelem(index1) = 0;
%end
bfelem = zeros(1,topologybnd);
for index1 = 1:topologybnd
    globalindex1 = elmatbnd(el_index,index1);
    boundarytype = myCFD.boundaries.p.(bndtag(globalindex1));
    if strcmp(boundarytype,"Fixed value")
        bfelem(index1) = 0;
    end
    if strcmp(boundarytype,"Zero gradient")
        sum_raH = 0;
        for i=1:topologybnd
            global_index = elmatbnd(el_index,i);
            sum_raH = sum_raH + 1/A(global_index)*...
                                (Hx(global_index)*nx+Hy(global_index)*ny)*...
                                (1+double(global_index==i));
            
        end
        bfelem(index1) = sum_raH*abs(Delta)/6;
    end
end
end

