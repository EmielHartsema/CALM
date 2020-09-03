function [bfelem] = pres_GenerateBoundaryElementVector(el_index,myCFD)
%GENERATEBOUNDARYELEMENTVECTOR Summary of this function goes here
%   Detailed explanation goes here

topologybnd = myCFD.Mesh.topology.boundary;
elmatbnd = myCFD.Mesh.BndElements;
bndtag = myCFD.Mesh.PhysicalTag;

Ax = myCFD.Residual.Ax;
Ay = myCFD.Residual.Ay;
Hx = myCFD.Residual.Hx;
Hy = myCFD.Residual.Hy;
rAHx = (1./Ax).*Hx;
rAHy = (1./Ay).*Hy;

x = myCFD.Mesh.Nodes(1,:);
y = myCFD.Mesh.Nodes(2,:);
xc = zeros(1,topologybnd);
yc = zeros(1,topologybnd);
for index1 = 1:topologybnd
	xc(index1) = x(elmatbnd(el_index,index1));
	yc(index1) = y(elmatbnd(el_index,index1));
end;

lek = sqrt((xc(2)-xc(1))^2+(yc(2)-yc(1))^2);
Delta = lek*2;

%for index1 = 1:topologybnd
%		bfelem(index1) = 0;
%end
[nx,ny] = FindNormal(elmatbnd(el_index,:),myCFD.Mesh);

bfelem = zeros(1,topologybnd);
for index1 = 1:topologybnd
    globalindex1 = elmatbnd(el_index,index1);
    boundarytype = myCFD.boundaries.p.(bndtag(globalindex1)).type;
    if strcmp(boundarytype,"Fixed value")
        bfelem(index1) = 0;
    end
    if strcmp(boundarytype,"Zero gradient")
        sum_raH = 0;
        for i=1:topologybnd
            global_index = elmatbnd(el_index,i);
            sum_raH = sum_raH + ...
                                (rAHx(global_index)*nx+rAHy(global_index)*ny)*...
                                (1+double(global_index==i));
            
        end
        bfelem(index1) = sum_raH*abs(Delta)/6;
    end
end
end

