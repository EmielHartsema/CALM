function [BMelem] = uy_GenerateBoundaryElementMatrix(~,myCFD)
%GENERATEBOUNDARYELEMENTMATRIX Summary of this function goes here
%   Detailed explanation goes here
clear xc
clear yc
clear BMelem

topologybnd = myCFD.Mesh.topology.boundary;
% elmatbnd = myCFD.Mesh.BndElements;
% x = myCFD.Mesh.Nodes(1,:);
% y = myCFD.Mesh.Nodes(2,:);
% 
% nu = myCFD.transport_prop.nu;
% bndtag = myCFD.Mesh.PhysicalTag;
% xc = zeros(1,topologybnd);
% yc = zeros(1,topologybnd);
% for index1=1:topologybnd
% 	xc(index1) = x(elmatbnd(el_index,index1));
% 	yc(index1) = y(elmatbnd(el_index,index1));
% end
% 
% %lek = sqrt((xc(2)-xc(1))^2 + (yc(2)-yc(1))^2);
% [nx,ny] = FindNormal(elmatbnd(el_index,:),myCFD.Mesh);
BMelem = zeros(topologybnd);
% for index1=1:topologybnd
%     globalindex1 = elmatbnd(el_index,index1);
%     for index2 = 1:topologybnd
%         boundarytype = myCFD.boundaries.Uy.(bndtag(globalindex1)).type;
%         
%         %BMelem(index1,index2) = -nu*(beta(index2)*nx+gamma(index2)*ny)*abs(Delta)/2; %Transfer_rate_Coeff*lek/6*(1+double(eq(index1,index2)));
%         if strcmp(boundarytype,"Fixed value")
%             BMelem(index1,index2) = 0;
%         end
%         if strcmp(boundarytype,"Zero gradient")
%             BMelem(index1,index2) = 0;
%         end
%     end
% end
end