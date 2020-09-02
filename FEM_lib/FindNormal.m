function [nx,ny] = FindNormal(bndelement,mesh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%bndelement = [1,9]; % testcase i know exists, element 475
% find the element index that also contains the boundary
checkindex1 = any(mesh.Elements(:,:)==bndelement(1),2);
checkindex2 = any(mesh.Elements(:,:)==bndelement(2),2);
index = find(and(checkindex1,checkindex2));
if length(index)>1
    warning('FindNormal: boundary element shares edge with more than one element in the mesh')
end

% calculate normal vector normal to boundary
vec1 = zeros(1,3);
vec2 = zeros(1,3);
zaxis = [0,0,1];

vec1(1:2) = mesh.Nodes(:,bndelement(1));
vec2(1:2) = mesh.Nodes(:,bndelement(2));

normal = cross(vec2-vec1,zaxis);

% check if other vertecies is ahead or behind the chosen normal, flip if
% needed
element = mesh.Elements(index,:);
ip = zeros(1,length(element));
for i=1:length(element)
    if any(bndelement==element(i))
        ip(i) = 0;
        continue % if element is checked agains one on the boundray skip it
    end
    testvec = mesh.Nodes(:,element(i));
    ip(i) = dot(normal(1:2),testvec);
end

if all(ip>=0)
    %normal is pointing in the wrong direction
    normal = normal*(-1);
elseif all(ip<=0)
    %no need to flip the normal
else
    error('FindNormal, this should not be possible')
end
nx = normal(1);
ny = normal(2);
end

