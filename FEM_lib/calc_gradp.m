function [gradpx,gradpy] = calc_gradp(myCFD)
% calc grad P
% cp = myCFD.Solution.p;
% 
% gradpx = zeros(size(myCFD.Mesh.Nodes,2),1);
% gradpy = zeros(size(myCFD.Mesh.Nodes,2),1);
% for el_index = 1:size(myCFD.Mesh.Elements,1)
%     topology = myCFD.Mesh.topology.element;
%     elmat = myCFD.Mesh.Elements;
%     x = myCFD.Mesh.Nodes(1,:);
%     y = myCFD.Mesh.Nodes(2,:);
% 
%     xc = zeros(1,topology);
%     yc = zeros(1,topology);
%     for index1 = 1:topology
%         xc(index1) = x(elmat(el_index,index1));
%         yc(index1) = y(elmat(el_index,index1));
%     end
%     B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);
%     
%     beta  = B_mat(2,1:3);
%     gamma = B_mat(3,1:3);
%     
%     for i=1:topology
%         global_index = elmat(el_index,i);
%         gradpx(global_index) = gradpx(global_index) + beta(i)*cp(global_index);
%         gradpy(global_index) = gradpy(global_index) + gamma(i)*cp(global_index);
%     end
% end
% %scale gradp accourding to connectivity
% connectivity = zeros(size(myCFD.Mesh.Nodes,2),1);
% for i=1:size(myCFD.Mesh.Nodes,2)
%     connectivity(i) = nnz(reshape(elmat,numel(elmat),1)==i);
% end
% gradpx = gradpx./connectivity;
% gradpy = gradpy./connectivity;

topology = myCFD.Mesh.topology.element;
elmat = myCFD.Mesh.Elements;
p = myCFD.Mesh.Nodes;
t = ones(size(myCFD.Mesh.Elements',1)+1,size(myCFD.Mesh.Elements',2));
t(1:3,:) = myCFD.Mesh.Elements';

gradpx = zeros(size(myCFD.Mesh.Nodes,2),1);
gradpy = zeros(size(myCFD.Mesh.Nodes,2),1);
[elgradpx,elgradpy] = pdegrad(p,t,myCFD.Solution.p);
for el_index = 1:size(myCFD.Mesh.Elements,1)
    for i=1:topology
         global_index = elmat(el_index,i);
         gradpx(global_index) = gradpx(global_index) + elgradpx(el_index);
         gradpy(global_index) = gradpy(global_index) + elgradpy(el_index);
     end
end
%scale gradp accourding to connectivity
connectivity = zeros(size(myCFD.Mesh.Nodes,2),1);
for i=1:size(myCFD.Mesh.Nodes,2)
    connectivity(i) = nnz(reshape(elmat,numel(elmat),1)==i);
end
gradpx = gradpx./connectivity;
gradpy = gradpy./connectivity;
end

