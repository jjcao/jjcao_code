function [sidx, edges] = edges_sharing_vertex_property(edges, faces, vidx, bShare)
% vidx: index of vertices owing special property.
% 
%
%   Copyright (c) 2012 Hui Wang, JJCAO

if nargin < 4
    bShare = true;
end

if isempty(edges)
    A = triangulation2adjacency(faces);
    [row, col] = find(A);
    tmp=row<=col;
    row=row(tmp);
    col=col(tmp);
    edges = [row, col];
else
    row = edges(:,1);
    col = edges(:,2);
end

%%
SIGN = zeros(max(max(faces)), 1);
SIGN(vidx) = 1;
fun1 = SIGN(row);
fun2 = SIGN(col);

sidx = zeros(size(fun1));
if bShare
    sidx = (fun1 + fun2) >1;
else
    sidx = (fun1 + fun2) >0;
end

return;

%% function [concave,edges] = compute_concave_edges(verts,faces)
% % old
% % warning, is not good!
% %
% % Copyright (c) 2012 Hui Wang, JJCAO
% A = compute_dual_graph(faces,verts);
% [normal,fnormal] = compute_normal(verts,faces);
% 
% [row,col] = find(A>0);
% tmp=row<=col;
% row=row(tmp);
% col=col(tmp);
% edges = [row, col];
% concave = zeros(size(row));
% tic;
% for i=1:length(row)
%     vids = intersect(faces(row(i),:),faces(col(i),:));   
%     edge = verts(vids(1),:)-verts(vids(2),:);
%     edge = edge/sqrt(sum(edge.^2));
%     edges(i,:) = vids;
%     concave(i) = dot(cross(fnormal(row(i),:), fnormal(col(i),:)), edge);
% end
% toc;
% % ConcaveWeight = 0.2 + 0 * edgeidx;
% % ConcaveWeight(find(edgeidx < 0)) = 1.0;