function [A,fverts,E1,E2] = compute_dual_graph(faces,verts, with_edge)

% compute_dual_graph - compute the dual graph of a given triangulation
%
%   [A,fverts] = compute_dual_graph(faces,verts);
%
%   'A': adjacency matrix of the abstract dual graph (recall that this graph link togeter adjacent faces
%   in the triangulation). A(i,i)=0
%
%   'fverts' the position of the verts of the dual graph, i.e. the centroids of the input faces.
%   
%   if with_edge == true, then
%       'E1' if A(i,j)=1, then E1(i,j) is 1st vertex index of edge shared by faces i & j in input mesh.
%       'E2' if A(i,j)=1, then E2(i,j) is 2nd vertex index of edge shared by faces i & j in input mesh.
%
%   changed by jjcao, 2012
%
%   Copyright (c) 2004 Gabriel Peyr?
if nargin < 3
    with_edge = false;
end

[verts,faces] = check_face_vertex(verts,faces);
nface = size(faces,2);
fring = compute_face_ring(faces);
fverts = compute_face_centers(verts, faces, nface);

%%
sprintf('compute_dual_graph: ')
tic;
A = sparse(nface,nface); % from zeros => sparse, jjcao
for i=1:nface
    ring = fring{i};
    for j=1:length(ring)
        A(i,ring(j)) = 1;
    end
end
toc;

A = max(A,A');

%%
if with_edge
    sprintf('compute_dual_graph, with edge info: ')
    tic;
    E1 = sparse(nface,nface);
    E2 = E1;
    for i=1:nface
        ring = fring{i};
        for j=1:length(ring)
            vids = intersect(faces(:,i),faces(:,ring(j))); 
            E1(i,ring(j)) = min(vids);
            E2(i,ring(j)) = max(vids);
        end
    end
    E1 = max(E1,E1');
    E2 = max(E2,E2');
    toc;
else
    E1=[]; E2=[];
end

%%
function fverts = compute_face_centers(verts, faces, nface)
% compute the center of the faces 
if isempty(verts)
    fverts = [];
else
    fverts = [   ...
    sum(reshape(verts(1,faces),[3 nface]), 1)/3; ...
    sum(reshape(verts(2,faces),[3 nface]), 1)/3; ...
    sum(reshape(verts(3,faces),[3 nface]), 1)/3 ];     
end

