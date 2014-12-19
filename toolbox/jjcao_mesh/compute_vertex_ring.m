function vring = compute_vertex_ring(face)

% compute_vertex_ring - compute the 1 ring of each vertex in a triangulation.
%
%   vring = compute_vertex_ring(face);
%
%   vring{i} is the set of vertices that are adjacent
%   to vertex i.
%
%   Copyright (c) 2004 Gabriel Peyr?

[tmp,face] = check_face_vertex([],face);

nverts = max(max(face));

A = triangulation2adjacency(face);
[i,j,s] = find(sparse(A));

% create empty cell array
vring{nverts} = [];

for m = 1:length(i)
    vring{i(m)}(end+1) = j(m);
end
end

%% old method, without adjacency
% function vring = compute_vertex_ring(face)
% ring{nverts} = [];
% 
% for i=1:nfaces
%     for k=1:3
%         ring{face(k,i)} = [ring{face(k,i)} face(mod((k+1),3),i) face(mod((k+2),3),i)];
%     end
% end
% 
% for i=1:nverts
%     ring{i} = unique(ring{i});
% end