function A = triangulation2adjacency(faces,verts)

% triangulation2adjacency - compute the adjacency matrix
%   of a given triangulation, i.e. A(i,j)=1 if edge e_ij exists.
%   A = triangulation2adjacency(faces);
% or for getting a weighted "adjacency" matrix using length of each edge,
% i.e. A(i,j)=dist(i,j) if edge e_ij exists.
%   A = triangulation2adjacency(faces,vertex);
%
%   Adapted by jjcao 2012
%   Copyright (c) 2005 Gabriel Peyr?

if nargin < 2; verts=[]; end;

[verts,faces] = check_face_vertex(verts,faces);
faces = double(faces)';
verts = double(verts)';

if nargin < 2
%     % the matlab operation should be faster!
%     if exist('adjacency_matrix', 'file')
%         [I J] = adjacency_matrix(faces);        
%         A = sparse(I,J,ones(size(I)));
%         % make sure that all edges are symmetric
%         A = max(A,A');        
%     else        
        A = sparse([faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
                   [faces(:,2); faces(:,3); faces(:,1); faces(:,3); faces(:,1); faces(:,2)], ...
                   1.0);
        % avoid double links
        A = double(A>0);
%     end
else
    nvert = max(max(faces));
    nface = size(faces,1);
    A = spalloc(nvert,nvert,3*nface);

    for i=1:nface
        for k=1:3
            kk = mod(k,3)+1;
            v = verts(faces(i,k),:)-verts(faces(i,kk),:);
            A(faces(i,k),faces(i,kk)) = sqrt( sum(v.^2) );    % euclidean distance
        end
    end 
    
    % make sure that all edges are symmetric
    A = max(A,A');
end