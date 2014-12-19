function A = compute_triangulation_angles(vertex,face)

% compute_triangulation_angles - compute the matrix of triangle angles
%
%   A = compute_triangulation_angles(vertex,face);
%
%   A(i,1:3) is the angles of face i.
%
%   example:
% [verts,faces] = read_mesh('/data/nefertiti_v299.off');
% A = compute_triangulation_angles(verts,faces);
% figure('name','hist of all angles'); set(gcf,'color','white');hist(A(:));
%
%
%   Copyright (c) 2008 Gabriel Peyre
%   changed 2012, Junjie Cao

[vertex,face] = check_face_vertex(vertex,face);

m = size(face,2);
A = zeros(3,m);

for i=1:3
    j1 = mod(i,3)+1;
    j2 = mod(i-2,3)+1;
    v1 = vertex(:,face(j1,:)) - vertex(:,face(i,:));
    v1 = v1 ./ repmat( sqrt(sum(v1.^2)), [3 1] );
    v2 = vertex(:,face(j2,:)) - vertex(:,face(i,:));
    v2 = v2 ./ repmat( sqrt(sum(v2.^2)), [3 1] );
    A(i,:) = acos( sum( v1.*v2, 1 ) );
end
%A = min(A); % min angle of each triangle
A = A'; % all angle of each triangle