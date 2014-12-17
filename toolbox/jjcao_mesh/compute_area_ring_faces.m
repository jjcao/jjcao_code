function areas = compute_area_ring_faces(verts,faces, rings)
% compute_area_ring_faces - compute area for 1-ring faces of each vertex.
%
%   areas = compute_area_ring_faces(verts,faces, rings);
%
%   Copyright (c) 2009 JJCAO

if exist('rings','var') ==0
    rings = compute_vertex_face_ring(faces);
end

areas = zeros(size(rings,2),1);
n = size(areas,1);
for i = 1:n
    ring = rings{i};
    A = cross(verts(faces(ring,2),:)- verts(faces(ring,1),:), verts(faces(ring,3),:)- verts(faces(ring,1),:));
    tmpAreas = 0.5 * sqrt(A(:,1).^2+A(:,2).^2+A(:,3).^2);
    areas(i) = sum(tmpAreas);
end

% %% old way
% areas1 = zeros(size(rings,2),1);
% for i = 1:n
%     ring = rings{i};
%     for b = ring
%         bf = faces(b,:);
%         vi = verts(bf(1),:); vj = verts(bf(2),:); vk = verts(bf(3),:);
%         areas1(i) = areas1(i) + 0.5 * norm(cross(vi-vk,vi-vj));
%     end
% end
% sum(abs(areas-areas1))
