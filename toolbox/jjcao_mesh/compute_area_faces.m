function areas = compute_area_faces(verts,faces)
% Compute the area of each triangle
%
% 
%
% Copyright (c) 2013 Baochang Han, Junjie Cao
%

A = cross(verts(faces(:,2),:)- verts(faces(:,1),:), verts(faces(:,3),:)- verts(faces(:,1),:));
areas = 0.5 * sqrt(A(:,1).^2+A(:,2).^2+A(:,3).^2);

% %% old
% function [area areas] = compute_area_faces(vertices,faces)
% %for trianglar mesh only
% areas = zeros(size(faces,1),1);
% n = size(areas,1);
% for i = 1:n
%     vi = vertices(faces(i,1),:); vj = vertices(faces(i,2),:); vk =vertices(faces(i,3),:);
%     areas(i) = 0.5 * norm(cross(vi-vk,vi-vj));
% end 
% area = sum(areas);