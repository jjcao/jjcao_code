%
%  
%
% Copyright (c) 2013 Junjie Cao
clc;clear all;close all;
addpath(genpath('../../'));
DEBUG = 1;

%% input
[parts, vertices]= read_parts_obj('CR1b.obj');
figure('Name','input'); set(gcf,'color','white');hold on;
% for i=1:length(parts)
%     h=trisurf(parts(i).faces,vertices(:,1),vertices(:,2),vertices(:,3), 'FaceVertexCData', vertices(:,3));     
% end
% set(h, 'edgecolor', 'none');colormap jet(256); axis off; axis equal; mouse3d;

%% display pair
pair = [2, 3];
faces1 = parts(pair(1)).faces;
faces2 = parts(pair(2)).faces;
verts1 = vertices(min(faces1):max(faces1),:);
verts2 = vertices(min(faces2):max(faces2),:);
figure('Name','parts'); set(gcf,'color','white');hold on;axis off; axis equal; mouse3d;
h=trisurf(parts(pair(1)).faces,vertices(:,1),vertices(:,2),vertices(:,3),'faceAlpha', 0.5);
set(h, 'FaceColor', 'g');
h=trisurf(parts(pair(2)).faces,vertices(:,1),vertices(:,2),vertices(:,3),'faceAlpha', 0.5);
set(h, 'FaceColor', 'b');
scatter3(verts1(:,1),verts1(:,2), verts1(:,3),50,'g','filled');
scatter3(verts2(:,1),verts2(:,2), verts2(:,3),50,'b','filled');

[center1, coeff1, axis1] = compute_pca_info(verts1);
[center2, coeff2, axis2] = compute_pca_info(verts2);
plot_axis(axis1);
plot_axis(axis2);

%% compute center
trans = (center1 - center2);
newvertices = vertices + repmat(trans, length(vertices), 1);
h=trisurf(parts(pair(2)).faces, newvertices(:,1),newvertices(:,2),newvertices(:,3));
set(h, 'FaceColor', 'b');
newverts2 = newvertices(min(faces2):max(faces2),:);
distance_between(verts1, newverts2)

center = center2 + trans*0.5;
centerAxis = axis2 + repmat(trans*0.5, size(axis2,1),1);
plot_axis(centerAxis);

%% rotate verts2 around centerAxis 
line = [center coeff2(:,1)'];
angle = pi*1;
rot = createRotation3dLineAngle(line, angle);
[axisR angle2] = rotation3dAxisAndAngle(rot);
angle2
rot1 = create_rotation3d_line_angle(center, coeff2(:,1)', angle);
[axisR angle1] = rotation3dAxisAndAngle(rot1);
angle1

newvertices = transformPoint3d(vertices,rot);
newvertices1 = transform_point3d(vertices,rot);

newverts2 = newvertices(min(faces2):max(faces2),:);
distance_between(verts1, newverts2)

figure('Name','After rotation'); set(gcf,'color','white');hold on;axis off; axis equal; mouse3d;
h=trisurf(parts(pair(1)).faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h, 'FaceColor', 'g');
h=trisurf(parts(pair(2)).faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h, 'FaceColor', 'b');
h=trisurf(parts(pair(2)).faces, newvertices(:,1),newvertices(:,2),newvertices(:,3));
set(h, 'FaceColor', 'r');