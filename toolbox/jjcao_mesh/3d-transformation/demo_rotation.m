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
if DEBUG
for i=1:length(parts)
    h=trisurf(parts(i).faces,vertices(:,1),vertices(:,2),vertices(:,3), 'FaceVertexCData', vertices(:,3));     
end
set(h, 'edgecolor', 'none');colormap jet(256); axis off; axis equal; mouse3d;
end
%% display pair
pair = [6,8];
[parts( pair(1)).name, '-', parts( pair(2)).name]

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

%% rotate verts2 around centerAxis 
line = [-0.0015451344444444477,-0.15137158888888888,0.34855222222222221,...
    0.0037052556183702194,0.0037047106307757718,-0.99998627300575715];
angle = 2.0943950666666669;
% angle = -0.5*pi;
rot = createRotation3dLineAngle(line, angle);
[axisR angle2] = rotation3dAxisAndAngle(rot);
angle2
verts2new = transformPoint3d(verts2,rot);
scatter3(verts2new(:,1),verts2new(:,2), verts2new(:,3),100,'r','filled');

[dist1, dist2] = dist_between_parts(verts1, verts2);
(dist1+dist2)*0.5
[dist1, dist2] = dist_between_parts(verts1, verts2new);
(dist1+dist2)*0.5