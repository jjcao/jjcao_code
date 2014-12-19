% test_vertex_area
%
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../'));

[verts,faces] = read_mesh('catHead_v131.off');
va = vertex_area(verts,faces,0);

figure('Name','2'); set(gcf,'color','white');hold off;
h=trisurf(faces, verts(:,1), verts(:,2), verts(:,3), 'FaceVertexCData', va);
axis off; axis equal; mouse3d;
colormap jet;
colorbar