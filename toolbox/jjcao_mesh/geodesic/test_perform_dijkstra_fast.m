%% test_perform_dijkstra_fast
%
% Copyright (c) 2014 Junjie Cao
clear;clc;close all;
addpath(genpath('../../'));
MYTOOLBOXROOT='../..';

[verts,faces] = read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
% [verts,faces] = read_mesh('/data/catHead_v131.off');
nverts = size(verts,1);
diagLength = 100;
[verts,diagLength] = normalize_vertex3d(verts,diagLength);

%% compute adjacency matrix and similarity matrix which's element is distance of each corresponding edge
A = triangulation2adjacency(faces);
ind = find(A>0);
[I, J] = ind2sub(size(A), ind);
dist2 = sum((verts(I,:) - verts(J,:)).^2, 2);
W = sparse(I,J, sqrt(dist2) );

start_points = [1,2];
tic
D = perform_dijkstra_fast(W, start_points)';
toc
D1 = perform_dijkstra_fast(A, start_points)';
sum(sum(abs(D - D1)))

%% show the result
i = 1;
figure('Name','shortest distance by Dijkstra using distance matrix'); set(gcf,'color','white');
h=trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
    'FaceVertexCData', D(:,i), 'edgecolor','none');
axis off;    axis equal;   set(gcf,'Renderer','OpenGL'); shading interp;
camorbit(0,0,'camera'); axis vis3d; view(-90, 0); mouse3d; colorbar;hold on;

scatter3(verts(start_points(i),1),verts(start_points(i),2), verts(start_points(i),3), ...
    50,[1,0,0],'filled');

figure('Name','shortest distance by Dijkstra using adjacency (all edges'' distance == 1) matrix'); set(gcf,'color','white');
h=trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
    'FaceVertexCData', D1(:,i), 'edgecolor','none');
axis off;    axis equal;   set(gcf,'Renderer','OpenGL'); shading interp;
camorbit(0,0,'camera'); axis vis3d; view(-90, 0); mouse3d; colorbar;hold on;
scatter3(verts(start_points(i),1),verts(start_points(i),2), verts(start_points(i),3), ...
    50,[1,0,0],'filled');
