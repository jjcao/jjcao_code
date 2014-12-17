function test_compute_boundary(file)
%   
%   Copyright (c) 2013 jjcao

if nargin < 1
    clear;clc;close all;
    addpath(genpath('../../'));
%     file = 'catHead_v131.off';
    file = 'cat5.off';
end

[vertices,faces]=read_mesh(file);

%%
tic
boundary=compute_boundary(faces);
toc
boundary=cell2mat(boundary)';
figure;
plot_mesh(vertices, faces);
shading 'faceted';
% display starting points
hold on;
ms = 25;
for i=1:length(boundary)
    cv = vertices(boundary(i),:);
    h = plot3( cv(1),cv(2), cv(3), 'r.');
    set(h, 'MarkerSize', ms);    
end
hold off;

%%
tic
[Mhe,Mifs] = to_halfedge(vertices, faces);
boundary = Mhe.boundary_vertices;
toc
if iscell(boundary)
    boundary=cell2mat(boundary)';
end
figure;
plot_mesh(vertices, faces);
shading 'faceted';
% display starting points
hold on;
ms = 25;
for i=1:length(boundary)
    cv = vertices(boundary(i),:);
    h = plot3( cv(1),cv(2), cv(3), 'r.');
    set(h, 'MarkerSize', ms);    
end
hold off;