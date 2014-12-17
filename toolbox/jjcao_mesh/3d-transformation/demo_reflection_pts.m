% demo_reflection
%  
%
% Copyright (c) 2013 Junjie Cao
clc;clear all;close all;
addpath(genpath('../../'));
DEBUG = 1;
USE_PAIR = 0;

%% input
verts0 = read_pts('GlobalReflectionSymmDetector-cpts_good.obj');
verts = read_pts('GlobalReflectionSymmDetector-cpts.obj');

figure('Name','input'); set(gcf,'color','white');hold on;
scatter3(verts(:,1),verts(:,2), verts(:,3),20,'b','filled');
axis off; axis equal; mouse3d;


%% 
verts1 = verts;
verts2 = verts;    

center = mean(verts);
% center = mean(verts0);
normal = [1,0,0];
% center = [0,2.77556e-17,0.57881];
% normal = [-1,2.67949e-08,0];
drawVector3d(center, normal, 'g');
scatter3(center(:,1),center(:,2), center(:,3),100,'r','filled');

normal = normal./sqrt(dot(normal,normal));
plane0 = createPlane(center, normal);
h = drawPlane3d(plane0, 'g');
set(h, 'faceAlpha', 0.5);
%%
verts2new = reflect_points(verts2, center, normal);
scatter3(verts2new(:,1),verts2new(:,2), verts2new(:,3),20,'r','filled');

%% distance 
[mean_dist, max_dist] = dist_between_parts(verts1, verts2new)

% p1 = min(verts1);
% p2 = max(verts1);
% sqrt( sum( (p1-p2).^2))
% 
% p1 = min(verts2);
% p2 = max(verts2);
% sqrt( sum( (p1-p2).^2))