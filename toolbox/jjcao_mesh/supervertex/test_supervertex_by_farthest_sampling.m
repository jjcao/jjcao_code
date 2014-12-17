% test_supervertex_by_farthest_sampling
%
%
% Copyright (c) 2013 Junjie Cao

clear;clc;close all;
addpath(genpath('../../'));


%%
[M.verts,M.faces] = read_mesh('/data/wolf0.off');
nverts = size(M.verts,1);
nbr_landmarks = 20; % number of points, eg. 400
[M.face_patch, patch_area,landmark] = supervertex_by_farthest_sampling(M.verts,M.faces,nbr_landmarks);
M.npatch = max(M.face_patch);

%% show result
figure('Name','Supervertex by farthest sampling'); movegui('southwest'); set(gcf,'color','white');
options.face_vertex_color = M.face_patch;
h = plot_mesh(M.verts, M.faces, options);view3d rot;
colormap(jet(M.npatch)); lighting none;