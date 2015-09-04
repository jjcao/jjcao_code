% supervertex by NCut
%
% a cluster may contain multi-patches, why?
% sometimes a cluster contains no face, why?
% 每次的结果可能不一样
%
% Copyright (c) 2012 Junjie Cao

clear;clc;close all;
addpath(genpath('../../'));

DEBUG = 1;
%% input
M.npatch = 200;
USE_CONCAVE_WEIGHT=0;
M.filename = '/data/fandisk.off'; %cube_f300, cube_f1200,fandisk,wolf0
[M.verts,M.faces] = read_mesh(M.filename);
% [M.verts,M.faces] = read_mesh();
[normal M.fnormal] = compute_normal(M.verts,M.faces);
M.fnormal = M.fnormal';
%%
[M.face_patch, M.npatch] = compute_supervertex_ncut(M.verts,M.faces,M.fnormal,M.npatch,USE_CONCAVE_WEIGHT);

%% show result
figure('Name','Supervertex by NCut'); movegui('southwest'); set(gcf,'color','white');
options.face_vertex_color = M.face_patch;
h = plot_mesh(M.verts, M.faces, options);view3d rot;
colormap(jet(M.npatch)); lighting none;
%%
figure('Name','Supervertex by NCut'); movegui('southeast'); set(gcf,'color','white');
options.face_vertex_color = M.face_patch;
h = plot_mesh(M.verts, M.faces, options);view3d rot;
colormap(colorcube(M.npatch)); lighting none;

%% save result
[pathstr, name, ext] = fileparts(M.filename);
save(sprintf('%s_p%d.mat', name, M.npatch), 'M');