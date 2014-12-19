% function [] = eg_skeleton_Laplacian(filename)
% extract curve skeleton from a  triangular mesh
%
% by: JJCAO @ 2013
%
%% setting
clc;clear all;close all;
addpath(genpath('../../'));

M.filename = 'armadillo_v4326.off';% which file we should run on

%%
tic
disp(sprintf('1. read mesh, compute 1-ring:'));
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);
M.faceArea = compute_area_faces(M.verts, M.faces);
M.verts= M.verts/ sqrt(sum(M.faceArea)); %πÈ“ªªØ
toc

%% Step 1: Contract point cloud by Laplacian
tic
disp(sprintf('Contraction:'));
[M.cverts] = contraction_by_mesh_laplacian(M.verts, M.faces);
toc
