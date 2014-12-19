%test_compute_HKS
%
%
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../'));
%% Load mesh
M.filename = 'wolf0.off';
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);

%% Laplacian eigen
nbasis = 100;
withAreaNormalization = true;
adjustL = false;
[M.eigvector,M.eigvalue, A]=compute_Laplace_eigen(M.verts,M.faces,nbasis, withAreaNormalization, adjustL);

%% HKS
tic
scale = true;
M.hks = HKS(M.eigvector, M.eigvalue, A, scale);
toc

%% display HKS of some vertices
shading_type = 'interp';
for i = 1:ceil(size(M.hks,2)/6):size(M.hks,2)
% for i = 1:5
    figure('Name', sprintf('eigen function: %d', i));
    h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
        'FaceVertexCData', M.hks(:,i), 'FaceColor',shading_type); 
    axis off; axis equal; set(h, 'edgecolor', 'none'); mouse3d;
end