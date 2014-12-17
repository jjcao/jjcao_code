%
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../../'));

%% load a mesh
test_file = {'/data/fandisk.off','/data/wolf0.off','/data/catHead_v131.off'};
M.filename = test_file{2};
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);

%% compute eigen value and function
tic
withAreaNormalization = 1;
adjustL = 0;
[eigvec, eigv]=compute_Laplace_eigen(M.verts,M.faces,100, withAreaNormalization, adjustL);
toc

%% display some eigen functions
shading_type = 'interp';
for i = 2:ceil(length(eigv)/8):length(eigv)
% for i = 2:5
    figure('Name', sprintf('eigen function: %d', i));    
    h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
        'FaceVertexCData', eigvec(:,i), 'FaceColor',shading_type); 
    axis off; axis equal; set(h, 'edgecolor', 'none'); mouse3d;
    options = []; options.edgecolor = 'black';
%     options.Nlines = 10;
    options.values = [0];
    plot_contours(M.faces,M.verts,eigvec(:,i), options);
    title(sprintf('eigen function: %d', i));
end



