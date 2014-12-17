% test_plot_mesh
% 
% If texture is not needed, use trisurf, see eg_trisurf.m
% the texture display is not right, still!!, no time to fix it.
%
% Copyright (c) 2013 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])


%% load a mesh
% test_file = {[MYTOOLBOXROOT '/data/Ball_dABF.obj'],[MYTOOLBOXROOT '/data/Cow_dABF.obj'],[MYTOOLBOXROOT '/data/Gargoyle_ABF.obj']};
% M.filename = test_file{1};
% [M.verts,M.faces, normal, M.uv] = read_mesh(M.filename);

M.filename = [MYTOOLBOXROOT '/data/catHead_v131.off'];
[M.verts,M.faces] = read_mesh(M.filename);
load([MYTOOLBOXROOT '/data/catHead_v131-combinatorial-circle.mat']);
M.uv = uv';

% if ~isempty(M.uv)
%     M.uv(:,1) = M.uv(:,1) - min(M.uv(:,1));
%     M.uv(:,2) = M.uv(:,1) - min(M.uv(:,2));
%     M.uv(:,1) = M.uv(:,1)/max(M.uv(:,1));
%     M.uv(:,2) = M.uv(:,2)/max(M.uv(:,2));
% end
M.nverts = size(M.verts,1);
M.edges = compute_edges(M.faces); 

%% display point cloud
figure;set(gcf,'color','white');hold on;
% scatter3(M.verts(:,1),M.verts(:,2), M.verts(:,3),10,'b','filled');
scatter3(M.verts(:,1),M.verts(:,2), M.verts(:,3),10,M.verts(:,2),'filled');
axis off;    axis equal;   set(gcf,'Renderer','OpenGL');
camorbit(0,0,'camera'); axis vis3d; view(-90, 0); mouse3d

%% display mesh
figure('Name','Mesh'); set(gcf,'color','white');hold on;
% options.face_vertex_color = M.verts(:,3);
options.face_vertex_color = M.verts;
options.alfa = 0.6;
h = plot_mesh(M.verts, M.faces, options);camorbit(0,0,'camera'); axis vis3d; view(-90, 0);mouse3d
set(h, 'edgecolor', 'none'); % cancel display of edge.
%colormap jet(256);
colorbar('off');

%% display bundary vertices
boundary=compute_boundary(M.faces);
scatter3(M.verts(boundary{:},1),M.verts(boundary{:},2), M.verts(boundary{:},3),50,'r','filled');

%% display texture, why not correct?
figure('Name','Texture'); set(gcf,'color','white');hold on;
options=[];
% options.texture = rgb2gray(imread([MYTOOLBOXROOT '/data/Texture.bmp']));
options.texture = imread([MYTOOLBOXROOT '/data/Texture.bmp']);
options.texture_coords = M.uv;
h = plot_mesh(M.verts, M.faces, options);camorbit(0,0,'camera'); axis vis3d; view(-90, 0);mouse3d
