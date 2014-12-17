% function test_compute_dual_graph
%
% 
% Copyright (c) 2012 Junjie Cao
%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_plot'])

%% load a mesh
test_file = {[MYTOOLBOXROOT 'data/wolf0.off'], [MYTOOLBOXROOT 'data/catHead_v131.off']};
M.filename = test_file{2};
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);
M.edges = compute_edges(M.faces); 

%% display mesh
figure('Name','Mesh'); set(gcf,'color','white');hold on;
% options.face_vertex_color = M.verts(:,3);
options.face_vertex_color = M.verts;
h = plot_mesh(M.verts, M.faces, options);
camorbit(0,0,'camera'); axis vis3d; view(-90, 0);view3d rot;colorbar('off');

%% display dual mesh
[A,vertex1] = compute_dual_graph(M.faces,M.verts);
figure('Name','Dual Mesh'); set(gcf,'color','white');hold on;
options.ps=10;
h = plot_graph(A,vertex1, options);
axis off;    axis equal; 
camorbit(0,0,'camera'); axis vis3d; view(-90, 0);view3d rot;
