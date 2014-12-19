% test_compute_concave_edges
%
% Copyright (c) 2012 JJCAO

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_plot'])

%% load a mesh
test_file = {[MYTOOLBOXROOT 'data/cube_602.off'],[MYTOOLBOXROOT 'data/corner_tris_with_hole.off'],[MYTOOLBOXROOT 'data/100.off']};
filename = test_file{3};
[verts,faces] = read_mesh(filename);
%%
[concave,edges] = compute_concave_edges(verts,faces);
cedgeidx = find(concave<-3e-1);
%% concave edge
figure('Name','Concave edge'); set(gcf,'color','white');hold on;
h = plot_mesh(verts,faces);
sizee = 2; color = [0,0,1];
h = plot_edges(edges(cedgeidx,:), verts, sizee, color);
view3d zoom;
