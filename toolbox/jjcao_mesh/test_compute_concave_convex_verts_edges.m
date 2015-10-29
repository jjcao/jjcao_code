% test_compute_concave_convex_verts_edges
%
% Copyright (c) 2012 JJCAO

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_plot'])

%% load a mesh
test_file = {[MYTOOLBOXROOT 'data/cube_f300.off'],[MYTOOLBOXROOT 'data/fandisk.off'],[MYTOOLBOXROOT 'data/100.off']};
filename = test_file{2};
[verts,faces] = read_mesh(filename);
%%
options.num_concavenhb = 1;
[concave_id,convex_id, concave_fun, concave_thres,convex_thres] = compute_concave_convex_vertices(verts,faces,options);
% options.concave_thres = concave_thres*5;
% [concave_id,convex_id, concave_fun, concave_thres,convex_thres] = compute_concave_convex_vertices(verts,faces,options);

%% display results
%concave verts
figure('Name','Concave vertices'); set(gcf,'color','white');hold on;
h = plot_mesh(verts,faces);
h=scatter3(verts(concave_id,1),verts(concave_id,2),verts(concave_id,3),80,'b','filled');
view3d zoom;
% convex verts
h=scatter3(verts(convex_id,1),verts(convex_id,2),verts(convex_id,3),80,'r','filled');

%% concave edge
bShare = 1;
[sidx, edges] = edges_sharing_vertex_property([], faces, concave_id, bShare);

figure('Name','Concave edge'); set(gcf,'color','white');hold on;
h = plot_mesh(verts,faces);
sizee = 2; color = [0,0,1];
h = plot_edges(edges(sidx,:), verts, sizee, color);
view3d zoom;

%% convex edge
[sidx, edges] = edges_sharing_vertex_property([], faces, convex_id, bShare);
sizee = 2; color = [1,0,0];
h = plot_edges(edges(sidx,:), verts, sizee, color);
% h = plot_edges(edges(sidx,:), verts, sizee, (1:length(sidx))');
