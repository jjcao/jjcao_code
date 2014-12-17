% test_compute_geodesic_mesh
%
%
% Copyright (c) 2012 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_common'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])

%%
[verts,faces] = read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
nverts = size(verts,1);
%%
sids = [2686]+1;% start points
eids = [1000, 2000];% end points
% perform the front propagation, Q contains an approximate segementation
[D,S,Q] = perform_fast_marching_mesh(verts, faces, sids);

options.method='continuous';
[path,vlist,plist] = compute_geodesic_mesh(D, verts, faces, eids, options);
% options.start_points = landmark;
plot_fast_marching_mesh(verts,faces, D, path, options);hold on;
view3d zoom;

%%
options.method='discrete';
[path,vlist,plist] = compute_geodesic_mesh(D, verts, faces, eids, options);
figure;
% options.start_points = landmark;
plot_fast_marching_mesh(verts,faces, D, path, options);hold on;
view3d zoom;


