% test_compute_voronoi_triangulation_mesh
%
%
% Copyright (c) 2012 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])
addpath ([MYTOOLBOXROOT '/jjcao_common'])

[verts,faces] = read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
%%
nverts = size(verts,1);
nbr_landmarks = 100; % number of points, eg. 400
options.W = ones(nverts,1); % the speed function, for now constant speed to perform uniform remeshing
landmark = perform_farthest_point_sampling_mesh(verts,faces,[],nbr_landmarks,options);% perform the sampling of the surface
% perform the front propagation, Q contains an approximate segementation
[D,S,Q] = perform_fast_marching_mesh(verts, faces, landmark);

%%
[vertex_voronoi,faces_voronoi] = compute_voronoi_triangulation_mesh(Q,verts,faces);
figure;
options.voronoi_edges = [];
options.alfa = 0.2;
plot_fast_marching_mesh(verts,faces, Q(:,1), [], options); 
hold on;
options.alfa = 1;
plot_mesh(vertex_voronoi',faces_voronoi',options);
view3d zoom;