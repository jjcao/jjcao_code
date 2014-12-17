% test_perform_fast_marching_mesh
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

%% read mesh
[verts,faces] = read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
nverts = size(verts,1);

%% compute geodesic distance of all vertices relative to the landmark vertices
landmark = [2686,4132]+1;% start points

% perform the front propagation, Q contains an approximate segementation
options.nb_iter_max = Inf; % trying with a varying number of iterations: 100, 1000, Inf
disp('running time for computing all geodesic distances from two landmark:');
tic
[D,S,Q] = perform_fast_marching_mesh(verts, faces, landmark, options);
toc

% display the result using the distance function
options.start_points = landmark;
figure('name', 'geodesic distance of all vertices relative to the landmark vertices');
plot_fast_marching_mesh(verts,faces, D, [], options);hold on;
% plot_fast_marching_mesh(verts,faces, Q, [], options);
scatter3(verts(landmark,1), verts(landmark,2),verts(landmark,3),40,'r','filled');
view3d zoom;



%% compute geodesic distance between a pair of vertices
start_id = 2686+1;% start points
end_id = 4132+1;
% perform the front propagation, Q contains an approximate segementation
options.end_points = end_id;

disp('running time for computing geodesic distance between a pair of vertices:');
tic
[D1,S,Q] = perform_fast_marching_mesh(verts, faces, landmark, options);
toc
% D1(end_id) is the distance from start_id to end_id;
abs( D(end_id) - D1(end_id))
