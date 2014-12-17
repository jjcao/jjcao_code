
%   Copyright (c) 2012 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_common'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
addpath ([MYTOOLBOXROOT '/accelerated_mds'])
DEBUG=1;
%% input
filename = 'data/Isis_dABF.obj';
[verts, faces] = read_obj(filename);
nface = size(faces,1);
if DEBUG
    figure('Name','Input model'); set(gcf,'color','white'); 
    h = plot_mesh(verts, faces);
    view3d rot; lighting none;
end
D = compute_geodesic_distance_matrix(verts, faces, 1:length(verts));

%%
ndim = 3;
% options.landmarks = 1:length(verts);
options.landmarks = round(rand(100,1)*length(verts));
options.method = 'classical_mds';
% options.method = 'rre';
verts1 = isomap_for_mesh(D,ndim,options);
GDD = compute_geodesic_distortion(verts1, D);
sprintf('mean geodesic distortion: %s', mean(GDD))
sprintf('max geodesic distortion: %s', max(GDD))
%%
if DEBUG
    figure('Name','parameterization'); set(gcf,'color','white'); 
%     verts1(:,3) = 0;
    h = plot_mesh(verts1, faces);
    view3d rot; lighting none;
end

%%
% use isomap.m of toolbox_graph
% A = triangulation2adjacency(faces,verts);
% A(A==0) = inf;
% verts1 = isomap(A,ndim,options);