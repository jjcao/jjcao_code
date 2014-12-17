% test_perform_farthest_point_sampling_mesh
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

nbr_landmarks = 100; % num of samples
DEBUG_ = 1;

[verts, faces]= read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
landmark = perform_farthest_point_sampling_mesh(verts,faces,[],nbr_landmarks);
    
if DEBUG_ 
    figure; set(gcf,'color','white');hold on;axis off; axis equal;
    h = plot_mesh(verts, faces);
    scatter3(verts(landmark,1), verts(landmark,2),verts(landmark,3),40,'b','filled');
    view3d rot;
end

   
