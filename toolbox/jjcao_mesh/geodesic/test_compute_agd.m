%% test_compute_agd
%
% Copyright (c) 2013 Shuhua Li, Junjie Cao
clear;clc;close all;
addpath(genpath('../../'));

%% compute agd
[verts,faces] = read_mesh('wolf0.off');
nbr_landmarks = 100; % the number of landmarks 
DEBUG=0;     % 1,view the voronoi cells and agd cells; 0,otherwise;
normalize=0; %1£¬normalize agd to 0-1£»0,otherwise
[agd,landmark]=compute_agd(verts,faces,nbr_landmarks,DEBUG,normalize);
%% display agd
figure('Name','agd'); set(gcf,'color','white'); 
options.face_vertex_color =agd;
h = plot_mesh(verts, faces, options);
colormap(jet(nbr_landmarks)); 
view3d rot; lighting none;