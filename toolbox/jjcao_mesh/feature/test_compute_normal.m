% function test_compute_normal
%
% 
% Copyright (c) 2012 Junjie Cao
%%
clear;clc;close all;
addpath(genpath('../../'));

[vertices,faces] = read_mesh('data/cube_f300.off');
[normals, fnormals] = compute_normal(vertices, faces);

%% display 
figure;set(gcf,'color','white');hold on;
plot_mesh(vertices, faces);

options.subsample_normal = 0.8;
options.normal_scaling = 0.2;
h = plot_face_normal(vertices, faces, fnormals', options);
