% test_compute_angle_defect
%
% Copyright (c) 2012 JJCAO
%% initialize & read mesh
clear;clc;close all;
addpath(genpath('../../'));

tau = 1.2;% options for display
test_file = {'data/cube_602.off','data/armadillo_v502.off', 'data/armadillo_v4326.off','data/657_7k.off', 'data/wolf0.off','data/dancer_v5k.off'};
filename = test_file{3};
[verts,faces] = read_mesh(filename);

options.bsaturate = 1;
%% curvature by angle defect & corresponding conformal factor
Cgauss = compute_angle_defect(verts, faces);
options.figname='Cgauss by angle defect';options.position='northwest';
plot_mesh_scalar(verts, faces, Cgauss, options);

options.Cgauss = Cgauss;
cf = compute_conformal_factor(verts, faces, options);  
options.figname='Conformal factor by angle defect';options.position='north';
plot_mesh_scalar(verts, faces, cf, options);

%% curvature by normal cycle & corresponding conformal factor
options.curvature_smoothing = 3; % defualt is 3
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss1,Normal] = compute_curvature(verts,faces,options);
options.figname='Cgauss by normal cycle';options.position='northeast';
plot_mesh_scalar(verts, faces, Cgauss1, options);

options.Cgauss = Cgauss1;
cf1 = compute_conformal_factor(verts, faces, options);  
options.figname='Conformal factor by normal cycle';options.position='southwest';
plot_mesh_scalar(verts, faces, cf1, options);

