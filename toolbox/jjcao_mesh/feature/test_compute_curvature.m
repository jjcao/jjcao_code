% test_curvature_confomal_factor
%
%     curvatures
%     conformal factor
%
% Copyright (c) 2012 JJCAO
%% initialize & read mesh
clear;clc;close all;
addpath(genpath('../../'));

tau = 1.2;% options for display
test_file = {[MYTOOLBOXROOT 'data/cube_602.off'],[MYTOOLBOXROOT 'data/fandisk.off'],[MYTOOLBOXROOT 'data/100.off'],'E:\jjcao_data\MeshsegBenchmark-1.0\data\off\99.off'};
filename = test_file{4};
[verts,faces] = read_mesh(filename);
%% curvature by normal cycle
options.curvature_smoothing = 3; % defualt is 3
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(verts,faces,options);

options.figname='Cmax';options.position='northwest';
plot_mesh_scalar(verts, faces, Cmax, options);

options.figname='Cmin';options.position='north';
plot_mesh_scalar(verts, faces, Cmin, options);

options.figname='Cmean';options.position='northeast';
plot_mesh_scalar(verts, faces, Cmean, options);

options.figname='Cgauss';options.position='southwest';
plot_mesh_scalar(verts, faces, Cgauss, options);

options.figname='abs(Cmin)+abs(Cmax)';options.position='south';
plot_mesh_scalar(verts, faces, abs(Cmin)+abs(Cmax), options);
