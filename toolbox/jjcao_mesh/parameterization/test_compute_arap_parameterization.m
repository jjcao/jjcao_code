
%   Copyright (c) 2012 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_common'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/parameterization/arap'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
DEBUG=1;
%% input
filename = 'data/Isis_dABF.obj';
[verts, faces, normal, vt] = read_obj(filename);
if DEBUG
    figure('Name','Input model'); set(gcf,'color','white'); 
    h = plot_mesh(verts, faces);
    view3d rot; lighting none;
end
D1 = compute_geodesic_distance_matrix(verts, faces, 1:length(verts));

%%
[vt iterations] = compute_arap_parameterization(verts,faces,vt, 0.001);
vt(:,3) = 0;
D2 = compute_geodesic_distance_matrix(vt, faces, 1:length(verts));
N = length(verts);
GDD = sqrt(sum((D1-D2).^2,2)./(N-1));
sprintf('mean geodesic distortion: %s', mean(GDD))
sprintf('max geodesic distortion: %s', max(GDD))

%% %%%%%% Show para result %%%%%%%%%%%%%%%%%%%%%
figure('Name','Input model'); set(gcf,'color','white'); 
h = plot_mesh(vt, faces);
view3d rot; lighting none;
title([num2str(iterations),' Iterations']);
