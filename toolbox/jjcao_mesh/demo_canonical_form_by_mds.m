%
% (C) Copyright JJCAO, 2012
% All rights reserved.

%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh/geodesic'])
addpath ([MYTOOLBOXROOT 'accelerated_mds'])

DEBUG=1;
KEEP_DETAIL = 1;

%% input
[verts,faces] = read_mesh([MYTOOLBOXROOT '/data/wolf0.off']);
nverts = size(verts,1);
if KEEP_DETAIL
    L = compute_mesh_laplacian(verts,faces,'conformal');%combinatorial;conformal;%spring
    delta_coords = L*verts; % keep details using Laplacian coordinates
else
    delta_coords = zeros(nverts,3); % remove details
end

%% landmark
nbr_landmarks = min(100,nverts); % number of subsamples

options.W = ones(nverts,1); % the speed function, for now constant speed to perform uniform remeshing
landmark = perform_farthest_point_sampling_mesh(verts,faces,[],nbr_landmarks,options);% perform the sampling of the surface
if DEBUG
    figure; set(gcf,'color','white');
    options.start_points = landmark;
    plot_fast_marching_mesh(verts,faces, [], [], options);hold on;view3d zoom;
end

options.X0 = [verts(landmark,1), verts(landmark,2), verts(landmark,3)];
D = compute_geodesic_distance_matrix(verts, faces, landmark);
options.method = 'rre';
options.cycles = 2;
options.iter = 10;
[X_] = mds(D,options);

%%
options.method = 'bi-harmonic';% 'hard', 'soft', 'bi-harmonic'
options.type = 'conformal'; % the other choice is 'combinatorial', 'distance', spring, 'conformal', or, 'authalic', mvc, 
options.normalize = 0;
options.symmetrize = 1; 
options.constraint_pos = X_;
constraint_weight = 10;
s_weights= ones(nverts,1)*1;
c_weights= ones(length(landmark),1)*constraint_weight;

[newverts, A] = compute_least_square_mesh(verts,faces,landmark,delta_coords,s_weights,c_weights,options);
figure;set(gcf,'color','white');
trisurf(faces,newverts(:,1),newverts(:,2),newverts(:,3)); axis off; axis equal; view3d rot;
