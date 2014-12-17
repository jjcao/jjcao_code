% test_compute_least_square_mesh
%
% 
% Copyright (c) 2012 Junjie Cao
%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_plot'])

%% options
DEBUG=1;
constraint_weight = 10;
KEEP_DETAIL = 0;
USE_FAREST_POINT_SAMPLING = 0;
NUM_SAMPLES = 100;

%% input
[verts,faces] = read_mesh([MYTOOLBOXROOT 'data/wolf0.off']);
% [verts,faces] = create_mesh(50); verts = verts'; faces = faces';
 
nverts = size(verts,1);
if KEEP_DETAIL
    options.use_c_implementation = 1;
    tic
    L = compute_mesh_laplacian(verts,faces,'conformal',options);%combinatorial;conformal;%spring
    toc
    delta_coords = L*verts; % keep details using Laplacian coordinates
else
    delta_coords = zeros(nverts,3); % remove details
end

%% constraints
if USE_FAREST_POINT_SAMPLING
    nbr_landmarks = min(NUM_SAMPLES,nverts); % number of subsamples
    options.W = ones(nverts,1); % the speed function, for now constant speed to perform uniform remeshing
    landmark = perform_farthest_point_sampling_mesh(verts,faces,[],nbr_landmarks,options);% perform the sampling of the surface
    if DEBUG
        figure; set(gcf,'color','white');
        options.start_points = landmark;
        plot_fast_marching_mesh(verts,faces, [], [], options);hold on;view3d zoom;
    end
else
    landmark = 1:100:nverts;  
%     landmark = 1;
end
%%
% options.constraint_pos = verts(landmark,:)+10*rand(length(landmark),3);
options.constraint_pos = verts(landmark,:);

s_weights= ones(length(verts),1)*1;
c_weights= ones(length(landmark),1)*constraint_weight;
options.normalize = 0;
options.symmetrize = 1; 
options.type = 'conformal'; % the other choice is 'combinatorial', 'distance', spring, 'conformal', or, 'authalic', mvc, 

%%
options.method = 'hard';% 'hard', 'soft', 'bi-harmonic'
[newVertices, A] = compute_least_square_mesh(verts,faces,landmark,delta_coords,s_weights,c_weights,options);
figure('Name', options.method);set(gcf,'color','white');
trisurf(faces,newVertices(:,1),newVertices(:,2),newVertices(:,3)); axis off; axis equal; view3d rot;hold on;
scatter3(options.constraint_pos(:,1),options.constraint_pos(:,2), options.constraint_pos(:,3),100,'r','filled');

options.method = 'soft';% 'hard', 'soft', 'bi-harmonic'
[newVertices, A] = compute_least_square_mesh(verts,faces,landmark,delta_coords,s_weights,c_weights,options);
figure('Name', options.method);set(gcf,'color','white');
trisurf(faces,newVertices(:,1),newVertices(:,2),newVertices(:,3)); axis off; axis equal; view3d rot;

options.method = 'bi-harmonic';% 'hard', 'soft', 'bi-harmonic'
[newVertices, A] = compute_least_square_mesh(verts,faces,landmark,delta_coords,s_weights,c_weights,options);
figure('Name', options.method);set(gcf,'color','white');
trisurf(faces,newVertices(:,1),newVertices(:,2),newVertices(:,3)); axis off; axis equal; view3d rot;

