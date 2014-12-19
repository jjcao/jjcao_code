% ----------------------  failed!  ------------------------------
% almost all of the eigenvectors are not "symmetric" about 0.
%
%   Copyright (c) 2012 Junjie Cao

%% 
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_common'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
addpath ([MYTOOLBOXROOT '/jjcao_point'])
%%
DEBUG=0;
USE_SAVED_DATA=1;
ADJUST_SIGN=1;
if USE_SAVED_DATA
    load('wolf0_landmarks_43.mat')
else
    %% preprocess
    M.filename = [MYTOOLBOXROOT '/data/wolf0.off'];   
    [M.verts,M.faces] = read_mesh(M.filename);    
    M.nverts = size(M.verts,1);
    
    M.nlandmarks = round(M.nverts*0.01);
%     M.nlandmarks = M.nverts;    
    if M.nlandmarks < M.nverts % farest point sampling
        options.W = ones(M.nverts,1); % the speed function, for now constant speed to perform uniform remeshing
        M.landmarks = perform_farthest_point_sampling_mesh(M.verts,M.faces,[],M.nlandmarks,options);% perform the sampling of the surface
        M.landmarks = sort(M.landmarks);        
        if DEBUG
            figure('Name','Input model'); set(gcf,'color','white'); 
            h = plot_mesh(M.verts, M.faces);hold on;
            scatter3(M.verts(M.landmarks,1), M.verts(M.landmarks,2),M.verts(M.landmarks,3),40,'b','filled');
            view3d rot; lighting none;
        end        
        [M.D, M.landmark2all_geodesic_matrix] = compute_geodesic_distance_matrix(M.verts, M.faces, M.landmarks);
    else
        M.landmarks = 1:M.nverts;
        [M.D, M.landmark2all_geodesic_matrix] = compute_geodesic_distance_matrix(M.verts, M.faces, M.landmarks);
    end    
    [pathstr, name, ext] = fileparts(M.filename);
    save(sprintf('%s_landmarks_%d.mat', name, M.nlandmarks), 'M');
end

%% Construct low-dimensional embeddings (Classical MDS) %%%%%
ndims = min(10, M.nlandmarks);
[landmark_vec, val] = compute_isomap(M.D, ndims);
figure; plot(val);
%% 
if M.nlandmarks < M.nverts
    %% mds interpolation on the full set of points
    % x = 1/2 * (L^T) * ( delta_n-delta_x )
    vec = zeros(M.nverts,ndims);
    % transpose of embedding
    LT = landmark_vec'; 
    for i=1:ndims
        LT(i,:) = LT(i,:) / val(i);
    end
    deltan = mean(M.D,2);
    for i=1:M.nverts
        deltax = M.landmark2all_geodesic_matrix(:,i).^2;
        vec(i,:) = 1/2 * ( LT * ( deltan-deltax ) )';
    end
else
    vec = landmark_vec;
end
clear landmark_vec;

%% plot eigenvector 
for i=1:6
    figure; set(gcf,'color','white'); 
    options.face_vertex_color = vec(:,i);
    h = plot_mesh(M.verts, M.faces,options);
    view3d rot; lighting none;colormap(jet);
end
