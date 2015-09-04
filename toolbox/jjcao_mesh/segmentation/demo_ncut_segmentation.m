clear;clc;close all;
%MYTOOLBOXROOT='C:\jjcao_code\toolbox';
MYTOOLBOXROOT='../..';
addpath(genpath(MYTOOLBOXROOT));

DEBUG=1;
%% input
filename = 'fandisk_p200.mat';% cube_f1200_p96, fandisk_p100,wolf0_p200
load(filename);
nface = size(M.faces,1);
if DEBUG
    figure('Name','Supervertex by NCut'); set(gcf,'color','white'); 
    options.face_vertex_color = M.face_patch;
    h = plot_mesh(M.verts, M.faces, options);
%     set(h, 'edgecolor', 'none');
    colormap(jet(M.npatch)); 
    view3d rot; lighting none;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.nbins = 10; % for feature histgram
M.nsegments = 14;
M.USE_CONCAVE_WEIGHT = 0;
M.thresDist = 0.1; % let adjaceny matrix more sparse

%% compute adjacency matrix by inner product of patch features
if ~isfield(M,'patch_adjancy')
    [M.patch_adjancy,M.patch_centers,M.patch_verts,M.patch_faces, ...
        M.verts_between_patch] = compute_face_patch_graph(M.faces,M.face_patch,M.verts,M.npatch); % adjacency matrix A, A(i,i)=0
    if DEBUG
        figure('Name','patch_graph'); set(gcf,'color','white');
        h = plot_graph(M.patch_adjancy,M.patch_centers);axis equal;view3d rot;
    end    
    M.patch_normal = compute_patch_angle(M.fnormal,M.patch_faces);
%     M.patch_curvature_hist = compute_patch_curvature_hist(M.verts,M.faces, M.patch_verts, M.nbins);
    save(filename, 'M');
end

options.USE_CONCAVE_WEIGHT = M.USE_CONCAVE_WEIGHT;
options.verts = M.verts;
options.faces = M.faces;
A=compute_random_walk_graph(M.patch_adjancy,M.patch_normal,M.patch_curvature_hist,M.verts_between_patch,options);

% Distances less than ThresDist are set to 0. 
if DEBUG
    ind = find(A); tmp = full(A(ind));
    sprintf('distances less thresDist(%f): %d', M.thresDist, sum(tmp<M.thresDist))
end
A(A<M.thresDist) = 0 ;
% Set diagonal elements to zeros
A(speye(size(A))~=0)=0 ;

%% NCut
tic;
[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A,M.nsegments);
toc;

segs = zeros(M.npatch,1);
for j=1:M.nsegments,
    pid = find(NcutDiscrete(:,j));      
    segs(pid) = j;
end
M.face_segments = segs(M.face_patch);

%% show result
% load armadillo_v4326_p40.mat;
figure('Name','Segment by NCut');movegui('southwest'); set(gcf,'color','white');
options.face_vertex_color = M.face_segments;
h = plot_mesh(M.verts, M.faces, options);
colormap(jet(M.nsegments));view3d rot;
% colormap(lines(M.nsegments));view3d rot;
lighting none;
%%
figure('Name','Segment by Ncut'); movegui('southeast');set(gcf,'color','white');
tmp = (M.nsegments:-1:1)';
options.face_vertex_color = tmp(M.face_segments);
h = plot_mesh(M.verts, M.faces, options);
% colormap(jet(M.nsegments));view3d rot;
colormap(colorcube(M.nsegments));view3d rot;
lighting none;

%% save result
[pathstr, name, ext] = fileparts(filename);
save(sprintf('%s_seg%d_ncut.mat', ['../../result/' name], M.nsegments), 'M');