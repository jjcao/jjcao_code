%% % demo_random_walks_segmentation
% seeds are selected by demo_smo_segmentation

clear;clc;close all;
%MYTOOLBOXROOT='C:\jjcao_code\toolbox';
MYTOOLBOXROOT='../..';
addpath(genpath(MYTOOLBOXROOT));

DEBUG=1;

%% input

filename = 'result\wolf2_p199_seg11_som_0.010000.mat';% cube_f1200_p96, fandisk_p100,wolf0_p200
load(filename);
nface = size(M.faces,1);
%%
% if DEBUG
%     figure('Name','Patches');movegui('southwest'); set(gcf,'color','white');    
%     for i=1:M.npatch
%         options.face_vertex_color = zeros(nface,1);
%         options.face_vertex_color((M.face_patch==i))=1;
% %         options.face_vertex_color = M.face_patch;
%         h = plot_mesh(M.verts, M.faces, options);
%         colormap(jet(2));view3d rot;lighting none;
%         delete(h);
%     end
% end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.nbins = 10; % for feature histgram
M.nsegments = 16; % fandisk_p100 (15=>16)
M.USE_CONCAVE_WEIGHT = 0;
M.thresDist = 0.1; % let adjaceny matrix more sparse

%%
options.USE_CONCAVE_WEIGHT = M.USE_CONCAVE_WEIGHT;
options.verts = M.verts;
options.faces = M.faces;
options.seed_id = M.rank_set(1:min(length(M.rank_set),M.nsegments));
[A B] =compute_random_walk_graph(M.patch_adjancy,M.patch_normal,M.patch_curvature_hist,M.verts_between_patch,options);
% Distances less than ThresDist are set to 0. 
% if DEBUG
%     ind = find(A); tmp = full(A(ind));
%     sprintf('distances less thresDist(%f): %d', M.thresDist, sum(tmp<M.thresDist))
% end
% A(A<M.thresDist) = 0 ;
% P = inv(A) * B;
P = A\B;
[C cluster_id] = max(P,[],2);

%% run clustering from rank set.
segs = zeros(M.npatch,1);
for j=1:M.nsegments
    pid = find(cluster_id==j);      
    segs(pid) = j;
end
M.face_segments = segs(M.face_patch);

%% show result
% load armadillo_v4326_p40.mat;
figure('Name','Segment by SMO');movegui('southwest'); set(gcf,'color','white');
options.face_vertex_color = M.face_segments;
h = plot_mesh(M.verts, M.faces, options);
colormap(jet(M.nsegments));view3d rot;
% colormap(lines(M.nsegments));view3d rot;
lighting none;
%%
figure('Name','Segment by SMO'); movegui('southeast');set(gcf,'color','white');
% tmp = (M.nsegments:-1:1)';
tmp = circshift(1:M.nsegments,2)';
options.face_vertex_color = tmp(M.face_segments);
h = plot_mesh(M.verts, M.faces, options);
% colormap(jet(M.nsegments));view3d rot;
colormap(colorcube(M.nsegments));view3d rot;
lighting none;

%% save result
[pathstr, name, ext] = fileparts(filename);
save(sprintf('%s_seg%d_random_walks.mat', ['../../result/' name], M.nsegments), 'M');