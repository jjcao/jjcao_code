%
% warning: isochart can not be obtained by random walk using developbility
% between patches as similarity.
%   the developbility between two patches may be high, even distance
% between normals are large!
%
% (C) Copyright jjcao, 2012

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../..';
addpath ([MYTOOLBOXROOT '/jjcao_mesh'])
addpath ([MYTOOLBOXROOT '/jjcao_point'])
addpath ([MYTOOLBOXROOT '/jjcao_io'])
addpath ([MYTOOLBOXROOT '/jjcao_plot'])
addpath ([MYTOOLBOXROOT '/jjcao_interact'])
addpath ([MYTOOLBOXROOT '/jjcao_common'])
addpath ([MYTOOLBOXROOT '/jjcao_mesh/geodesic'])
addpath ('div_rank')

DEBUG=1;
USEFILE=1;
%% input
filename = 'fandisk_p100.mat';% cube_f1200_p96, fandisk_p100,wolf0_p200
load(filename);
nface = size(M.faces,1);
if DEBUG
    figure('Name','Supervertex by NCut'); set(gcf,'color','white'); 
    options.face_vertex_color = M.face_patch;
%     options.alfa=0.5;
    h = plot_mesh(M.verts, M.faces, options);
%     set(h, 'edgecolor', 'none');
    colormap(jet(M.npatch)); 
    view3d rot; lighting none;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.alfa = 1;
M.nsegments = 20; % fandisk_p100 (15=>16)
M.USE_CONCAVE_WEIGHT = 0;
M.constZ = 0.01;% The contant for ground conductance0.01
M.thresDist = 0.1; % let adjaceny matrix more sparse
M.merge_option = 0; % 0: not merge; 1: merge segments by potential; 2: you have to set M.thresMerge
M.thresMerge = 0.1; % merge segments by reducing number of nsegments, thresMerge is threshold for difference of possibility of a face belong to different segments.

% Select a method to solve linear eq Ax=b.
% 0: incremental inverse. % 1: LU decomposition. % Others: naive backslash. 
M.LIN_MODE = 0 ;

%% compute adjacency matrix by inner product of patch features
if ~USEFILE
    [M.patch_adjancy,M.patch_centers,M.patch_verts,M.patch_faces, ...
        M.verts_between_patch] = compute_face_patch_graph(M.faces,M.face_patch,M.verts,M.npatch); % adjacency matrix A, A(i,i)=0
    if DEBUG
        figure('Name','patch_graph'); set(gcf,'color','white');
        h = plot_graph(M.patch_adjancy,M.patch_centers);axis equal;view3d rot;
    end
    M.A_iso=compute_patch_iso_similarity(M.patch_adjancy,M.verts,M.faces,M.patch_verts, M.patch_faces,M.face_patch);
    save(filename, 'M');
end
% %% i==6,下面；i=10,上平面
% if DEBUG
%     figure('Name','Segment by SMO');movegui('southwest'); set(gcf,'color','white');
%     for i=1:M.npatch
%         tmp = find(M.patch_adjancy(i,:)>0);
%         full([tmp; M.patch_adjancy(i,tmp); M.A_iso(i,tmp)])
%         for j=tmp
%             options.face_vertex_color = zeros(size(M.face_patch));
%             options.face_vertex_color(M.face_patch==i)=10;
%             options.face_vertex_color(M.face_patch==j)=20;
%             h = plot_mesh(M.verts, M.faces, options);
%             colormap(jet);view3d rot;
%             lighting none;        
%             delete(h);        
%         end
%     end
% end
%%
M.A_angle=compute_patch_angle_similarity(M.patch_adjancy,M.fnormal,M.patch_faces);
A = compute_patch_similarity(M.patch_adjancy,M.A_iso, M.A_angle, M.alfa);

%% Run Diversity Ranking and Clustering
% conductance to dummy ground
% z =  M.constZ*mean(A(A>0)); % simple Z
% Z = z*ones(M.npatch,1);

Z = std(A,[],2).*mean(A,2);% better Z for image segmentation (std(A,[],2).*mean(A,2);),
Z = M.constZ*Z*10;

%%
% Distances less than ThresDist are set to 0. 
if DEBUG
    ind = find(A); tmp = full(A(ind));
    sprintf('distances less thresDist(%f): %d', M.thresDist, sum(tmp<M.thresDist))
end
A(A<M.thresDist) = 0 ;
% Set diagonal elements to zeros
A(speye(size(A))~=0)=0 ;

% Each SP is weighted by its number of faces, todo should be replaced with
% area
massSP = zeros(M.npatch,1);
for i=1:M.npatch
    massSP(i) = length(M.patch_faces{i});
end

%% (1) Lazy greedy
[M.rank_set rank_set_val] = run_div_rank_lazy_greedy(A, M.nsegments, M.LIN_MODE, Z,1:M.npatch,massSP);
% (2) naive greedy
% [M.rank_set rank_set_val] = run_div_rank_naive_greedy(A, M.nsegments, M.LIN_MODE, Z,1:M.npatch,massSP);
tmp = diff(rank_set_val);
figure; plot(tmp,'DisplayName','tmp','YDataSource','tmp');

%% determine acutal segmentation number
options.merge_option=M.merge_option;
options.nclusters=M.nsegments; 
options.thresMerge=M.thresMerge; 
[cluster_id,cluster_info,probMat] = compute_actual_clusters(A,M.rank_set,rank_set_val,options);
% compute marginal temperature gain
% [obj_val rank_set_gt] = compute_obj_val(A, M.nsegments, z, M.LIN_MODE);
% tmp=sum(obj_val);
% tmp = diff(tmp);
M.nsegments = length(cluster_info);
%% run clustering from rank set.

sprintf('num of segments: %d, num of patches of each segment:', M.nsegments)
cluster_info'

segs = zeros(M.npatch,1);
for j=1:M.nsegments
    pid = find(cluster_id==j);      
    segs(pid) = j;
end
M.face_segments = segs(M.face_patch);
%%
% if DEBUG
%     figure('Name','Segment by SMO');movegui('southwest'); set(gcf,'color','white');    
%     for i=1:M.nsegments
%         options.face_vertex_color = zeros(nface,1);
%         options.face_vertex_color((M.face_segments==i))=1;
%         h = plot_mesh(M.verts, M.faces, options);
%         colormap(jet(2));view3d rot;lighting none;
%         delete(h);
%     end
% end    
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
save(sprintf('%s_seg%d_som_%f.mat', ['result/' name], M.nsegments,M.constZ), 'M');