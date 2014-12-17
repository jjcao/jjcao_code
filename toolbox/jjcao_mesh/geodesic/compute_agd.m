function [agd,landmark]=compute_agd(verts,faces,nbr_landmarks,DEBUG,normalize)
%% Compute average geodesic distance(agd) of surface,
% Based on the paper "Topology Matching for Fully Automatic Similarity 
% Estimation of 3D Shapes",but different in computing b_i and area(b_i).
%
% agd(v)=sum｛g(v,b_i)*area(b_i)｝
% where b={b_1,...,b_n)} are the base vertices which are sccttered almost equally on the surface S; 
% area(b_i) is the area that b_i occupied, and sum｛area(b_i)｝equals area(S); 
% g(v,b_i) is the geodesic distance between v and b_i on S。
% 
% Input: verts
%        faces
%        nbr_landmark， the number of landmarks
%        DEBUG， 1,view the voronoi cells and agd cells; 0,otherwise;
%        normalize，1，normalize agd to 0-1；0,otherwise
%         
% Output: agd
%         landmark, the index of b={b_1,...,b_n)}
%
% Copyright (c) 2013 Shuhua Li, Junjie Cao

if nargin < 4
    DEBUG = 0;
    normalize = 0;
end
if nargin < 5
    normalize = 0;
end

%%
[face_patch, patch_area,landmark]=supervertex_by_farthest_sampling(verts,faces,nbr_landmarks);

%% compute agd(v)=sum｛g(v,b_i)*area(b_i)｝
D=zeros(size(verts,1),nbr_landmarks); %D(:,i)为所有点到b_i测地距
options.nb_iter_max = Inf; % trying with a varying number of iterations: 100, 1000, Inf
for i=1:nbr_landmarks
    [D(:,i),S,Q] = perform_fast_marching_mesh(verts, faces, landmark(i), options);
end
agd=D*patch_area;  %agd(v)=对i求和｛g(v,b_i)*area(b_i)｝

%% 归一化
if normalize
    agd=(agd-min(agd))/max(agd);
end
