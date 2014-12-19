function [A B]=compute_random_walk_graph(A,patch_normal,patch_curvature_hist,verts_between_patch,options)
%
% computation of similarity of two faces: refer to tog12_Variational Mesh Decomposition
% computation of concave: refer to tvcg11_Mesh Segmentation with concavity-aware fields
%
%   Copyright (c) 2012 JJCAO

DEBUG=false;
USE_CONCAVE_WEIGHT=getoptions(options,'USE_CONCAVE_WEIGHT',0);
verts = getoptions(options,'verts', []);
faces = getoptions(options,'faces',[]);
seed_id = getoptions(options,'seed_id',[]);
npatch = length(A);
%% % use curvature hist
% tmp = sum( (M.patch_curvature_hist(row,:) - M.patch_curvature_hist(col,:)).^2, 2);
% A = sparse(row,col,exp(-0.5*tmp/mean(tmp)));

%% use patch normal
% compute eta by concave & convex verts shared by two patches
ind = find(A);
if USE_CONCAVE_WEIGHT
    eta = compute_eta(verts, faces, verts_between_patch, ind, DEBUG)*0.5;
else
%     eta = ones(length(ind), 1)*0.5;
    eta = ones(length(ind), 1);
%     tmp = A(ind);    
%     eta = max(tmp)./tmp; % eta=tmp/max(tmp);
%     eta = eta./max(eta)*0.5;    
end

%%
[row,col] = ind2sub(size(A),ind);
d1 = eta.*sum( (patch_normal(row,:) - patch_normal(col,:)).^2, 2);
w = exp(-d1./mean(d1)*4); % w = exp(-d1./mean(d1)*4);% for cad models
A = sparse(row,col,w);
%%
if isempty(seed_id)
    B = [];
else    
    A=sparse(1:npatch,1:npatch,-1.0./sum(A,2))*A;    
    A = A + speye(size(A));
    B =  -A(:,seed_id);
    B(seed_id,1:length(seed_id)) = 0;
%     B = zeros(size(A,1),length(seed_id));
    for i=1:length(seed_id)
        B(seed_id(i),i) = 1;
        A(seed_id(i),:) = 0;
        A(seed_id(i),seed_id(i)) = 1;
    end
end
return;

%%
function eta = compute_eta(verts, faces, verts_between_patch, ind, DEBUG)
%
% % bad
% concave_weight = 0.1;
% convex_weight = 0.2;
% normal_weight = 1;

% good
concave_weight = 1;
convex_weight = 0.8;
normal_weight = 0.1;

%%
eta = ones(length(ind), 1)*normal_weight;
[concave_id,convex_id] = compute_concave_convex_vertices(verts,faces);
concidx = zeros(size(verts,1),1);
convidx = concidx;
concidx(concave_id)=1;
convidx(convex_id)=1;
for i=1:length(ind);
    tmp = verts_between_patch{ind(i)};
    tmpconcave =sum(concidx(tmp));
    tmpconvex = sum(convidx(tmp));
    eta(i)=(tmpconcave*concave_weight+tmpconvex*convex_weight+length(tmp)-tmpconcave-tmpconvex)/length(tmp);
end

if DEBUG    
    figure('Name','Concave vertices'); set(gcf,'color','white');hold on;
    h = plot_mesh(verts,faces);
    h=scatter3(verts(concave_id,1),verts(concave_id,2),verts(concave_id,3),80,'b','filled');
    h=scatter3(verts(convex_id,1),verts(convex_id,2),verts(convex_id,3),80,'r','filled');
    view3d zoom;
end