%
% Implementatoin of: tog14_Mesh Saliency via Spectral Processing
%
% now singal scale saliency seems correct for wolf2_v940, but failed for wolf2_v1565
% and the multi-one seems wrong!
%
% There are some details not todl by the paper. 
% I comment them with "jjcao: pay attention!!"
%
%jjcao @ 2014
% 

clc;clear all;close all;
addpath(genpath('../../../'));
options.bDebug = 1;
options.bSymmetrize=1;
options.bNormalize=0;
options.diagLength = 800;
options.n = 9;
epsilon = 0.002 * options.diagLength; % see the paper
t = (1:1:5)*epsilon^2; % see the paper: (1:1:5)*epsilon^2; % 5+(1:2:10)*epsilon^2£»
dist_const = 2.5; % see the paper

%% input
M.filename = 'Bimba_cvd_30K_R22.off';
M.filename = 'wolf2_v753.off';
% 85_v1089.off, 85_v1814.off, wolf2_v753.off
% armadillo_v502.off, gargoyle_v502.off,dragon_v1257.off,v2.off
% sphere1.obj,cube_f1200.off,torus_v500.off
% Bimba_cvd_30K_R22_reduce1001.off
[M.verts,M.faces] = read_mesh(M.filename); M.nverts = size(M.verts,1);
[M.verts,options.diagLength] = normalize_vertex3d(M.verts,options.diagLength);
M.nverts = size(M.verts, 1);
tree = kdtree_build(M.verts);

%% compute k by eq 17, not understand it still!
A = triangulation2adjacency(M.faces); % adjacency matrix 
ind = find(A>0);
[I, J] = ind2sub(size(A), ind);
dist2 = sum((M.verts(I,:) - M.verts(J,:)).^2, 2);
W = sparse(I,J, 1.0\(dist2+0.0000001) );
d = sum(W,2);
M.L = speye(M.nverts) - diag(d.^(-1)) * W; % bSymmetrize==0 && bNormalize==1

c = mean(sqrt(dist2));% average distance of all edges
sumDistPerVert = sum( sparse(I,J, sqrt(dist2)), 2);
k = c*sum(A,2)./sumDistPerVert + 1;

j = 1;
for i = t
    % Mt
    [verts1, nneigh] = gaussian_smoothing(M.verts, M.verts, dist_const*sqrt(i)*ones(M.nverts,1), i*ones(M.nverts,1), tree);    
    % SMt: saliency of Mt
    SMt = log_spectral_saliency(verts1, M.faces, options);
    % Mkt
    [verts2, nneigh] = gaussian_smoothing(M.verts, M.verts,dist_const*sqrt(i*k), i*k, tree);
    % SMkt: saliency of Mkt
    SMkt = log_spectral_saliency(verts2, M.faces, options);   
    % absolute difference of them
    S(:,j) = abs(SMkt - SMt);
    j = j + 1;
    sprintf('mean of |neighbors|: %f, if it is too small, it means that sigma2 is too low!', mean(nneigh))
    % map back to original mesh
    % ...    
end

Saliency = sum(S,2);
[Saliency, nneigh] = gaussian_smoothing(M.verts, Saliency, dist_const*sqrt(t(end)*k), t(end)*k, tree);
logSaliency = log(Saliency);
Saliency = (Saliency - min(Saliency))/(max(Saliency) - min(Saliency));
if options.bDebug
    figure;
    trisurf(M.faces,verts1(:,1),verts1(:,2),verts1(:,3), ...
    'FaceVertexCData', Saliency, 'FaceColor','interp','edgecolor', 'none');
    axis off;axis equal;colorbar;mouse3d; 
    figure;
    trisurf(M.faces,verts1(:,1),verts1(:,2),verts1(:,3), ...
    'FaceVertexCData', logSaliency, 'FaceColor','interp','edgecolor', 'none');
    axis off;axis equal;colorbar;mouse3d;     
end 
kdtree_delete(tree);

%% smooth Saliency with Laplaican smoothing 10 times?
for i = 1:10
    Saliency = M.L * Saliency;
%     Saliency = M_orign.L * Saliency;
end
logSaliency = log(Saliency);
if options.bDebug
    figure;
    trisurf(M.faces,verts1(:,1),verts1(:,2),verts1(:,3), ...
    'FaceVertexCData', Saliency, 'FaceColor','interp','edgecolor', 'none');
    axis off;axis equal;colorbar;mouse3d; 
    figure;
    trisurf(M.faces,verts1(:,1),verts1(:,2),verts1(:,3), ...
    'FaceVertexCData', logSaliency, 'FaceColor','interp','edgecolor', 'none');
    axis off;axis equal;colorbar;mouse3d;     
end 






