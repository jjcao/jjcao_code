function A = compute_face_features(verts,faces,fnormal,USE_CONCAVE_WEIGHT)
%
% computation of similarity of two faces: refer to tog12_Variational Mesh Decomposition
% computation of concave: refer to tvcg11_Mesh Segmentation with concavity-aware fields
%
%   Copyright (c) 2012 JJCAO
DEBUG=false;
nface = size(faces,1);

%% compute adjacency matrix by inner product of face normals
if USE_CONCAVE_WEIGHT
    [A fverts E1 E2] = compute_dual_graph(faces, verts, true); % adjacency matrix A, A(i,i)=0
else
    [A fverts] = compute_dual_graph(faces, verts, false);
end

if DEBUG
    figure('Name','Dual Graph'); set(gcf,'color','white');
    options.ps = 1;
    h = plot_graph(A,fverts,options);axis equal; hold on;view3d rot;    
end

% compute eta by concave & convex edge
ind = find(A);
if USE_CONCAVE_WEIGHT
    eta = compute_eta(verts, faces, E1, E2, ind, DEBUG)*0.5;
else
    eta = ones(length(ind), 1)*0.5;
end

%
[row,col] = ind2sub(size(A),ind);
d1 = eta.*sum( (fnormal(row,:) - fnormal(col,:)).^2, 2);
w = exp(-d1./mean(d1)*4);
A = sparse(row,col,w);
% A = A + speye(nface,nface);
return;
