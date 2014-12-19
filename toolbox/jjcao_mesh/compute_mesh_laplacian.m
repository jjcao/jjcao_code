function [L, D, W] = compute_mesh_laplacian(vertex,face,type,options)

% compute_mesh_laplacian - compute a laplacian matrix L, L = -(Laplacian operator).
%
%   L = compute_mesh_laplacian(vertex,face,type,options);
%
%   If options.symmetrize=1 and options.normalize=0 then 
%       L = D-W
%   If options.symmetrize=1 and options.normalize=1 then 
%       L = eye(n)-D^{-1/2}*W*D^{-1/2}
%   If options.symmetrize=0 and options.normalize=1 then 
%       L = eye(n)-D^{-1}*W.
%   where D=diag(sum(W,2)) and W is the unormalized weight matrix 
%   (see compute_mesh_weight).
%
%   type can 'combinatorial', 'distance', 'Laplace-Beltrami', 'Mean_curvature', 'conformal', 'MVC'.
%
%   See also compute_mesh_weight.
%
%   Changed by JJCAO
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
normalize = getoptions(options, 'normalize', 0);
symmetrize = getoptions(options, 'symmetrize', 1);
ttype = type;
switch lower(type)
    case {'combinatorial','graph'}
        ntype = 0;
    case 'distance'
        ntype = 1;
    case 'spring'
        ntype = 2;        
    case {'conformal','dcp','mean_curvature', 'laplace-beltrami','manifold-harmonic'}
        ntype = 3;      
        ttype = 'dcp';
    case {'mvc'}
        ntype = 6;
    otherwise
        error('Unknown type.')
end
    
use_c_implementation = getoptions(options, 'use_c_implementation', 1);% use fast C-coded version if possible
if exist('perform_mesh_weight', 'file') && use_c_implementation
    W = perform_mesh_weight(vertex',face',ntype, options);
else
    W = compute_mesh_weight(vertex,face,ttype,options);
end

switch lower(type)
    case {'mean_curvature', 'laplace-beltrami'}
        vert_areas = vertex_area(vertex,face);
        W = diag(vert_areas)*W;
    case {'manifold-harmonic'}
        vert_areas = vertex_area(vertex,face);
        [I,J,V]  = find(W);
        V = V ./ sqrt(vert_areas(I).*vert_areas(J));
        W = sparse(I,J,V);
end    
%%
D = diag(sum(W,2));
if symmetrize==1 && normalize==0
    L = D - W;
elseif symmetrize==1 && normalize==1
    L = speye(n) - D.^(-1/2) * W * D.^(-1/2);
elseif symmetrize==0 && normalize==1
    L = speye(n) - D.^(-1) * W;
else
    error('Does not work with symmetrize=0 and normalize=0');    
end