function [eigvector,eigvalue,A]=compute_Laplace_eigen(verts,faces,k,withAreaNormalization, adjustL)
%
% compute Laplacian matrix and its eigenstructure
%
%  eigvector:  ith each column in this matrix is the ith eigenfunction of the Laplace-Beltrami operator
%  eigvalue:  ith element in this vector is the ith eigenvalue of the Laplace-Beltrami operator
%  A:      ith element in this vector is the area associated with the ith vertex
%
% Copyright (c) 2013 Shuhua Li, Junjie Cao

if nargin < 4
    withAreaNormalization = 0;
    adjustL = 0;
end
if nargin < 5
    adjustL = 0;
end

%%
options.use_c_implementation = 1;
[L, D, W]= compute_mesh_laplacian(verts,faces,'conformal',options);%Manifold-harmonic,Mean_curvature, conformal, combinatorial
if adjustL
    L = D - 0.99*W;
end
%%
if withAreaNormalization
    A = vertex_area(verts, faces);
    tmp = sparse(1:length(A),1:length(A), A);
    [eigvector,eigvalue] = eigs(L,tmp,k,'sm');% compute the first k smallest eigenvalue
else
    A = [];
    [eigvector,eigvalue] = eigs(L,k,'sm');
end
eigvalue = diag(eigvalue);

%% The increasing order
[lamda, index] = sort(eigvalue);
eigvalue=lamda;
eigvector = eigvector(:, index);
