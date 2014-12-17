function [concave_id,convex_id, concave_fun, concave_thres,convex_thres] = compute_concave_convex_vertices(verts,faces,options)
%   refer to: tvcg11_Mesh Segmentation with concavity-aware fields
%   num_concavenhb = 0 or 1, default is 1, tvcg11 use 0;
%   concave_fun: is for each edge
%   concave_id: is index for concave vertex
%   convex_id: is index for convex_id vertex
%
%   Copyright (c) 2012 Hui Wang, JJCAO

if nargin < 3
    options = [];
end

A = triangulation2adjacency(faces);
normal = compute_normal(verts, faces)';

[row, col] = find(A);
edge = verts(row, :) - verts(col, :);
edgeLength = sqrt(sum(edge .^ 2, 2));
edgeLengthR = sparse(1 : length(edgeLength), 1 : length(edgeLength), edgeLength .^ -1) * edge;
edgeNormal = normal(col, :) - normal(row, :);

concave_fun = dot(edgeLengthR,edgeNormal,2);

%%
concave_thres = getoptions(options, 'concave_thres', mean(concave_fun(concave_fun>0)));
convex_thres = getoptions(options, 'convex_thres', mean(concave_fun(concave_fun<0)));
num_concavenhb = getoptions(options, 'num_concavenhb', 1);

matrix = sparse(row, col, concave_fun > concave_thres);
concave_id = find(sum(matrix, 2) > num_concavenhb);
matrix = sparse(row, col, concave_fun < convex_thres);
convex_id = find(sum(matrix, 2) > num_concavenhb);