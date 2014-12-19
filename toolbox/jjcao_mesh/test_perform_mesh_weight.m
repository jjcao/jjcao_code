% test_perform_mesh_weight
%
% 
% Copyright (c) 2013 Junjie Cao
%%
clc;clear all;close all;
addpath(genpath('../'));

% mex -g -largeArrayDims -I"../../include/eigen-3.1.3" perform_mesh_weight.cpp
% mex -g -largeArrayDims -I"D:/eigen-3.1.3" perform_mesh_weight.cpp
% mex -g -largeArrayDims -DNAN_EQUALS_ZERO -I"D:/eigen-3.1.3" perform_mesh_weight.cpp

%%
[verts,faces] = read_mesh('cube_f300.off');%regular_triangle % cube_f300
options.use_c_implementation = 1;
L1 = compute_mesh_laplacian(verts,faces,'combinatorial',options);
options.use_c_implementation = 0;
L10 = compute_mesh_laplacian(verts,faces,'combinatorial',options);
if sum(sum(L1-L10))
    tmp = full(max(max(abs(L1-L10))));
    sprintf('max difference is: %f', tmp)
end

options.use_c_implementation = 1;
L3 = compute_mesh_laplacian(verts,faces,'conformal',options);
options.use_c_implementation = 0;
L30 = compute_mesh_laplacian(verts,faces,'conformal',options);
if sum(sum(L3-L30)) > 1e-10
    tmp = full(max(max(abs(L3-L30))));
    sprintf('max difference is: %f', tmp)
end

%%%%%%%%%%%%%%%%%%
options.use_c_implementation = 1;
L2 = compute_mesh_laplacian(verts,faces,'laplace-beltrami',options);
options.use_c_implementation = 0;
L20 = compute_mesh_laplacian(verts,faces,'laplace-beltrami',options);
if sum(sum(L2-L20)) > 1e-10
    tmp = full(max(max(abs(L2-L20))));
    sprintf('max difference is: %f', tmp)
end

options.use_c_implementation = 1;
L4 = compute_mesh_laplacian(verts,faces,'manifold-harmonic',options);
options.use_c_implementation = 0;
L40 = compute_mesh_laplacian(verts,faces,'manifold-harmonic',options);
if sum(sum(L4-L40)) > 1e-10
    tmp = full(max(max(abs(L4-L40))));
    sprintf('max difference is: %f', tmp)
end
