% function test_differential_coords(fid)
% %% 
% if nargin < 1
%     fid = 1;
% end
clear;clc;close all;
addpath(genpath('../../'));

%% test 1, 1-ring vertices of center vertex are co-planar. Dcp has no
%% tangential component, and Combinatorial has tangential component.
fid = 2;
sel=1:7;%sel=3;

% %% test 2
% fid = 2;
% sel=1:8;

%%
test_file = {'hexagon_1.off', 'hexagon_2.off'};
[vertices,faces]=read_mesh(test_file{fid});
% [normalv, normalf] = compute_normal(vertices, faces);
% normalv = normalv';
[normalv, normalf] = compute_vertex_normal(vertices, faces, 0);

options.rings = compute_vertex_face_ring(faces);
options.normalize = 1; options.symmetrize=1;

%% dcp and combinatorial
dcpL = compute_mesh_laplacian(vertices, faces, 'dcp',options);%dcp, combinatorial
dcp_delta_coords = dcpL*vertices;
dcp_bc = (speye(size(dcpL)) - dcpL) * vertices;

cL = compute_mesh_laplacian(vertices, faces, 'combinatorial',options);
c_delta_coords = cL*vertices;
c_bc = (speye(size(cL)) - cL) * vertices;

%% display
figure;hold on;
plot_mesh(vertices, faces);

h = plot3(dcp_bc(sel,1),dcp_bc(sel,2), dcp_bc(sel,3), 'rx');
h = plot3(c_bc(sel,1),c_bc(sel,2), c_bc(sel,3), 'g.');
set(h, 'MarkerSize', 10);

h = quiver3(dcp_bc(sel,1),dcp_bc(sel,2),dcp_bc(sel,3),dcp_delta_coords(sel,1),...
    dcp_delta_coords(sel,2), dcp_delta_coords(sel,3),0, '-r');
h = quiver3(c_bc(sel,1),c_bc(sel,2),c_bc(sel,3),c_delta_coords(sel,1),...
    c_delta_coords(sel,2), c_delta_coords(sel,3),0, '-g');
h = quiver3(c_bc(sel,1),c_bc(sel,2),c_bc(sel,3),normalv(sel,1),...
    normalv(sel,2), normalv(sel,3),0, '-b');
view3d zoom;
sprintf('dcp diff coord: %d, %d, %d', dcp_delta_coords(3,:))
sprintf('combinatorial diff coord: %d, %d, %d', c_delta_coords(3,:))

% if fid == 1
%     assert( abs(dcp_delta_coords(3,1))<1e-04);
%     assert( abs(dcp_delta_coords(3,3))<1e-04);
%     assert( abs(c_delta_coords(3,1))>1e-04);
%     assert( abs(c_delta_coords(3,3))>1e-04);  
% end