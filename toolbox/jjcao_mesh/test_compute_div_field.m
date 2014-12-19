% function test_compute_div_field(fid)
% %% 
% if nargin < 1
%     fid = 1;
% end
clear all;
%close all;
%% test 1
fid = 1;
sel=1:7;%sel=3;
% %% test 2
% fid = 2;
% sel=1:8;

%%
test_file = {'hexagon_1.off', 'cube.off'};
[vertices,faces]=read_mesh(test_file{fid});

%% display mesh
figure('name', 'blue is laplacian, green is div');hold on;
plot_mesh(vertices, faces); 
options.rings = compute_vertex_face_ring(faces);
options.normals = compute_face_normals(vertices, faces);
options.areas = area_1_ring(vertices, faces, options.rings);
%% display, laplacian
options.normalize = 0; options.symmetrize=1;
L = compute_mesh_laplacian(vertices, faces, 'Laplace-Beltrami',options);%dcp, combinatorial,Laplace-Beltrami
delta_coords = diag(options.areas)*L*vertices;
NL = diag(diag(L))-L; 
nbc = diag(diag(L).^-1) * NL * vertices;
h = plot3(nbc(sel,1),nbc(sel,2), nbc(sel,3), 'gx');
set(h, 'MarkerSize', 10);
h = quiver3(nbc(sel,1),nbc(sel,2),nbc(sel,3),delta_coords(sel,1),...
    delta_coords(sel,2), delta_coords(sel,3),0, '-b');

%% display, div
div_coords = zeros(size(vertices));
for i = 1:3
    vfs = compute_grad_field(vertices(:,i), vertices, faces, options.normals);
    divs = compute_div_field(vfs, vertices, faces, options);
    div_coords(:,i) = divs;
end
h = quiver3(nbc(sel,1),nbc(sel,2),nbc(sel,3),div_coords(sel,1),...
    div_coords(sel,2), div_coords(sel,3),0,'-g');
view3d rot;

%% reconstruct mesh from delta_coords
handle_id = [3];
s_weights= ones(length(vertices),1)*1;
c_weights= ones(length(handle_id),1)*100;
options.method = 'soft';%soft or 'hard' or bi-harmonic
options.type = 'dcp'; % the other choice is 'combinatorial', 'distance', spring, 'conformal', mvc, 
nvs = compute_least_square_mesh(vertices,faces,handle_id,delta_coords,s_weights,c_weights,options);

figure('name', 'new mesh-laplacian');hold on;
plot_mesh(nvs, faces, options); colormap jet(256);
h = plot3(nvs(handle_id,1),nvs(handle_id,2), nvs(handle_id,3), 'r.');
set(h, 'MarkerSize', 10);
view3d rot;

%% reconstruct mesh from div
nvs = compute_least_square_mesh(vertices,faces,handle_id,div_coords,s_weights,c_weights,options);
figure('name', 'new mesh-div');hold on;
plot_mesh(nvs, faces, options); colormap jet(256);
h = plot3(nvs(handle_id,1),nvs(handle_id,2), nvs(handle_id,3), 'r.');
set(h, 'MarkerSize', 10);
view3d rot;