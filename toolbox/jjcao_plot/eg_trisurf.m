% eg_trisurf
%
% If texture is needed, use Gabriel Peyre's plot_mesh in toolbox_graph 
%
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../'));

%% load a mesh
test_file = {['/data/fandisk.off'],['/data/wolf0.off'],[ '/data/catHead_v131.off']};
M.filename = test_file{3};
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);
M.edges = compute_edges(M.faces); 

%% display point cloud
figure('Name', 'point set');set(gcf,'color','white');hold on;
% scatter3(M.verts(:,1),M.verts(:,2), M.verts(:,3),10,'b','filled');
scatter3(M.verts(:,1),M.verts(:,2), M.verts(:,3),...
    10,M.verts(:,1),'filled');
axis off;    axis equal;   set(gcf,'Renderer','OpenGL');
camorbit(0,0,'camera'); axis vis3d; view(-90, 0); mouse3d

%% display mesh with specified face and edge color
figure('Name','specified face and edge color'); set(gcf,'color','white');hold off;
h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
    'FaceColor', 'cyan',  'edgecolor',[1,0,0], 'faceAlpha', 0.9); axis off; axis equal; mouse3d

%% display mesh with specified face color
figure('Name','specified face color'); set(gcf,'color','white');hold off;
h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
     'FaceColor', 'cyan', 'edgecolor','none'); axis off; axis equal; mouse3d
L1 = light('Position', [-1, -1, -1]);
L2 = light('Position', [1, 1, 1]);
lighting phong
set(h, 'AmbientStrength', 0.75);
set(h, 'DiffuseStrength', 0.5);


 %% display mesh with specified vertex color using flat shading
figure('Name','specified vertex color using flat shading'); set(gcf,'color','white');hold off;
shading_type = 'flat';
h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
    'FaceVertexCData', M.verts(:,1), 'FaceColor',shading_type, 'faceAlpha', 0.9); axis off; axis equal; mouse3d

%% display mesh with specified vertex color using interp shading
figure('Name','specified vertex color using interp shading'); set(gcf,'color','white');
% shading_type = 'interp';
% h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
%     'FaceVertexCData', M.verts(:,1), 'FaceColor',shading_type, 'faceAlpha', 0.9); axis off; axis equal; mouse3d
% set(h, 'edgecolor', 'none'); % cancel display of edge.

h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
    'FaceVertexCData', M.verts(:,1), 'edgecolor','none','faceAlpha', 0.9); axis off; axis equal; mouse3d
shading interp
colormap jet(256);colorbar
% colorbar('off');

%% display bundary vertices
boundary=compute_boundary(M.faces);
hold on;
scatter3(M.verts(boundary{:},1),M.verts(boundary{:},2), M.verts(boundary{:},3),50,'r','filled');
