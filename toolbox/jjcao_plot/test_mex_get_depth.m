% test_mex_get_depth
% mex mex_get_depth.cpp -lOpenGL32


clear;clc;close all;
addpath(genpath('../'));

test_file = {['/data/fandisk.off'],['/data/wolf0.off'],[ '/data/catHead_v131.off']};
M.filename = test_file{3};
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);
[M.verts, newlen] = normalize_vertex3d(M.verts, 1);

figure('Name','specified face and edge color'); set(gcf,'color','white');hold off;
h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
    'FaceColor', 'cyan','edgecolor','none'); axis off; axis equal; mouse3d
set(gcf,'Renderer','OpenGL');mouse3d;
set(gca, 'CameraTarget', [0,0,0]);
set(gca, 'Projection', 'orthographic');
set(gca,'cameraViewAngle',30);
set(gca, 'CameraPosition', [1,1,1]);%set(gca, 'CameraPosition', [0.5,0.5,0.5]);
figure(1);

%%
depthData=mex_get_depth;
figure
imshow(depthData);

binData = (depthData<1);
imshow(binData);

figure;
[B,L] = bwboundaries(binData);
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

