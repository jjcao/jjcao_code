clear;clc;close all;
path(path,'toolbox') ;
%%
filename='../data/horse_v11912.off';
[M.verts,M.faces] = read_mesh(filename);
[M.normals] = compute_vertex_normal(M.verts,M.faces, true);
figure;set(gcf,'color','white');movegui('northwest');set(gcf,'Renderer','OpenGL');
hold on; plot_mesh(M.verts, M.faces);
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(00,90);view3d rot;

bbox = [min(M.verts(:,1)), min(M.verts(:,2)), min(M.verts(:,3)), max(M.verts(:,1)), max(M.verts(:,2)), max(M.verts(:,3))];
rs = bbox(4:6)-bbox(1:3);
radius = sqrt(dot(rs,rs))*0.5;

type = 'random';
% type = 'gaussian';
 p3 = 3*0.01;
%p3 = radius/1280*1.35;%radius/1280 too much of noise
% type = 'salt & pepper';
% p3 = 0.01;%0.01 is too bad

overts = pcd_noise(M.verts, type, p3);
figure;set(gcf,'color','white');movegui('northeast');set(gcf,'Renderer','OpenGL');hold on; 
options.face_vertex_color = GS.PC_COLOR;
h = plot_mesh(overts, M.faces, options);
set(h, 'edgecolor', 'none');
set(h, 'FaceColor', GS.PC_COLOR);

%colormap jet(256);
colorbar('off');
%scatter3(overts(:,1),overts(:,2),overts(:,3),60,'.b');hold on;
axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(00,90);view3d rot;

ofilename = sprintf('%s_%s_%f.off',filename(1:(end-4)), type, p3);
write_mesh(ofilename, overts,M.faces,[]);
ofilename = sprintf('%s_%s_%f.xyz',filename(1:(end-4)), type, p3);
write_mesh(ofilename, overts,[],M.normals);