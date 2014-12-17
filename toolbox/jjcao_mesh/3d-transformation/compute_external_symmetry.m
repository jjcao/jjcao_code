%compute_exteranl_symmetry

%Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../../'));
%% input mesh
filename='wolf0.off';
[verts,faces] = read_mesh(filename);
nverts = size(verts,1);
faceArea = compute_area_faces(verts, faces);
verts= verts/ sqrt(sum(faceArea)); %πÈ“ªªØ
%% display input mesh
 figure('name','input'); set(gcf,'color','white');hold on;
 h=trisurf(faces, verts(:,1), verts(:,2), verts(:,3), 'FaceColor', [0.1,0.1,0.1],  'faceAlpha', 0.4); axis off; axis equal; mouse3d
 set(h, 'edgecolor', 'none');     h=[];
%% compute reflect points
center = mean(verts);
normal = [1,0,0];
drawVector3d(center, normal, 'b');
scatter3(center(:,1),center(:,2), center(:,3),100,'r','filled');
normal = normal./sqrt(dot(normal,normal));

plane0 = createPlane(center, normal);
drawPlane3d(plane0, 'g');

verts_reflect = reflect_points(verts, center, normal);
%% compute symmetry points on the surface
tic

tree = kdtree_build(verts);
[symIdx,dist]=kdtree_nearest_neighbor(tree,verts_reflect);
kdtree_delete(tree);

verts_sym=verts(symIdx,:);

toc
error_dist=sum(dist);
%% display sparse symmetry 
idx=randi(size(verts,1),50,1);
sparseIdx=[idx,symIdx(idx)];
h = plot_edges(sparseIdx, verts);
scatter3(verts(idx,1),verts(idx,2), verts(idx,3),50,'r','filled');
scatter3(verts(symIdx(idx),1),verts(symIdx(idx),2), verts(symIdx(idx),3),50,'r','filled');
%% save
sym.filename=filename;  sym.symIdx=symIdx;  sym.dist=dist; 
[pathstr,name,ext]=fileparts(filename);
save(sprintf('%s_external_sysmetry.mat',name),'sym');