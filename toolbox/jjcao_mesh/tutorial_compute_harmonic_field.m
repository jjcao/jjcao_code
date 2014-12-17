% tutorial_compute_harmonic_field
%
% 
% Copyright (c) 2013 Junjie Cao
%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh\geodesic'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'jjcao_plot'])

%% options
DEBUG=1;

%% input
[verts,faces] = read_mesh([MYTOOLBOXROOT 'data/wolf0.off']);
nverts = size(verts,1);
constraint_id = 1;
[D,S,Q] = perform_fast_marching_mesh(verts, faces, constraint_id);
[tmp, I] = max(D);
constraint_id(end+1) = I(1);
constraint_value = [0, 1];
    
if DEBUG
    figure('Name','constraints & geodesic distance'); set(gcf,'color','white');hold on;
    shading_type = 'interp';% 'flat';
    h=trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
        'FaceVertexCData', D,  'FaceColor',shading_type,'edgecolor','none'); axis off; axis equal; mouse3d;
    scatter3(verts(constraint_id,1), verts(constraint_id,2),verts(constraint_id,3),40,'r','filled');
    options = []; options.Nlines = 10; options.edgecolor = 'black';
    plot_contours(faces,verts,D,options);
end

%% compute
options = [];
options.type = 'conformal';%combinatorial;conformal;%spring
options.solver = 3;
tic
fid = compute_mesh_harmonic_field(verts, faces, constraint_id, constraint_value, options.type, options);
toc

%% output & plot
figure('Name','constraints & harmonic field'); set(gcf,'color','white');hold on;
shading_type = 'interp';% 'flat';
h=trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
    'FaceVertexCData', fid, 'FaceColor',shading_type,'edgecolor','none'); axis off; axis equal; mouse3d;
scatter3(verts(constraint_id,1), verts(constraint_id,2),verts(constraint_id,3),40,'r','filled');
options = []; options.Nlines = 20; options.edgecolor = 'black';
plot_contours(faces,verts,fid,options);