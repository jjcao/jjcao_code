%
%  
%
% Copyright (c) 2013 Junjie Cao
clc;clear all;close all;
addpath(genpath('../../'));
DEBUG = 1;

%% input
[parts, vertices]= read_parts_obj('CS2b.obj');
figure('Name','input'); set(gcf,'color','white');hold on;
for i=4:9
    faces1 = parts(i).faces;
    verts1 = vertices(min(faces1):max(faces1),:);
    h=trisurf(faces1,vertices(:,1),vertices(:,2),vertices(:,3), 'FaceVertexCData', vertices(:,3));     
end
set(h, 'edgecolor', 'none');colormap jet(256); axis off; axis equal; mouse3d;

%% rotate verts2 around centerAxis 
center = [-0.0268528,-0.038966,0.0378716];
direction = [-0.016064,-0.0371969,0.999179];
line = [center direction];
drawLine3d(line);