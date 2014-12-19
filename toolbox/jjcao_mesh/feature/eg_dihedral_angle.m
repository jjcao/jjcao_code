% compute signed dihedral angle
% 
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../../'));

% verts = [0,0,0; 1,0,0; 1,1,0; 0,1,1];% ravine
verts = [0,0,0; 1,0,0; 1,1,0; 0,1,-1];% ridge
faces = [1,2,3;1,3,4];
[sda, n1, n2] = signed_dihedral_angle(verts)


figure; set(gcf,'color','white');hold on;
shading_type = 'flat';
h=trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
    'FaceVertexCData', [0,1,0;0,0,1], 'FaceColor',shading_type); axis off;axis equal; mouse3d;
drawAxis3d(1, 0.02);%Ox vector is red, Oy vector is green, and Oz vector is blue.
text('Position',verts(1,:),'String','v1','FontSize',30);
text('Position',verts(2,:),'String','v2','FontSize',30);
text('Position',verts(3,:),'String','v3','FontSize',30);
text('Position',verts(4,:),'String','v4','FontSize',30);

pos = [mean(verts(faces(1,:),:));mean(verts(faces(2,:),:))];
h = quiver3(pos(1, 1), pos(1, 2), pos(1, 3), ...
    n1(1, 1), n1(1, 2), n1(1, 3), 0.5,'g','LineWidth',4);
h = quiver3(pos(2, 1), pos(2, 2), pos(2, 3), ...
    n2(1, 1), n2(1, 2), n2(1, 3), 0.5,'b','LineWidth',4);
