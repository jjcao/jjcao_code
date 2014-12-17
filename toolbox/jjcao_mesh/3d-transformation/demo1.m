%
%  
%
% Copyright (c) 2013 Junjie Cao
clc;clear all;close all;
addpath(genpath('../../'));
DEBUG = 1;

%% input
p1 = [-0.32849099043605101,0.35927207722296428, 0.016297145120710760];
p2 = [-0.32947573275434283,	0.40803202649574011, 1.3997502849072203];
p3 = [-0.53508738872598993,0.47649070331453369, -0.0059969902861965602];
p4 = [-0.53266738872598984,	0.47602370331453370, -1.28882397685329653];

v1 = [p1;p2];
v2 = [p3;p4];
c1 = mean(v1);
c2 = mean(v2);
edge1 = [p1 p2];
edge2 = [p3 p4];
edge = [c1 c2];

% prepare a figure for drawing
figure(1); clf; hold on;axis off; axis equal; mouse3d;set(gcf, 'renderer', 'opengl');
drawEdge3d(edge1, 'color', 'r', 'linewidth', 2);
drawEdge3d(edge2, 'color', 'r', 'linewidth', 2);
drawEdge3d(edge, 'color', 'g', 'linewidth', 2);

n1 = cross( p1 - p2, c1 - c2);
n2 = cross( p3 - p4, c1 - c2);
if norm(n1) > norm(n2)
    n = n1;
else
    n = n2;
end

plane0 = createPlane(c1, n);
drawPlane3d(plane0, 'g');

p5 = [-0.53906522021842873,	-0.46784248574973292,-0.0056450998623272935];
scatter3(p5(:,1),p5(:,2), p5(:,3),50,'g','filled');

%%
c = 0.5*(c1+c2);
figure('Name','input'); set(gcf,'color','white');hold on;
scatter3(v1(:,1),v1(:,2), v1(:,3),50,'g','filled');
scatter3(c1(:,1),c1(:,2), c1(:,3),100,'g','filled');
scatter3(v2(:,1),v2(:,2), v2(:,3),50,'b','filled');
scatter3(c2(:,1),c2(:,2), c2(:,3),100,'b','filled');
scatter3(c(:,1),c(:,2), c(:,3),100,'r','filled');
axis off; axis equal; mouse3d;

%% rotate verts2 around center Axis 
center = [0.00318914,0.223123,0.819125];
direction =[-0.0020372,0.10127,0.994856];
angle = 3.14159;

figure('Name','c++'); set(gcf,'color','white');hold on;
scatter3(v1(:,1),v1(:,2), v1(:,3),50,'g','filled');
scatter3(c1(:,1),c1(:,2), c1(:,3),100,'g','filled');
scatter3(v2(:,1),v2(:,2), v2(:,3),50,'b','filled');
scatter3(c2(:,1),c2(:,2), c2(:,3),100,'b','filled');
scatter3(center(:,1),center(:,2), center(:,3),100,'r','filled');
axis off; axis equal; mouse3d;

%%
line = [center direction];
rot = createRotation3dLineAngle(line, angle);
[axisR angle2] = rotation3dAxisAndAngle(rot);
angle2

newvertices = transformPoint3d(v2,rot);
sum( sqrt(sum((v1 - newvertices).^2,2)) )/size(v1,1)

scatter3(newvertices(:,1),newvertices(:,2), newvertices(:,3),50,'r','filled');