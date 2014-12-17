% A. M. Bronstein, M. M. Bronstein, R. Kimmel, Numerical geometry of
% non-rigid shapes, Springer, 2008
%
% http://tosca.cs.technion.ac.il/book
%
% TUTORIAL
% Demonstration of multidimensional scaling using SMACOF algorithm
%
% This tutorial is based on MDS RRE code, distributed as part of 
% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
%
% Adapted by jjcao
%
% (C) Copyright Michael Bronstein, 2008
% All rights reserved.

%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='./';
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
%%
options.method = 'rre';
options.cycles = 2;
options.iter = 10;

%% load and display man dataset
load([MYTOOLBOXROOT 'data/man1']);
figure(1), set(gcf,'color','white');
h = trisurf(man.TRIV,man.X,man.Y,man.Z);
shading interp, colormap([1 1 1]*0.9), lighting phong, camlight head
axis image; axis off; view3d rot; title('Man, pose 1');
drawnow; pause;

%man = init_surface(man);
options.X0 = [man.X,man.Y,man.Z];
[X_] = mds(man.D,options);

figure(2), set(gcf,'color','white');
h = trisurf(man.TRIV,X_(:,1),X_(:,2),X_(:,3)); 
shading interp, colormap([1 1 1]*0.9), lighting phong, camlight head
axis image; axis off; view3d rot;axis([-150 150 -150 150 -150 150])
title('Canonical form, pose 1');
drawnow;pause;

%%
load([MYTOOLBOXROOT 'data/man2']);
figure(3), set(gcf,'color','white');
h = trisurf(man.TRIV,man.X,man.Y,man.Z); 
shading interp, colormap([1 1 1]*0.9), lighting phong, camlight head
axis image; axis off;  view3d rot;title('Man, pose 2');
drawnow;pause;

%man = init_surface(man);
options.X0 = [man.X,man.Y,man.Z];
[X_] = mds(man.D,options);

figure(4), set(gcf,'color','white');
h = trisurf(man.TRIV,X_(:,1),X_(:,2),X_(:,3)); 
shading interp, colormap([1 1 1]*0.9), lighting phong, camlight head
axis image; axis off;  view3d rot;title('Canonical form, pose 2');
drawnow;