clear;clc;close all;
extension='.off';
filename = 'E:\jjcao_paper\PCD_Orientation\result\3-2-1.58-25';

tic
P.filename = [filename extension];% point set
[P.pts,P.faces, P.normal1] = read_mesh(P.filename);
P.npts = size(P.pts,1);
P.pts = GS.normalize(P.pts);

figure('Name','Original point cloud');movegui('northwest');set(gcf,'color','white');
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),10,'.','MarkerEdgeColor',GS.PC_COLOR); axis off;axis equal; hold on;
view3d rot;