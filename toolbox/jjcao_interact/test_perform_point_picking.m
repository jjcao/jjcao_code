clc;clearvars;close all;
addpath ../jjcao_plot;

PointCloud = rand(3,100)*100;
perform_point_picking( PointCloud );
%view3d rot; hold on;
