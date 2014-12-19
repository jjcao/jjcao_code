% test_texturepatch
% 
% plot_mesh use patcht 
%
% Copyright (c) 2013 Junjie Cao

clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])

%%
load([MYTOOLBOXROOT 'data/texturedata.mat']); 
figure, patcht(FF,VV,TF,VT,I);axis off; mouse3d