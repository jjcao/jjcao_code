%%
%
% jjcao @ 2013
%
clear;clc;close all;

in.inputdir = 'F:\jjcao_data\MSRA_1000/';%MSRA_1000, test_img,cvpr07supp
out.outputdir = 'F:\jjcao_data\MSRA_1000_result/';%MSRA_1000_result%test_img
in.imprefix = 'Imgs_';
IS_DEBUG = 1;
%%
for i = 1:1
    close all;
    in.name = sprintf('%s%04d', in.imprefix, i);
    in.img = imread([in.inputdir in.name '.jpg']);  
    if IS_DEBUG, figure;set(gcf,'color','white');movegui('northwest'); imshow(in.img); end
    
    
end
