clear;
clc;
close all;
%% loadind
I = imread('D:\MSRA_1000\Imgs_0104.jpg');
%labels = ReadDAT([size(I,1),size(I,2)],'D:\Segmentation\superpixel\super_pixel\dat\0_0_272_SLIC.dat');
%load('D:\MSRA_1000\Imgs_LabelsFine_0352.mat');
%count = getCount(labels);
%%
%% features
img = I;
pfeats.ratio = 0.9;
pfeats.binnum = 12;
img = fix(im2double(img)*255);
[row,col,depth] = size(img);
img = reshape(img,1,row*col*depth);
spImg = ones(1,row*col);
numPatches = 1;
pixelNumEachPatch = row*col;
% [L,a,b] = RGB2Lab(I);
% I(:,:,1) = L;
% I(:,:,2) = a;
% I(:,:,3) = b;
% words
[pfeats.pixelNumRGB, pfeats.R, pfeats.G, pfeats.B, pfeats.colorNumEachPatch] = ...
    quantize_patch_color_no_optimization(img, spImg, pixelNumEachPatch, numPatches, pfeats.ratio, pfeats.binnum);      
R = pfeats.R;
G = pfeats.G;
B = pfeats.B;
colorNumEachPatch = pfeats.colorNumEachPatch;
pixelNumRGB = row*col;
% 
% Image = reshape(img,1,row*col*depth);
% numPatches = 1;
% ratio = 0.9;
%  binnum = 12;
%   bin = fix((256-0.1)/binnum);
% [pixelNumRGB, R, G, B, colorNumEachPatch,showImage] = show_patchHistogram(Image, spImg, pixelNumEachPatch, numPatches, ratio, binnum);  
% 
% sp = showImage;
% showImage =(showImage+1)*bin;
% showImage = im2uint8(showImage./255);
% showImage = reshape(showImage,row,col,depth);
% figure(1)
% imshow(showImage,[]);
binnum = 12;
bin = fix((256-0.1)/binnum);
img = fix(img./bin);
sp = nearcolor(img, R, G, B,colorNumEachPatch,row,col,depth);
img = fix(sp.*(bin));
img = im2uint8(img./255);
img = reshape(img,row,col,depth);
figure(1);
imshow(img,[]);
%% document
documentNum = 500;
 X = randomDocument(documentNum,sp, R, G, B,row, col,depth)

%% plsa
topic = 2;


%[pz pdz pwz pzdw]=plsa(X,topic);
[pwz,pdz,pz,Li] = pLSA_EM(X,2);

% Learn.Verbosity = 1;
% Learn.Max_Iterations = 200; 
% Learn.heldout = .1; % for tempered EM only, percentage of held out data
% Learn.Min_Likelihood_Change = 1;
% Learn.Folding_Iterations = 20;
%[Pw_z,Pz_d,Pd,Li,perp,beta,Learn] = pLSA(X,[],topic,Learn);
 
%% saliency
pzw = getpzw(pwz',pz');
topicnum =1;
I = saliency(pzw,sp,row,col,R,G,B,topicnum);
figure(2);
imshow(I,[]);
topicnum2 =2;
I2 = saliency(pzw,sp,row,col,R,G,B,topicnum2);
figure(3);
imshow(I2,[]);
