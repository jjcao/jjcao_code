clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'data'])

img = imread('lena.jpg');
imshow(img)

figure;
im = rgb2gray(img);
imshow(im);

level = graythresh(im)
binData = im2bw(im, level);
imshow(binData);

figure;
[B,L] = bwboundaries(binData);
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end