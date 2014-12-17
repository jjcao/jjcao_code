% http://www.klab.caltech.edu/~xhou/
% cvpr07_Saliency Detection A Spectral Residual Approach
%
% changed by jjcao @ 2014
%

clear;clc;close all;
addpath(genpath('../../../'));
%% Read image from file 
inImg = im2double(rgb2gray(imread('curve.jpg')));
inImg = imresize(inImg, 64/size(inImg, 2));
figure(1);
subplot(1,3,1); imshow(inImg);
%% Spectral Residual
myFFT = fft2(inImg); 
myAmplitude = abs(myFFT);
myLogAmplitude = log(myAmplitude);
smoothedLogAmplitude = imfilter(myLogAmplitude, fspecial('average', 3), 'replicate');
mySpectralResidual = myLogAmplitude - smoothedLogAmplitude; 

figure(2);
subplot(1,3,1);imagesc(myAmplitude);colorbar;
subplot(1,3,2);imagesc(myLogAmplitude);colorbar;
subplot(1,3,3);imagesc(mySpectralResidual);colorbar;

myPhase = angle(myFFT);
saliencyMap = abs(ifft2(exp(mySpectralResidual + i*myPhase))).^2;
rImg = exp(smoothedLogAmplitude + i*myPhase);
figure(1); subplot(1,3,2); imshow(rImg);
%% After Effect
saliencyMap = mat2gray(imfilter(saliencyMap, fspecial('gaussian', [10, 10], 2.5)));
figure(1); subplot(1,3,3); imshow(saliencyMap);
