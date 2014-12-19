%
% basic implementation of the paper, not including the two potentional
% improvements.
%
% Improvement 1: will be done later
% Improvement 2: no plan for it since it may be not improtant in practice?
%
% 11.jpg and 12b.jpg make the algorithm fail, why? 
% anwser: LMS = log10(LMS); used in RGB2lAlphaBeta leads to the problem, so use
% LMS = log10(1+LMS);. In the same time, we need change the inverse transformation
% lAlphaBeta2RGB (LMS = 10.^LMS-1;).
%
% warning: move computation to GPU not improve the speed, some times the
% speed is very slow.
%
% reference: cga01_Color transfer between images
%
% Copyright (c) 2014 Junjie Cao

close all;clear all; % clc;
addpath(genpath('../../'));
USE_GPU = false;

src=imread('11.jpg');
if USE_GPU
    src = gpuArray(src);
end
src = im2double(src);

tgt=im2double(imread('12b.jpg'));
tgt = im2double(tgt);
if USE_GPU
    tgt = gpuArray(tgt);
end
figure,subplot(131);imshow(src);title('Source');
subplot(132);imshow(tgt);title('Target');

%%
tic
% for i = 1:1000
lAlphaBetaT = RGB2lAlphaBeta(tgt);
lAlphaBetaS = RGB2lAlphaBeta(src);

lAlphaBetaN = channels_transfer(lAlphaBetaT, lAlphaBetaS);
im=lAlphaBeta2RGB(lAlphaBetaN, size(src));
% end
toc

%%
subplot(133);imshow(im);title('Result');
