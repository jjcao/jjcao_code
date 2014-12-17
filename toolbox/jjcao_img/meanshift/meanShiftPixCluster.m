function [y, MS] = meanShiftPixCluster(x,hs,hr,th,plotOn)
% FUNCTION: meanShiftPixCluster implements the classic mean shift pixel
% clustering algorithm introduced in Cmaniciu etal.'s PAMI paper 
% "Mean shift: a robust apporach toward feature space analysis", 2002.
% -------------------------------------------------------------------------
% Input: 
%       x = an input image (either gray or rgb, please expect long time processing if image size is large)
%      hs = the bandwidth of spatial kernel (see Eq.(35) in the cited paper)
%      hr = the bandwidth of feature kernel (see Eq.(35) in the cited paper)
%      th = the threshod of the convergence criterion (default = .25)
%  plotOn = switch on/off the image display of intermediate results (default = 1)
%
% Output:
%       y = the output pixel clustered image
%      MS = the output of averaged mean shift
% -------------------------------------------------------------------------
% Demo (please cut and paste):
%      clear all 
%      close all
%      clc
% 
%      x = double( imread('test.png') );
%      y = meanShiftPixCluster(x,20,32);
%      
%      sample = zeros(size(x,1),size(x,2));
%      sample(1:3:end,1:3:end) = 1;
% 
%      R = x(:,:,1); Rx = R(sample==1); Rn = randn( numel(Rx),1 )/3;
%      G = x(:,:,2); Gx = G(sample==1); Gn = randn( numel(Rx),1 )/3;
%      B = x(:,:,3); Bx = B(sample==1); Bn = randn( numel(Rx),1 )/3;
%      figure, 
%      subplot(221), imshow(uint8(x)), axis image; title('input image')
%      subplot(223), imshow(uint8(y)), axis image; title('output image')
%      subplot(222)
%      scatter3( Rx(:)-Rn, Gx(:)-Gn, Bx(:)-Bn, 3, [ Rx(:), Gx(:), Bx(:) ]/255 );
%      title('Pixel Distribution Before Meanshift')
%      xlim([0,255]),ylim([0,255]),zlim([0,255]);axis square
% 
%      R = y(:,:,1); Ry = R(sample==1); Rn = randn( numel(Rx),1 )/3;
%      G = y(:,:,2); Gy = G(sample==1); Gn = randn( numel(Rx),1 )/3;
%      B = y(:,:,3); By = B(sample==1); Bn = randn( numel(Rx),1 )/3;
%      subplot(224)
%      scatter3( Ry(:)-Rn, Gy(:)-Gn, By(:)-Bn, 3, [ Rx(:), Gx(:), By(:) ]/255 );
%      title('Pixel Distribution After Meanshift')
%      xlim([0,255]),ylim([0,255]),zlim([0,255]);axis square
% -------------------------------------------------------------------------
% Additional Info:
% This is a toy code of mean shift pixel clustering, but it did cover the
% core idea od mean shift. One may extend the current code for
% multiresolution analysis, or replace new kernel functions instead of
% gaussian kernels used in this implementations. Details of parameters and
% their influence on clustering results and the convergence proof of the
% mean shift algorithm please refer the cited paper. 
%
% Please play with this code with other images and parameters. Contact me
% if you find any bugs.
%
% by Yue Wu
% rex.yue.wu@gmail.com
% 11/30/2013
% 
% ------------------------------------------------------------------------- 

%% Argument Check
if nargin<3
    error('please type help for function syntax')
elseif nargin == 3
    th = 1/100; plotOn = 1;
elseif nargin == 4
    if th<0 || th >255
        error('threshold should be in [0,255]')
    else
        plotOn = 1;
    end
elseif nargin == 5
    if sum(ismember(plotOn,[0,1])) == 0
        error('plotOn option has to be 0 or 1')
    end
elseif nargin>5
    error('too many input arguments')
end
%% initialization
x = double(x);
[height,width,depth] = size(x);
y = x;
done = 0;
iter = 0;
if plotOn 
    figure(randi(1000)+1000); 
end
% padding image to deal with pixels on borders
xPad = padarray(x,[height,width,0],'symmetric');
% build up look up table to boost computation speed
weight_map = exp( -(0:255^2)/hr^2 );
MS = [];
%% main loop
while ~done
    weightAccum = 0;
    yAccum = 0;
    % only 99.75% area (3sigma) of the entire non-zero Gaussian kernel is considered
    for i = -hs:hs
        for j = -hs:hs
            if ( i~=0 || j~=0 )
                % spatial kernel weight 
                spatialKernel = 1;
                % uncomment the following line to active Gausian kernel
                %spatialKernel = exp(-(i^2+j^2)/(hs/3)^2/2);
                xThis =  xPad(height+i:2*height+i-1, width+j:2*width+j-1, 1:depth);
                xDiffSq = (y-xThis).^2;
                % feature kernel weight
                intensityKernel = repmat( prod( reshape( weight_map( xDiffSq+1 ), height, width, depth) , 3 ), [1,1, depth]);
                % mixed kernel weight
                weightThis = spatialKernel.*intensityKernel;
                % update accumulated weights
                weightAccum = weightAccum+ weightThis;
                % update accumulated estimated ys from xs 
                yAccum = yAccum+xThis.*weightThis;
            end
        end
    end
    % normalized y (see Eq.(20) in the cited paper)
    yThis = yAccum./(weightAccum+eps);
    % convergence criterion
    yMS = mean(abs(round(yThis(:))-round(y(:))));
    y = round(yThis);    
    MS(iter+1) = yMS;
    %xPad = padarray(y,[height,width,0],'symmetric');
    if plotOn
        subplot(121), imshow(uint8(y)),axis image, title(['iteration times = ' num2str(iter) '; averaged mean-shift = ' num2str(yMS)]);
        subplot(122), plot(0:iter, MS ), xlabel('iteration #'), ylabel('averaged mean shift');axis square
        drawnow
    end
    if yMS <= th % exit if converge
        done = 1;
    else % otherwise update estimated y and repeat mean shift
        iter = iter+1;
    end
end
