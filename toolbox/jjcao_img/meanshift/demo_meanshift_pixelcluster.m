 clear all 
 close all
 clc

 x = double( imread('test.png') );
 y = meanShiftPixCluster(x,20,32);

 sample = zeros(size(x,1),size(x,2));
 sample(1:3:end,1:3:end) = 1;

 R = x(:,:,1); Rx = R(sample==1); Rn = randn( numel(Rx),1 )/3;
 G = x(:,:,2); Gx = G(sample==1); Gn = randn( numel(Rx),1 )/3;
 B = x(:,:,3); Bx = B(sample==1); Bn = randn( numel(Rx),1 )/3;
 figure, 
 subplot(221), imshow(uint8(x)), axis image; title('input image')
 subplot(223), imshow(uint8(y)), axis image; title('output image')
 subplot(222)
 scatter3( Rx(:)-Rn, Gx(:)-Gn, Bx(:)-Bn, 3, [ Rx(:), Gx(:), Bx(:) ]/255 );
 title('Pixel Distribution Before Meanshift')
 xlim([0,255]),ylim([0,255]),zlim([0,255]);axis square

 R = y(:,:,1); Ry = R(sample==1); Rn = randn( numel(Rx),1 )/3;
 G = y(:,:,2); Gy = G(sample==1); Gn = randn( numel(Rx),1 )/3;
 B = y(:,:,3); By = B(sample==1); Bn = randn( numel(Rx),1 )/3;
 subplot(224)
 scatter3( Ry(:)-Rn, Gy(:)-Gn, By(:)-Bn, 3, [ Rx(:), Gx(:), By(:) ]/255 );
 title('Pixel Distribution After Meanshift')
 xlim([0,255]),ylim([0,255]),zlim([0,255]);axis square