function im = lAlphaBeta2RGB(lAlphaBeta, imsize)
%
% reference: cga01_Color transfer between images
% 
% jjcao @ 2014
%

d2=[1 1 1;1 1 -1;1 -2 0];
c2=[1/sqrt(3) 0 0;0 1/sqrt(6) 0;0 0 1/sqrt(2)];
LMS = d2 * c2 * lAlphaBeta;
LMS = 10.^LMS-1;

b2=[4.4687   -3.5887    0.1196
   -1.2197    2.3831   -0.1626
    0.0585   -0.2611    1.2057];
rgb = b2 * LMS;

% im = zeros(imsize);
im(:,:,1) = reshape( rgb(1,:), imsize(1), imsize(2) );
im(:,:,2) = reshape( rgb(2,:), imsize(1), imsize(2) );
im(:,:,3) = reshape( rgb(3,:), imsize(1), imsize(2) );
