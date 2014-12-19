function lAlphaBeta= RGB2lAlphaBeta(im)
%
% reference: cga01_Color transfer between images
% 
% jjcao @ 2014
%
[m, n, c] = size(im);
if c < 3
    error('input image should be an RGB image!');
end

rgb=[reshape(im(:,:,1),1,m*n);reshape(im(:,:,2),1,m*n);reshape(im(:,:,3),1,m*n)];

b1 = [.3811 .5783 .0402; .1967 .7244 .0782; .0241 .1288 .8444];
LMS = b1 * rgb;
LMS = log10(1+LMS);
c1 = [1/sqrt(3) 0 0;0 1/sqrt(6) 0;0 0 1/sqrt(2)];
d1 = [1 1 1;1 1 -2;1 -1 0];
lAlphaBeta = c1*d1*LMS;

% l = lAlphaBeta(1,:);
% alpha = lAlphaBeta(2,:);
% beta = lAlphaBeta(3,:);
