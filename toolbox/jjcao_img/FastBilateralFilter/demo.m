% Fast Bilateral Filter Using Raised Cosines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  f0             :  input grayscale image
%  f               : bilateral filter output
%  sigmas      : width of spatial Gaussian
%  sigmar      : width of range Gaussian
%  [-w, w]^2  : domain of spatial Gaussian
%  tol             : truncation error (kernel)
%
%  Author:    Kunal N. Chaudhury.
%  Date:        March 1, 2012.
%  Modified:  July 13, 2015.
%
%  References:
%  [1] K.N. Chaudhury, D. Sage, and M. Unser, "Fast O(1) bilateral
%  filtering using trigonometric range kernels," IEEE Trans. Image Proc.,
%  vol. 20, no. 11, 2011.
%
% [2] K.N. Chaudhury, "Acceleration of the shiftable O(1) algorithm for
% bilateral filtering and non-local means,"  IEEE Transactions on Image Proc.,
% vol. 22, no. 4, 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% load test image
clc, clear, close all force;
f0  =  double( imread('./images/1.png') );
if size(f0,3) 
    f0 = f0(:,:,1);
end
[m, n] = size(f0);
% filter parameters
sigmas = 2;
sigmar = 25;
tol   = 0.01;
% window for spatial Gaussian
w  = round(6*sigmas);
if (mod(w,2) == 0)
    w  = w+1;
end
% call bilateral filter
tic;
[f, param] =  shiftableBF(f0, sigmas, sigmar, w, tol);
toc;
% results
T  = param.T;
N  = param.N;
M  = param.M;
gamma  =  1 / (sqrt(N) * sigmar);

s  = linspace(-T, T, 200);
g  = exp( -s.^2 / (2 * sigmar *sigmar) );
gApprox  = cos(gamma * s).^N;

warning('off'); %#ok<WNOFF>
if M ==  0
    gTrunc = gApprox;
else
    if N < 50
        twoN   =  2^N;
        gTrunc = zeros( 1, length(s) );
        for k = M : N - M
            gTrunc = gTrunc + (nchoosek(N,k)/twoN) * ...
                cos((2*k - N)*gamma*s);
        end
    else
        gTrunc = zeros( 1, length(s) );
        for k = M : N - M
            % use Sterling's approximation
            factor = exp(logfactorial(N) - logfactorial(k) ...
                - logfactorial(N-k) - N*log(2));
            gTrunc = gTrunc + factor * ...
                cos((2*k - N)*gamma*s  );
        end
    end
end

figure('Units','normalized','Position',[0 0.5 1 0.5]);
plot(s, g, 'b');
hold on,
plot(s, gApprox, 'm'),
hold on,
plot(s, gTrunc, 'r');
axis('tight'), grid('on'),
legend('Gassian','Raised cosine','Truncated raised cosine','FontSize', 10);
title('Comparison of the range kernels', 'FontSize', 10),

%%
figure('Units','normalized','Position',[0 0.5 1 0.5]);
colormap gray,
subplot(1,5,1), imshow(uint8(f0)),
title('Input', 'FontSize', 10), axis('image', 'off');
subplot(1,5,2), imshow(uint8(f)),
title('Output','FontSize', 10),
axis('image', 'off');

f1 = f0-f;
f1 = f1/max(max(f1))*255;
subplot(1,5,3), 
% imagesc(uint8(f1)),
imshow(uint8(f1)),
title('Output','FontSize', 10),
axis('image', 'off');

%%%%%%%%%%%%%%%%%%%
f = imgaussfilt(f0, 5);
subplot(1,5,4), imshow(uint8(f)),
title('Output','FontSize', 10),
axis('image', 'off');

f1 = f0-f;
f1 = f1/max(max(f1))*255;
subplot(1,5,5), imshow(uint8(f1)),
title('Output','FontSize', 10),
axis('image', 'off');

