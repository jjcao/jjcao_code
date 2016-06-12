function [ f , param ]  =  shiftable_jointBF(f0, f1, sigmas, ...
    sigmar, w, tol)
%
%   joint (or cross) bilateral filtering
%
%   f0   : MxNxB image to filter
%   f1   : MxN image used to for pixel similarity
%   f     : filtered image of size MxNxB
%   sigmas     : width of spatial Gaussian
%   sigmar     : width of range Gaussian
%   [-w, w]^2  : domain of spatial Gaussian
%   tol        : truncation error
%
%   Modified by: Derek Hoiem
%   Date: May 18, 2012
%   Modifications: Allow inImg to be multiband;
%   allow image for computing pixel similarity to be
%  different than filtered image.
%
% create spatial filter
filt     = fspecial('gaussian', [w w], sigmas);
% set range interval and the order of raised cosine
T  =  maxFilter(f1, w);
N  =  ceil( 0.405 * (T / sigmar)^2 );
gamma    =  1 / (sqrt(N) * sigmar);
twoN     =  2^N;
% compute truncation
if tol == 0
    M = 0;
else
    if sigmar > 40
        M = 0;
    elseif sigmar > 10
        sumCoeffs = 0;
        for k = 0 : round(N/2)
            sumCoeffs = sumCoeffs + nchoosek(N,k)/twoN;
            if sumCoeffs > tol/2
                M = k;
                break;
            end
        end
    else
        M = ceil( 0.5 * ( N - sqrt(4*N*log(2/tol)) ) );
    end
end
% main filter
warning off; %#ok<WNOFF>
[m, n, b]  =  size(f0);
fnum      =  zeros(m, n, b);
fdenom  =  zeros(m, n, b);
f   =  zeros(m, n, b);
if N < 50
    for k = M : N - M
        bk = nchoosek(N,k) / twoN;
        temp1  = cos( (2*k - N) * gamma * f1);
        temp2  = sin( (2*k - N) * gamma * f1);
        if size(f0, 3) > 1
            temp1 = repmat(temp1, [1 1 size(f0, 3)]);
            temp2 = repmat(temp2, [1 1 size(f0, 3)]);
        end
        phi1 =  imfilter(f0 .* temp1, filt);
        phi2 =  imfilter(f0 .* temp2, filt);
        phi3 =  imfilter(temp1, filt);
        phi4 =  imfilter(temp2, filt);
        fnum = fnum + bk * ( temp1 .* phi1 +  temp2 .* phi2 );
        fdenom = fdenom + bk * ( temp1 .* phi3 +  temp2 .* phi4 );
    end
else
    for k = M : N - M
        bk = exp(logfactorial(N) - logfactorial(k) ...
            - logfactorial(N-k) - N*log(2));
        temp1  = cos( (2*k - N) * gamma * f1);
        temp2  = sin( (2*k - N) * gamma * f1);
        if size(f0, 3) > 1
            temp1 = repmat(temp1, [1 1 size(f0, 3)]);
            temp2 = repmat(temp2, [1 1 size(f0, 3)]);
        end
        phi1 =  imfilter(f0 .* temp1, filt);
        phi2 =  imfilter(f0 .* temp2, filt);
        phi3 =  imfilter(temp1, filt);
        phi4 =  imfilter(temp2, filt);
        fnum = fnum + bk * ( temp1 .* phi1 +  temp2 .* phi2 );
        fdenom = fdenom + bk * ( temp1 .* phi3 +  temp2 .* phi4 );
    end
end
idx1 = find( fdenom < 1e-3);
idx2 = find( fdenom > 1e-3);
f( idx1 ) = f0( idx1 );
f( idx2 ) = fnum( idx2 ) ./ fdenom (idx2 );
% save parameters
param.T  = T;
param.N  = N;
param.M  = M;
end

