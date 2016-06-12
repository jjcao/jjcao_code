function [ f , param ]  =  shiftableBF(f0, sigmas, sigmar, w, tol)
% create spatial filter
filt     = fspecial('gaussian', [w w], sigmas);
% set range interval and the order of raised cosine
T  =  maxFilter(f0, w);
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
        M = ceil( 0.5*( N - sqrt(4*N*log(2/tol)) ) );
    end
end
% filter
h = waitbar(0, 'Computing filter ...');
warning('off'); %#ok<WNOFF>
[m, n]   =  size(f0);
fnum  =  zeros(m, n);
fdenom  =  zeros(m, n);
f  =  zeros(m, n);
ii = sqrt(-1);
if N < 50
    for k = M : N - M
        waitbar(k / N - 2*M + 1, h);
        omegak = (2*k - N)*gamma;
        bk = nchoosek(N,k) / twoN;
        H  = exp(-ii*omegak*f0);
        G  = conj(H);
        F  = G.*f0;
        barF  = imfilter(F, filt);
        barG = imfilter(G, filt);
        fnum =  fnum + bk * H .* barF;
        fdenom  = fdenom + bk * H .* barG;
    end
    close(h);
else
    for k = M : N - M
        waitbar(k / N - 2*M + 1, h);
        omegak = (2*k - N)*gamma;
        % use Sterling's approximation
        bk = exp(logfactorial(N) - logfactorial(k) ...
            - logfactorial(N-k) - N*log(2));
        H  = exp(-ii*omegak*f0);
        G  = conj(H);
        F  = G.*f0;
        barF  = imfilter(F, filt);
        barG = imfilter(G, filt);
        fnum =  fnum + bk * H .* barF;
        fdenom  = fdenom + bk * H .* barG;
    end
    close(h);
end
% check: avoid division by zero
idx1 = find( fdenom < 1e-3);
idx2 = find( fdenom > 1e-3);
f( idx1 ) = f0( idx1 );
f( idx2 ) = real(fnum( idx2 )./fdenom (idx2 ));
% save parameters
param.T  = T;
param.N  = N;
param.M  = M;
end

