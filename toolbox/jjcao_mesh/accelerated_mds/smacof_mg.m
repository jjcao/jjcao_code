% smacof_mg    Performs least-squares multidimensional scaling using SMACOF 
%               algorithm with multigrid acceleration (V-cycle).
%
% Usage:
%
%  called from mds
%
% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
% Version: 1.0
%
% (C) Copyright Michael Bronstein, 2008
% All rights reserved.
%
% License:
%
% ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCES. 
% ANY COMMERCIAL USE PROHIBITED. PLEASE CONTACT THE AUTHORS FOR 
% LICENSING TERMS. PROTECTED BY INTERNATIONAL INTELLECTUAL PROPERTY 
% LAWS. PATENTS PENDING.
function [X,hist] = smacof_mg(D,X,cycles,iter,verbose,xhistory,rtol,atol,IND,UPMTX,DOWNMTX,lambda)
if nargin < 12,
    error('Incorrect number of arguments, exiting.')     
end


if size(D,1) ~= size(D,2),
    error('Matrix D must be square, exiting.')     
end

if size(D,1) ~= size(X,1),
    error('X0 and D dimensions mismatch, exiting. X0 must be a size(D,1)*dim matrix.')     
end

if rtol < 0
    error('rtol must be non-negative, exiting.')     
end


% check input and output flags
switch lower(verbose),
    case 'cycle',
        VERBOSE = 1;
    case 'iter',
        VERBOSE = 2;        
    otherwise
        VERBOSE = 0;
end

if nargout == 2,
    HISTORY = 1;
else
    HISTORY = 0;
end

if strcmp(lower(xhistory),'on'),
    XHISTORY = 1;
else
    XHISTORY = 0;
end

% construct data hierarchy
for k = 1:length(IND),
   D_{k} = D(IND{k},IND{k}); 
end

% construct number of iterations
for k = 1:length(IND),
   NITER(k) = iter;
   NITER_(k) = iter;
end


% initialize
[S,dF] = mg_lsmds(X,D_{1},zeros(size(X)), lambda);

% initialize history
if HISTORY,
    hist.s(1) = S;
end

if XHISTORY
   hist.X{1} = X; 
end


if VERBOSE,
    fprintf(1,'cycle        stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n',S) 
end



for k = 1:cycles
    tt = cputime;

    T = zeros(size(X));
    [X] = vcycle_rec(X,lambda, D_,T,IND,UPMTX,DOWNMTX, NITER,NITER_);

    TT = cputime - tt;
    if XHISTORY
        hist.X{k} = X;
    end
    
    % history
    if HISTORY
        [S,dF] = mg_lsmds(X,D_{1},zeros(size(X)), 0);
        hist.time(k) = TT;
        hist.s(k) = S;
    end

    if VERBOSE & HISTORY,
        fprintf(1,'%4d   %12.3g   %10.3g\n',k,S,TT);
    elseif VERBOSE & ~HISTORY
        [S,dF] = mg_lsmds(X,D_{1},zeros(size(X)), 0);
        fprintf(1,'%4d   %12.3g   %10.3g\n',k,S,TT);
    end

    if S < atol,
        fprintf(1,'stress reached %g, exiting\n', atol)  
        return
    end

    if k>1 & hist.s(k-1)>hist.s(k) & hist.s(k-1)/hist.s(k)-1 < rtol,
        fprintf(1,'stress changes by less than %g, exiting\n', rtol)  
        return
    end

end


return



% basic MG recursion
function [X] = vcycle_rec(X,lambda, D_,T,IND,UPMTX,DOWNMTX, NITER,NITER_)

D = D_{1};
N = size(D,1);
m = size(X,2);

X = mg_lsmds_relax(X,D,T, lambda, NITER(1)); %% (N^2(3+25) + N*m)*NITER
[dummyv1,G,dummyv2] = mg_lsmds(X,D, T, lambda); %(N^2(3+25) + N*m)

XX = DOWNMTX{1}*X;
DD = D_{2};  %UPMTX{1}'*D*UPMTX{1};    %D(IND{2},IND{2});
TT = DOWNMTX{1}*T;    %T(IND{2},:);

NN = size(DD,1);

[dummyv1,GG,dummyv2] = mg_lsmds(XX,DD, TT, lambda);

RR = GG - DOWNMTX{1}*G;
XX_ = XX;

if length(IND)<=2, %last resolution
    XX = mg_lsmds_relax(XX,DD,RR, lambda, NITER(2));

else
    for r = 1:length(IND)-1,
        INDD{r} = IND{r+1};
    end
    for r = 1:length(IND)-2,
        UPMTXX{r} = UPMTX{r+1};
    end
    for r = 1:length(IND)-2,
        DOWNMTXX{r} = DOWNMTX{r+1};
    end    
    for r = 1:length(IND)-1,
        DD_{r} = D_{r+1};
    end    
    
    [XX] = vcycle_rec(XX,lambda, DD_,RR,INDD, UPMTXX, DOWNMTXX, NITER(2:end),NITER_(2:end));
end

E = UPMTX{1}*(XX - XX_);
X = X + E;
X = mg_lsmds_relax(X,D,T, lambda, NITER_(1));




function [f,dF,dummyvar] = mg_lsmds(X,D,T,lambda)

N = size(D,1);
m = length(X(:))/N;

D_ = calc_D(X); 
B  = calc_B(D_,D); 

if nargout == 3, % N^2(3+25) + N*m
    dF = 2*( (N*X - repmat(sum(X,1),[N,1])) - B*X) - T + ...
         2*lambda*repmat(sum(X,1),[N,1]); %penalty
    f = NaN;
    dummyvar = 1;
    
elseif nargout <=2 
    d = triu((D_ - D).^2,1);
    f = sum(d(:)) - sum(sum(X.*T)) + ...
        lambda*sum(sum(X,1).^2); % penalty

    if nargout == 2,
        dF = 2*( (N*X - repmat(sum(X,1),[N,1])) - B*X) - T + ...
             2*lambda*repmat(sum(X,1),[N,1]); %penalty
    end
    
end


function [X] = mg_lsmds_relax(X,D, T, lambda, n)
N = size(D,1);
m = length(X(:))/N;

for k = 1:n,
    [dummy1,dFt,dummy2] = mg_lsmds(X,D, T, lambda); % N^2(3+25) + N*m
    X = X - 0.5*dFt/N;
end

return



function [D] = calc_D (X)
D = zeros(size(X,1));
for k=1:size(X,1),
    xk = repmat(X(k,:),size(X,1),1);
    D(:,k) = sqrt(sum((X - xk).^2, 2));
end    
return;


function [B] = calc_B (D_,D)
B = zeros(size(D));
i = find(D_(:) ~= 0);
B(i) = - D(i)./D_(i);
B = B - diag(diag(B));
d = sum(B);
B = B - diag(d);
return