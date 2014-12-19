% smacof    Performs least-squares multidimensional using SMACOF algorithm.
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
function [X,hist] = smacof(D,X0,iter,verbose,xhistory,rtol,atol)

% check input correctness
if nargin < 6,
    error('Incorrect number of arguments, exiting.')     
end


if size(D,1) ~= size(D,2),
    error('Matrix D must be square, exiting.')     
end

if size(D,1) ~= size(X0,1),
    error('X0 and D dimensions mismatch, exiting. X0 must be a size(D,1)*dim matrix.')     
end

if rtol < 0
    error('rtol must be non-negative, exiting.')     
end


% check input and output flags
if strcmp(lower(verbose),'iter'),
    VERBOSE = 1;
else
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


% initialize
iii = 1;
Z = X0;
X = X0;
D_ = calc_D (X);

% initialize history
if HISTORY | VERBOSE,
    %hist.time = zeros(iter);
    %hist.s = zeros(iter);
    hist.s(1) = calc_stress(X0,D);
end

if XHISTORY
   hist.X{1} = X0; 
end

if HISTORY & VERBOSE,
    fprintf(1,'iter         stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n', hist.s(1)) 
end


while (iii <= iter),
    t = cputime;       

    B_ = calc_B(D_,D);
    X = B_*Z/size(D,1);
    D_ = calc_D (X);
    
    S = calc_S (D,D_);
    Z = X;

    % add history
    if HISTORY | VERBOSE,    
        hist.time(iii) = cputime-t;
        hist.s(iii) = calc_stress(X,D);
    end
    
    if XHISTORY
        hist.X{iii} = X; 
    end

    if HISTORY & VERBOSE,
        fprintf(1,'%4d   %12.3g   %10.3g\n', iii,hist.s(iii),hist.time(iii)) 
    end
    
    
    % check stopping conditions
    if S < atol,
        fprintf(1,'atol=%g reached, exiting\n',atol)
        return 
    end

    if (iii > 1) & (HISTORY | VERBOSE),
        if (hist.s(iii-1)/hist.s(iii)-1) < rtol,
            fprintf(1,'rtol=%g reached, exiting\n',rtol)
            return 
        end
    end

    iii = iii+1;
        
end    



% SERVICE FUNCTIONS

% compute the stress 
function [S] = calc_stress (X,D)
D_ = calc_D (X);
S  = calc_S (D,D_);
return

function [S] = calc_S (D,D_)
d = triu((D - D_).^2,1);
S = sum(d(:));
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
return;

