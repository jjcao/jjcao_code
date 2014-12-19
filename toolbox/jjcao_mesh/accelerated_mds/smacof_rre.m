% smacof_rre    Performs least-squares multidimensional using SMACOF 
%               algorithm with RRE acceleration.
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
function [X,hist] = smacof_rre(D,X0,cycles,iter,verbose,xhistory,rtol,atol)

% check input correctness
if nargin < 7,
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

% initialize
iii_ = 1;
X = X0;
D_ = calc_D (X);
S = calc_stress(X0,D);

% initialize history
if HISTORY,
    hist.s(1) = S;
end

if XHISTORY
   hist.X{1} = X0; 
end


if (VERBOSE == 2),
    fprintf(1,'cycle          iter         stress   time (sec)\n') 
    fprintf(1,'INIT   ------------   %12.3g   ----------\n',S) 
elseif (VERBOSE == 1),
    fprintf(1,'cycle        stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n',S) 
end


while (iii_ <= cycles),
    t_ = cputime;       

    % internal iteration - SMACOF
    XX = zeros(prod(size(X)),iter);
    Z = X;
    for iii = 1:iter,
        t = cputime;
        B_ = calc_B(D_,D);
        X = B_*Z/size(D,1);
        XX(:,iii) = X(:);
        D_ = calc_D (X);
    
        S = calc_S (D,D_);
        Z = X;
        
        T = cputime-t;
        if HISTORY
            hist.time((iii_-1)*iter + iii + 1) = T;
            hist.s((iii_-1)*iter + iii + 1) = S;    
        end

        if XHISTORY
            hist.X{(iii_-1)*iter + iii + 1} = X; 
        end        
        
        if (VERBOSE == 2),
            fprintf(1,' ...           %4d   %12.3g   %10.3g\n',iii,S,T) 
        end
        
        % check stopping conditions
        if S < atol,
            fprintf(1,'atol=%g reached, exiting\n',atol)
            return 
        end
    end
        
    % extrapolate
    X_ = reshape(extrap(XX),size(X));
    
    % safeguard
    S_ = calc_S (D,calc_D(X_));
    if VERBOSE == 2,
        fprintf(1,'Extrap. stress: %12.3g, SMACOF stress: %12.3g\n',S_,S) 
    end
    
    if S_ <= S 
        X = X_;
        S = S_;
    elseif VERBOSE == 2,
        fprintf(1,'Safeguard: using SMACOF solution\n') 
    end
    
    
    
    
    T_ = cputime-t_;

    if (VERBOSE == 2),
        fprintf(1,'%4d   ------------   %12.3g   %10.3g\n',iii_,S,T_) 
    elseif (VERBOSE == 1),
        fprintf(1,'%4d   %12.3g   %10.3g\n',iii_,S,T_) 
    end


    % rtol stopping condition
    if (iii_ > 1) & HISTORY,
        if (hist.s((iii_-1)*iter + 1)/S-1) < rtol,
            fprintf(1,'rtol=%g reached, exiting\n',rtol)
            return 
        end
    end

    iii_ = iii_+1;
        
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

