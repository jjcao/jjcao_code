function xy = isomap_for_mesh(D,ndims,options)

% isomap - computes Isomap embedding using the original algorithm of 
% sci00_A global geometric framework for nonlinear dimensionality reduction
% or the adancve algorithm using landmarks for speedup: nips02_Global versus local methods in Nonlinear Dimensionality Reduction
%
%   xy = isomap(X,ndims,options); 
%       or 
%   xy = isomap(D,ndims,options); 
%
%   'X' is a D x N matrix (D = dimensionality, N = #points)
%   'D' is either a local distance matrix (with Inf when no connection),
%       or a global one (i.e. it contains already computed geodesic
%       distances between pair of points).
%   OPTIONAL:
%   'ndims' is the number of output dimensions (default =2).
%
%   options.method  - algorithm [classical_mds | rre | smacof | mg] (default: classical_mds)
%
%   Modified by Junjie Cao from the original code : Isomap code -- (c) 1998-2000 Josh Tenenbaum
%   & Gabriel Peyre's code
%
options.null = 0;

if nargin<2
    ndims = 2;
end
method = getoptions(options, 'method', 'classical_mds');
landmarks = getoptions(options, 'landmarks', []);
if isempty(landmarks)
    use_landmarks = 0;
else
    use_landmarks = 1;
    landmarks = sort(landmarks);
    landmarks = unique(landmarks);    
    nbr_landmarks = length(landmarks);
end
N = size(D,1);

if use_landmarks
    Dfull = D;
    D = D(landmarks,landmarks);
    Nfull = N;
    N = nbr_landmarks;
end
    
switch method
    case 'classical_mds'
        %%%%% Construct low-dimensional embeddings (Classical MDS) %%%%%
        opt.disp = 0; opt.isreal = 1; opt.issym = 1; 
        M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
        [xy, val] = eigs(M, ndims, 'LR', opt); 

        for i=1:ndims
            xy(:,i) = xy(:,i)*sqrt(val(i,i));
        end
    otherwise
        options.dim = ndims;
        options.method = method;
        xy = mds(D,options);
end


if use_landmarks
    % interpolation on the full set of points
    % x = 1/2 * (L^T) * ( delta_n-delta_x )
    xy1 = zeros(Nfull,ndims);
    % transpose of embedding
    LT = xy'; 
    for i=1:ndims
        LT(i,:) = LT(i,:) / val(i,i);
    end
    deltan = mean(D,2);
    for x=1:Nfull
        deltax = Dfull(landmarks,x).^2;
        xy1(x,:) = 1/2 * ( LT * ( deltan-deltax ) )';
    end
    xy = xy1;
end

