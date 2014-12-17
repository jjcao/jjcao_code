function weights = fn_weights_gauss_ls(feat, indA, indB, LOC_K)
%
%  Compute Gaussian similarity (with local scaling)
%  Local scaling: Zelnik-Manor and Perona, Self-Tuning Spectral Clustering, NIPS 2004
%
% -------
% INPUTS:
%
%   feat        : The Feature vectors ([# of SP x feat dim])
%   indA        : The first list of indices for feature computation.
%   indA        : The second list of indices for feature computation.
%   LOC_K       : K-th neighbor for local scaling computation 
%
% --------
% OUTPUTS:
%
%   weights:    An Mx2 weight list of M edges.
%
% ---------------
% Gunhee Kim
%

%Constants
EPS = 1e-5;

if ~isempty(indA) || ~isempty(indB), 
    [rA cA] = size(indA) ;
    [rB cB] = size(indB) ;
    % error check
    if rA ~= rB || cA ~= cB, 
        error('XX: indA and indB should have the same dimensions.') ;
    end
    indA = indA(:) ; 
    indB = indB(:) ;

    % compute the distances of pairs only.
    featDist = sum((feat(indA,:) - ...
        feat(indB,:)).^2,2) ;

    indA = reshape(indA, [rA cA]) ;
    indB = reshape(indB, [rB cB]) ;

    featDist = reshape(featDist, [rA cA]) ;
% if edge is empty, compute pairwise diatance
else
    featDist = fn_dist_l2_sqrt(feat, feat) ;
end

% local scaling 

% sort distance
sorted_D = sort(featDist, 2, 'ascend') ;

% local scales
% +1 is necessary because the most closest one is itself. 
loc_scale = sqrt(sorted_D(:,min(LOC_K+1,size(sorted_D,2)))) ;

% If local scale is less than 0.004, just set 0.004 to it. 
% I just follow the author's implementation.
loc_scale(loc_scale<0.001) = 0.001 ; 

% compute Gaussian affinity
if ~isempty(indA) || ~isempty(indB),
    weights = exp(- (featDist ./ (loc_scale(indA) .* loc_scale(indB)))) ; 
else
    weights = exp(- (featDist ./ (loc_scale*loc_scale'))) ; 
end

end


% normalize to [0,1]
function data = normalize(data)

    sr = size(data,1) ;
    sc = size(data,2) ;

    % normalize to [0,1]
    minDist = min(data(:)) ;
    maxDist = max(data(:)) ;
    data  = data - minDist ;
    diff = maxDist - minDist ; 

    if diff==0, 
        weights = ones(sr, sc) ;
    else
        data = data / diff ;
    end

end


