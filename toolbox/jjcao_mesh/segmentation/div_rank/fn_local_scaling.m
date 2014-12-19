function A_LS = fn_local_scaling(D, LOC_K)
%
%	Locally scaled affinity matrix based on 
%   Zelnik-Manor and Perona, Self-Tuning Spectral Clustering, NIPS 2004
%
% -------
% INPUTS:
%   D           : Distance matrix
%   LOC_K       : K-th neighbor for local scaling computation 
%
% --------
% OUTPUTS:
%
%   A_LS        : Scaled affinity matrix
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%


if nargin<2, 
    LOC_K = 7 ;
end

% sort D
sorted_D = sort(D,2,'ascend') ;

% local scales
% +1 is necessary because the most closest one is itself. 
loc_scale = sqrt(sorted_D(:,LOC_K+1)) ;

% pairwise local scales
loc_scale = loc_scale*loc_scale' ;

% If local scale is less than 0.004, just set 0.004 to it. 
% I just follow the author's implementation.
% loc_scale(loc_scale<0.004) = 0.004 ; 

% compute Gaussian affinity
A_LS = exp(- D ./ loc_scale) ; 
