function [clusterID probMat] = run_cluster_from_rankset(adj_mat, srcInd, clID)
%
%	Clustering from Source Points (ie. ranked items).
%
% -------
% INPUTS:
%
%	adj_mat 		: Similarity network ([N x N])
%   srcInd          : Source points. ([1 x K])
%   clID            : cluster IDs of source points [1 x K] 
%
% --------
% OUTPUTS:
%
%   clusterID       : cluster id of each point ([N x 1])
%   probMat         : probability matrix ([N x K])
%
% -------------
%
% Gunhee Kim
%

% Laplacian
L = diag(sum(adj_mat,2)) - adj_mat ;

% number of points
N=length(L);
% number of clusters
K = max(clID) ;

% Px1 Vector specifying the "boundary" nodes
srcMat = zeros(K, K) ;
for i=1:K,
    s_ind = find(clID==i) ;
    srcMat(i,s_ind) = 1 ;
end

% indices that are not sources
plainInd = 1:N ;
plainInd(srcInd) = [] ;

% Find RHS
b = -L(plainInd,srcInd) * srcMat ;

% Solve linear system
x = L(plainInd,plainInd) \ b ;

% Probability
probMat = zeros(size(srcMat)) ;
probMat(srcInd,:) = srcMat ;
probMat(plainInd,:) = x ;

% cluster ID
[foo clusterID] = max(probMat,[],2) ;

% if prob = 0, assign cluster ID by K+1. 
clusterID(sum(probMat,2)<0.1) = K+1 ;
%clusterID(isnan(sum(probMat,2))) = K+1 ;

probMat(probMat<0) = 0 ;
