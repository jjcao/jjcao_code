function [T_i] = disp_get_potential_per_source(adj_mat, source_set, Z)
%
%   Compute temperature contributed by each source
%
% -------
% INPUTS:
%
%	adj_mat         : Similarity network ([N x N]) 
%                     where N is the number of vertices 
%   source_set      : sources ([S x 1]) where S is the number of sources
%   Z               : ground conductance
%
% --------
% OUTPUTS:
%
%   T_i             : Temperature contributed by each souce 
%                     ([N x S])
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

% number of super-pixels
numSP = size(adj_mat,1) ;

if nargin<3,
    Z = 0.001*ones(numSP,1) ;
end

% number of sources
num_source = length(source_set) ;

% Add dummy ground and normalize it.
adj_mat = [[adj_mat Z]; [ones(1,numSP) 0]] ;
adj_mat = fn_normout(adj_mat) ;
neg_pole = numSP + 1 ;

% compute temperature for each source
T_i = zeros(numSP+1, num_source) ;
for j=1:num_source,
    T_i(:,j) = get_potential_lineq(adj_mat, [], source_set(j), ...
                neg_pole) ;
end
T_i(end,:) = [] ;
