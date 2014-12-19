function [potentials] = get_potential_mult_lineq_lu(cond_mat, pos_set, neg_pole)
%
%	Compute Potentials for multiple sets
%   Solution by direct linear method (Incomplete LU decomposition)
%   See the Page 21 of 'Random Walks and electric networks'
%       by  P.G. Doyle and J.L. Snell
%
% -------
% INPUTS:
%
%	cond_mat        : conductance network
%   pos_set         : 1 if nodes are attached to sources. Otherwise, 0
%                     [N x M] where N: # of nodes, M: # of sets
%   neg_pole        : nodes are attached to ground ([L x 1])
%                     where L is the number of negative nodes.
%                     Unlike pos_set, neg_pole fixed in all sets.
%
% --------
% OUTPUTS:           
%   potentials      : potential  ([N x M])
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%


% if neg_pole is not assigned, report an error.
if isempty(neg_pole)
    error('Negative pole should be given.');
end

% N: number of vertices
% M: number of sets
[N M] = size(pos_set) ;

% output potential 
potentials = zeros(N,M) ;

% Laplace
lap_mat = full(diag(sum(cond_mat,2)) - cond_mat) ;


% Parameters for 'ilu'
ilu_setup=[];
ilu_setup.type = 'ilutp';

% for each pos set, compute the potential. 
for i=1:M,
    potentials(:,i) = get_potential_lineq_lu([], lap_mat, ...
        find(pos_set(:,i)), neg_pole, ilu_setup);
end

end
