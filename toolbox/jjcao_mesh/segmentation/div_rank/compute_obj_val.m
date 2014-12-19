function [obj_val rank_set] = compute_obj_val(adj_mat, K, Z, lin_mode)
%
%   Compute Object Function Values (ie. potentials)
%
% -------
% INPUTS:
%
%	adj_mat     : Similarity network ([N x N]) 
%                 where N is the number of vertices 
%   K           : number of ranked sets
%   Z           : ground conductance
%   lin_mode    : methods to solve linear eq Ax=b.
%                 0: incremental inverse. 1: LU decomposition.
%                 Others: naive backslash. 
%
% --------
% OUTPUTS:
%
%   obj_val         : object function values (ie. potentials) ([N x K]) 
%   rank_set        : selected sources ([1 x K]) 
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%

% number of vertex
num_vtx = size(adj_mat,1) ;

% Add dummy ground and normalize it.
adj_mat = [[adj_mat Z*ones(num_vtx,1)]; [ones(1,num_vtx) 0]] ;
adj_mat = fn_normout(adj_mat) ;
neg_pole = num_vtx + 1 ;

% output
rank_set = [] ;
obj_val = zeros(num_vtx, K) ;

% Greedily select K-1 more nodes 
% Picking out the node causing max potential increase one by one.
% Once picked out, the item turns into voltage sources.
for i=1:K,
    % Get the sets to be evaluated. 
    pos_set = sparse(1:num_vtx, 1:num_vtx, ones(1,num_vtx), num_vtx+1, num_vtx) ;
    pos_set(rank_set,:) = 1 ;
    neg_set = neg_pole ;
    % (1) incremental inverse matrix
    if lin_mode==0,
        [potential] = get_potential_mult_lineq_inc_inv(adj_mat, pos_set, neg_set) ;
    % (2) LU decomposition
    elseif lin_mode == 1,
        [potential] = get_potential_mult_lineq_lu(adj_mat, pos_set, neg_set) ;
    % (3) Naive method using backslash.
    else
        [potential] = get_potential_mult_lineq_naive(adj_mat, pos_set, neg_set) ;
    end

    % object function values
    obj_val(:,i) = sum(potential)' ;

    % save the vertex to maximuze the marginal temperature gain.
    [max_v max_i] = max(sum(potential)) ;
    rank_set = [rank_set; max_i] ;

    %{
    % display the result
    t2 = toc ;
    str_rank_set = []; 
    for j=1:length(rank_set),
        str_rank_set = [str_rank_set num2str(rank_set(j)) ' '] ;
    end
    display(['[ ' str_rank_set '] are chosen as sources. ' ...
        num2str(t2) ' sec.']) ;
    %}
end

end
