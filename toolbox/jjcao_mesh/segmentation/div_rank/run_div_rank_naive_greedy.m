function [rank_set rank_set_val] = run_div_rank_naive_greedy(...
                        adj_mat, K, lin_mode, Z, EvalList, massPnt)
%
%   Run diversity ranking using naive greedy method.
%
% -------
% INPUTS:
%
%	adj_mat     : Similarity network ([N x N]) 
%                 where N is the number of vertices 
%   K           : number of ranked sets
%   lin_mode    : methods to solve linear eq Ax=b.
%                 0: incremental inverse. 1: LU decomposition.
%                 Others: naive backslash. 
%   Z           : ground conductance ([N x 1]) 
%   massPnt     : mass (weights) of each point ([N x 1])
%                 (default: ones(N,1)).
%   EvalList    : The vertices that are evaluated to search for max gain.
%                 As default, we evaluate all vertices. 
%
% --------
% OUTPUTS:
%
%   rank_set        : selected rank_sets ([1 x K]) 
%   rank_set_val    : potentials ([1 x K]) 
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%

% number of vertex
num_vtx = size(adj_mat,1) ;

% default parameters
if nargin<3,    lin_mode = 0 ;              end
if nargin<4,    Z = 0.001*ones(num_vtx,1) ; end
if nargin<5,    EvalList = 1:num_vtx ;      end
if nargin<6,    massPnt = ones(num_vtx,1) ; end

% number of evaluation points.
nEvalList = length(EvalList) ;

% Add dummy ground and normalize it.
adj_mat = [[adj_mat Z]; [ones(1,num_vtx) 0]] ;
adj_mat = fn_normout(adj_mat) ;
neg_pole = num_vtx + 1 ;
% all indices except negative pole.
ind_ord = 1:num_vtx ;

% output
rank_set = [] ; rank_set_val = [] ; 

% Greedily select K-1 more nodes 
% Picking out the node causing max potential increase one by one.
% Once picked out, the item turns into voltage sources.
for i=1:K,
    tic
    % Get the sets to be evaluated. 
    pos_set = sparse(EvalList, 1:nEvalList, ones(1,nEvalList), num_vtx+1, nEvalList) ;
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

    % potential of SP = computed potential * number of pixels of SP.
    potential(ind_ord,:) = potential(ind_ord,:) .* repmat(massPnt, [1 nEvalList]) ;

    %{
    % (DEBUG)
    [potential1] = get_potential_mult_lineq_naive(adj_mat, pos_set, neg_set) ;
    max(max(abs(potential-potential1)))
    %}

    clear pos_set neg_set ;

    % save the vertex to maximuze the marginal temperature gain.
    [max_v max_i] = max(sum(potential)) ;
    max_i = EvalList(max_i) ;
    rank_set = [rank_set max_i] ;
    rank_set_val = [rank_set_val max_v] ;

    % display the result
    t2 = toc ;
    str_rank_set = []; 
    for j=1:length(rank_set),
        str_rank_set = [str_rank_set num2str(rank_set(j)) ' '] ;
    end
    display(['[ ' str_rank_set '] are chosen as sources. ' ...
        num2str(t2) ' sec.']) ;
end

%{
% I implemented three different ways to compute potnetials
% (1) power method, (2) Using '\' (3) LU decomposition
potential2 = get_potential_power(adj_mat, rank_set, num_vtx+1) ;
potential = get_potential_lineq(adj_mat, [], rank_set, num_vtx+1) ;
potential = get_potential_lineq_lu(adj_mat, [], rank_set, num_vtx+1) ;
%}

end
