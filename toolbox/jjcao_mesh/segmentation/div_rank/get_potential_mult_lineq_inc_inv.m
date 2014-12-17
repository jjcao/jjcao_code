function [potentials] = get_potential_mult_lineq_inc_inv(cond_mat, ...
                        pos_set, neg_pole)
%
%	Compute Potentials for multiple sets
%   Solution by incremental inverse
%   Refer to http://en.wikipedia.org/wiki/Block_LU_decomposition
%
% -------
% INPUTS:
%
%	cond_mat        : conductance network
%   pos_set         : 1 if nodes are attached to sources. Otherwise, 0
%                     [N x M] where N: # of nodes, M: # of sets
%                     Each set is assumed to have the same number of positive poles
%                     except the already selected vertices 
%   neg_pole        : nodes are attached to ground ([1 x L])
%                     where L is the number of negative nodes.
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

% Positive poles shared by all sets
shared_pos_pole = find(sum(pos_set,2)==M)' ;
pos_set(shared_pos_pole,:) = 0 ;
[rr cc] = find(pos_set) ;

% check pos_set is OK.
if max(hist(cc,1:M)) > 1,
    error('Each set is assumed to have at most (num of shared poles)+1 positive poles.');
end

% output potential 
potentials = zeros(N,M) ;

% Laplace
lap_mat = full(diag(sum(cond_mat,2)) - cond_mat) ;

% indices of unlabeled points
uind = 1:N ;
uind([shared_pos_pole neg_pole]) = [];
nuind = length(uind) ;

% if the matrix is singular, use the lu decomposition method. 
if det(lap_mat(uind,uind))==0, 
    pos_set(shared_pos_pole,:) = pos_set(shared_pos_pole,:) + 1 ;
    [potentials] = get_potential_mult_lineq_lu(cond_mat, pos_set, neg_pole);
    return ;
end

% compute inverse of Laplace
inv_ulap_mat = inv(lap_mat(uind,uind)) ;

% for each pos set, 
for i=1:M,

    % the positive pole of this set. 
    cur_pos_pole = find(pos_set(:,i)) ;

    % skip if there's no new pole.
    if isempty(cur_pos_pole), 
        loc_inv = inv_ulap_mat ;
        loc_uind = uind ; 
    else
        % index of 'cur_pos_pole' in 'uind'
        puind = find(uind==cur_pos_pole) ;
        % local uind
        tmp_ind = 1:nuind ; tmp_ind(puind) = [] ;
        loc_uind = uind(tmp_ind) ;    

        % updated inverse by removing 'puind'-th row and column. 
        loc_inv = inv_ulap_mat(tmp_ind,tmp_ind) - ...
            (( inv_ulap_mat(tmp_ind,puind) / inv_ulap_mat(puind,puind) ) * ...
            inv_ulap_mat(puind,tmp_ind) ) ;     

    end

    % Find RH
    pos_pole = [shared_pos_pole cur_pos_pole] ;
    b = -sum(lap_mat(loc_uind, pos_pole),2) ;

    % compute potentials. 
    pot = loc_inv * b ;

    % compute potential
    potentials(loc_uind,i) = pot ; 
    potentials(pos_pole,i) = 1 ;
    potentials(neg_pole,i) = 0 ;
end

end
