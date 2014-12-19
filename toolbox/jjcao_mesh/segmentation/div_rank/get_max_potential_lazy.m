function [max_vtx max_val margin_benefit] = get_max_potential_lazy( lap_mat, ...
        EvalList, margin_benefit, rank_set, neg_pole, prev_max_obj, lin_mode, massSP)
        
%
%   Lazy greedy algorithm (refer to the slide in http://submodularity.org/)
%
% -------
% INPUTS:
%
%	lap_mat         : Laplace matrix ([N x N])
%   EvalList        : The vertices that are evaluated to search for max gain.
%                     As default, we evaluate all vertices. 
%   margin_benefit  : marginal benefit ([1 x N]) 
%   prev_max_obj    : maximum objective value in the previous step.
%   rank_set        : selected items
%   neg_pole        : negative pole
%   lin_mode        : methods to solve linear eq Ax=b.
%                     0: incremental inverse. 1: LU decomposition.
%                     Others: naive backslash. 
%   massSP          : number of pixels per SP. ([N x 1])
%
% --------
% OUTPUTS:           
%   max_vtx         : next source to be selected.
%   max_val         : max objective value.
%   margin_benefit  : marginal benefit ([1 x N]) 
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%


if nargin<8, massSP=[] ; end

% indices except negative pole
ind_ord = 1:size(lap_mat,1) ;
ind_ord(neg_pole) = [] ;

% sort margin benefit
[foo eval_order] = sort(margin_benefit,'descend') ;

% do some pre-processing for (1) incremental inverse matrix
if lin_mode==0,

    shared_pos_pole = rank_set ;
    N = size(lap_mat,1) ;

    % indices of unlabeled points
    uind = 1:N ;
    uind([shared_pos_pole neg_pole]) = [];
    nuind = length(uind) ;

    % if the matrix is singular, use the lu decomposition method. 
    if det(lap_mat(uind,uind))==0, 
        lin_mode = 1 ;
    end

    % compute inverse of Laplace
    inv_ulap_mat = inv(lap_mat(uind,uind)) ;
end
% Parameters for 'ilu'
if lin_mode==1, 
    ilu_setup=[];
    ilu_setup.type = 'ilutp';
end

% 1: newly computed, 0: not-yet computed.
bUpdate = zeros(1,length(margin_benefit)) ;
for j=1:length(eval_order),

    % the positive pole of this set. 
    cur_pos_pole = EvalList(eval_order(j)) ;

    % (1) incremental inverse matrix
    if lin_mode==0,

        potential = zeros(N,1) ;

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
        potential(loc_uind) = pot ; 
        potential(pos_pole) = 1 ;
        potential(neg_pole) = 0 ;

    % (2) LU decomposition
    elseif lin_mode == 1,
        [potential] = get_potential_lineq_lu([], lap_mat, ...
            [rank_set cur_pos_pole], neg_pole, ilu_setup);
    % (3) Naive method using backslash.
    else
        [potential] = get_potential_lineq([], lap_mat, ...
            [rank_set cur_pos_pole], neg_pole);
    end

    if ~isempty(massSP), 
        potential(ind_ord) = potential(ind_ord) .* massSP ;
    end

    % update marginal benefit.
    margin_benefit(eval_order(j)) = sum(potential) - prev_max_obj ;

    % change the flag to 1. 
    bUpdate(eval_order(j)) = 1 ;

    % check wheter the vertex with the maximum margin is the one
    % that is newly updated. If it is, it's the maximum. 
    [max_v_j max_i_j] = max(margin_benefit) ;
    if bUpdate(max_i_j)==1, 
        break; 
    end
end


%{
% display
display(['Lazy Greedy: Run only ' num2str(j) ' out of ' num2str(length(eval_order)) ...
     '.']) ;
%}

%{
% (DEBUG)
mark1 = '*r' ; mark2='.b' ;
figure ; hold on;
plot(find(bUpdate==0), margin_benefit(bUpdate==0), mark1);
plot(find(bUpdate==1), margin_benefit(bUpdate==1), mark2);
%}

% the vertex that has maximum gain
max_vtx = EvalList(max_i_j) ; 
% the max potential value
max_val = max_v_j + prev_max_obj ;

end
