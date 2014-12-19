function [fa_location, fa_assig,K] = submodular_salient_region_detection(aff_mat,prior ,K, lambda, pos_label_inds,sp_num,C_param)

%GuangyuZhong 05/10/2013
% 
ground_cond = 0.1.*ones(sp_num, 1); % tune this parameters
aff_mat = [[aff_mat, ground_cond]; [ground_cond', 0]] ;
% aff_mat = [[aff_mat, ground_cond]; [ones(1,sp_num), 0]] ;

% C = get_cij(aff_mat, pos_label_inds, ones(numel(pos_label_inds), 1));

C = get_multi_cij(aff_mat,prior,pos_label_inds,sp_num,C_param);
set_i = 1:sp_num;
num_pos = numel(pos_label_inds);
% [C] = get_cij(aff_mat, pos_label_inds, ones(num_pos,1));

fa_location = [];
fa_assig = zeros(sp_num,1);
pi_cur = zeros(sp_num,1);
obj_val = [0];
if K > num_pos;
    K = num_pos;
    disp(['new K  number: ' num2str(K)]);
end


for  i = 1:K
    if i==1
        for m = 1:num_pos
            H(m) = sum(C(:,m)) - lambda*i -obj_val(i);
        end  
        ulabel = pos_label_inds;
    else
        for u = 1:numel(ulabel)
            index = find(pos_label_inds==ulabel(u));
            H(u) = sum(max(pi_cur,C(:,index)) )- lambda*i - obj_val(i);
        end
    end
    
    [max_val, max_ind] = max(H);
    max_val = max_val+obj_val(i);
    if max_val<=obj_val(i)
        break;
    end
    clear H;
    obj_val = [obj_val max_val];
    fa_location = [fa_location ulabel(max_ind)];
%     ulabel = setdiff(pos_label_inds, fa_location);
ulabel(max_ind) = [];
    pi_cur(fa_location(i)) = 1;
    set_i_A = setdiff(set_i, fa_location);
    index = find(pos_label_inds==fa_location(i));
    
    for n = 1:sp_num-i        
        if pi_cur(set_i_A(n))<C(set_i_A(n),index)
            pi_cur(set_i_A(n)) = C(set_i_A(n),index);
            fa_assig(set_i_A(n)) = fa_location(i);
        end
    end
    fa_assig(fa_location) = fa_location;
end

figure;plot(obj_val);



function [C] = get_cij(aff_mat, label_inds, label_ranks)
if size(label_inds, 1) == 1
    label_inds = label_inds';
end
if size(label_ranks, 1) == 1
    label_inds = label_ranks';
end
num_label = numel(label_inds);

for j = 1:num_label
    ind = label_inds(j);
    rank = label_ranks(j);
    C(:,j) = solve_les_with_dirichlet_no_prior(aff_mat, rank, ind);

end

function [C] = get_multi_cij(aff_mat,prior,label_inds,sp_num,C_param)
for j = 1:numel(label_inds)
    ind = label_inds(j);
    
 C(:,j) = get_cij_prior(aff_mat,prior,ind,sp_num,C_param);
end


