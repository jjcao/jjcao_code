function h = get_cij_prior(aff_mat,prior,label_inds,sp_num,C_param)
C_param = 0.05;
 all_node = 1:sp_num;
 unlabel_inds = setdiff(all_node, label_inds);
 label = label_inds;
 label_inds = [label_inds sp_num+1];
 D = diag(1.0./sum(aff_mat,2));
 P = D*aff_mat;
 P_UU = P(unlabel_inds,unlabel_inds);
 P_UL = P(unlabel_inds,label_inds);
 h_L = [1,0];
 h_Q = 1-prior(unlabel_inds);
 I_UU = eye(numel(unlabel_inds));
 h_U = inv(I_UU-(1-C_param)*P_UU)*((1-C_param)*P_UL*h_L'+C_param*h_Q);

 h = zeros(sp_num,1);
 h(unlabel_inds) = h_U;
 h(label) = 1;