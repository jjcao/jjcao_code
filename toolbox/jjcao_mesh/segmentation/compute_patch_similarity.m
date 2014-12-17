function A = compute_patch_similarity(patch_adjancy, A_iso, A_angle, alfa)
%
% (C) Copyright jjcao, 2012
ind = find(patch_adjancy);
avg_iso = mean(A_iso(ind));
avg_angle = mean(A_angle(ind));

[row,col] = ind2sub(size(A_iso),ind);
w1 = exp(-A_iso(ind)./avg_iso*4);
w2 = exp(-A_angle(ind)./avg_angle*4);
% A = sparse(row,col,alfa*w1+(1-alfa)*w2);
A = sparse(row,col, min(w1,w2));

% A = exp(-(alfa*A_iso./avg_iso + (1-alfa)*A_angle./avg_angle*4)); 
