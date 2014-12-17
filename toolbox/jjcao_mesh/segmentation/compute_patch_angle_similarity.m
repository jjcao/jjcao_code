function A =compute_patch_angle_similarity(patch_adjancy,fnormal,patch_faces)
% compute patch similarity by patch's normal
%
% (C) Copyright jjcao, 2012

patch_normal = compute_patch_angle(fnormal, patch_faces);

ind = find(patch_adjancy);
[row,col] = ind2sub(size(patch_adjancy),ind);
d1 = sum( (patch_normal(row,:) - patch_normal(col,:)).^2, 2);
A = sparse(row,col,d1);