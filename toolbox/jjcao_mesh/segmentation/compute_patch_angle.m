function patch_normal = compute_patch_angle(fnormal, patch_faces)
% compute patch normals, normalized
%
% (C) Copyright jjcao, 2012

npatch = length(patch_faces);
patch_normal = zeros(npatch,3);
for i=1:npatch
    pfaces = patch_faces{i};
    patch_normal(i,:) = mean(fnormal(pfaces,:));
end
d = sqrt( sum(patch_normal.^2,2) ); 
d(d<eps)=1;
patch_normal = patch_normal ./ repmat( d, 1,3 );