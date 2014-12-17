function spAffinityMat = gene_weight( spAdjacentMat, sp_fea, sigma)
% 
% for cvpr13_Submodular Saliency
%
%

ind = find(spAdjacentMat>0);
[row,col] = ind2sub(size(spAdjacentMat), ind); 
valDistances = sum( (sp_fea(row,:) - sp_fea(col,:)).^2, 2)+eps;
beta = sigma(row).*sigma(col);
beta = 1./beta;
weights=exp(-beta.*valDistances);
spAffinityMat = zeros(size(spAdjacentMat));
spAffinityMat(ind) = weights;

