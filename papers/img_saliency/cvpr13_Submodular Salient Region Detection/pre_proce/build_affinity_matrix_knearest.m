function [adj_mat, sigma] = build_affinity_matrix_knearest(sp_fea, K)
%

sp_num = size(sp_fea,1);

tree = kdtree_build(sp_fea);
edges = zeros(sp_num * (K-1), 2);
sigma = zeros(sp_num,1);
for jj = 1:sp_num
%     jj
    en = jj * (K-1);
    st = en - K + 2;
    edges(st:en,1) = jj;
    [tmp, dist] = kdtree_k_nearest_neighbors(tree,sp_fea(jj,:),K);
    edges(st:en,2) = tmp(2:end);
    sigma(jj) = dist(end);
end
kdtree_delete(tree);

adj_mat = zeros(sp_num);
ind = sub2ind(size(adj_mat), edges(:,1), edges(:,2));
adj_mat(ind) = 1;
%     adj_mat = adj_mat + adj_mat'; adj_mat(adj_mat>0) = 1;

