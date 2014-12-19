function Dc_mat = compute_Dc(sregion,sp_fea)
num_region = length(sregion);
Dc_mat = zeros(num_region);
num_sp = size(sp_fea,1);

fea_dist_mat = zeros(num_sp);
ind = 1:num_sp^2;
[row, col] = ind2sub(size(fea_dist_mat), ind);
fea_dist_mat(ind) = sqrt( sum( (sp_fea(row,:) - sp_fea(col,:)).^2, 2));

for i = 1:num_region
    sri = sregion{i};
    for j = 1:num_region
        srj = sregion{j};
        tmp = fea_dist_mat(sri, srj);        
        Dc_mat(i,j) = mean(tmp(:));
    end
end