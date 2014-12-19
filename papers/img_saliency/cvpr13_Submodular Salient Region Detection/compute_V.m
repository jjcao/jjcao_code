function V = compute_V(sregion, sp_center)

num_region = length(sregion);
region_center = zeros(num_region,2);
for i = 1:num_region    
    ri = sregion{i};
    spc_region = sp_center(ri,:);
    region_center(i,:) = mean(spc_region);
end

V = zeros(num_region,1);
for i = 1:num_region    
    ri = sregion{i};    
    spc_region = sp_center(ri,:);    
    for k =1:num_region
        rc = region_center(k,:);        
        rc = repmat(rc, size(spc_region,1), 1);
        V(i) = V(i) + mean( sqrt(sum((rc - spc_region).^2, 2)) );        
    end
end