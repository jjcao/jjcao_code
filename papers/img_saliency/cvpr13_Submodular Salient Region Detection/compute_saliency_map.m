function [f, sregion] = compute_saliency_map(fa_location, fa_assig, sp_fea, sp_center)

num_region = length(fa_location);
sregion = cell(num_region,1);
sp_nums = zeros(num_region,1);    
for i = 1:num_region
    sregion{i} = find( fa_assig == fa_location(i));   
    sp_nums(i) = length(sregion{i});
end

%
Dc_mat = compute_Dc(sregion,sp_fea);
V = compute_V(sregion, sp_center);
maxV = max(V);

%
fc = ones(num_region,1);
fs = fc;
for i = 1:num_region        
    fc(i) = sum(Dc_mat(i,:).*sp_nums');
    fs(i) = 1 - V(i)/maxV;
end
fc=(fc-min(fc))/(max(fc)-min(fc));
fs=(fs-min(fs))/(max(fs)-min(fs));
f = fc.*fs; % f is saliency for each region