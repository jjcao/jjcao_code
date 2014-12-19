function S = compute_S(fa_location, f, sp_fea, sigma)

num_region = length(f);
num_sp = size(sp_fea,1);
W = zeros(num_sp, num_region);
ind = 1:num_sp*num_region;
[row, col] = ind2sub(size(W), ind);

valDistances = sum( (sp_fea(row,:) - sp_fea(fa_location(col),:)).^2, 2)+eps;

% if length(sigma)>1
% if isempty(sigma)    
%     beta = sigma(row).*sigma(col);
%     beta = 1./beta;
    beta = 1.0/max(min(valDistances), sigma*mean(valDistances));
    W(ind)=exp(-beta.*valDistances);
% else    
%     minVal = min(valDistances);
%     valDistances=(valDistances-minVal)/(max(valDistances)-minVal);%Normalize to [0,1]
%     W(ind)=exp(-sigma*valDistances);        
% end

S = W * f;

%% % old
% function S = compute_S(f, sp_fea, sregion, theta)
% 
% num_region = length(f);
% num_sp = size(sp_fea,1);
% W = zeros(num_sp);
% ind = 1:num_sp^2;
% % W = zeros(num_sp, num_region);
% % ind = 1:num_sp*num_region;
% [row, col] = ind2sub(size(W), ind);
% 
% tmp = sum( (sp_fea(row,:) - sp_fea(col,:)).^2, 2);
% valDistances = sqrt(tmp)+eps;
% valDistances=normalize(valDistances); %Normalize to [0,1]
% W(ind)=exp(-theta*valDistances);
% 
% F = zeros(num_sp,1);
% for i=1:num_region
%     F(sregion{i}) = f(i);
% end
% 
% S = W * F;