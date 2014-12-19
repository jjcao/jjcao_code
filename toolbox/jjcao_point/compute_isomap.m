function [vec val] = compute_isomap(D, ndim)
%
%   Copyright (c) 2012 Junjie Cao

N = size(D,1);
opt.disp = 0; opt.isreal = 1; opt.issym = 1; 
ND = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));

[vec, val] = eigs(ND, ndim, 'LR', opt);
val = diag(val);
for j=1:ndim
    vec(:,j) = vec(:,j)*sqrt(val(j));
end  