% function [newfunc, nneigh] = gaussian_smoothing(verts, func, neighDist, sigma2, kdtree)
% smoothing a scalar or vector function func defined on verts, using
% Gaussian with sigma^2 = sigma2, and find neighbor vertices within
% neighDist via kdtree.
%
% verts: n*m matrix, which are n m-dimentional points, such as n*3 in 3d space
% func: n*m1 matrix, which is a function defined on verts, you can set func=verts to smooth the point set
% neighDist: n*1 distance vector, neighDist(i) is used to collect neighbors of vertex i
% sigma2: n*1 vector, sigma2(i) is sigma^2 for the Gaussian filter of vertex i
%
% nneigh: numbers of neighbors selected for each vertex
%
% jjcao @ 2014
%

% %% matlab implementation
% function [newfunc, nneigh] = gaussian_smoothing(verts, func, neighDist, sigma2, kdtree)
% 
% newfunc=func;
% nneigh = ones(size(verts,1), 1);
% for i=1:size(verts,1)
%     [idxs, dists] = kdtree_ball_query( kdtree, verts(i,:), neighDist(i));
%     nneigh(i) = length(idxs);
%     wvec = exp(-dists.^2/(2*sigma2(i))) / sqrt(2*pi*sigma2(i));
%     wmat = repmat(wvec, 1, size(func,2));
%     newfunc(i,:) = sum(func(idxs,:).*wmat,1)/sum(wvec);
% end