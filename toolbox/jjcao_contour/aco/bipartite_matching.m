%
% Compute the bipartite matching between two shapes according to the
% Hungarian algorithm
%
% Input -
%   - Y1, Y2: input contours or general sets of points to be matched.
%   Y1 is a matrix of dimensions <n1 x d> and Y2 is a matrix of
%   dimensions <n2 x d>, where 'n1' and 'n2' are the number of
%   vertices/points in each shape and d = 2 is the dimension of the
%   points
%   - S: a matrix of dimensions <n1 x n2> representing the similarity
%   matrix for the specified descriptor. S(i, j) is the similarity
%   between the i-th vertex/point of shape 1 and the j-th vertex/point
%   of shape 2
%
% Output -
%   - K: a matrix of dimensions <n1 x 2> which contains the computed
%   matching: vertex K(i, 1) on shape 1 is matched to vertex K(i, 2) on
%   shape 2
%   - best_cost: cost of computed matching
%
function [K, best_cost] = bipartite_matching(Y1, Y2, S)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Call bipartite matching code
[C, best_cost] = hungarian(S.value);

% Translate matching to the common format
K = zeros(size(Y1, 1), 2);
K(:, 1) = 1:size(Y1, 1);
K(:, 2) = C;
