%
% Compute an order preserving matching between two shapes according to
% the COPAP algorithm
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
function [K, best_cost] = order_preserving_matching(Y1, Y2, S)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Normalize the two contours with respect to their enclosed areas
Y1 = area_normalize(Y1);
Y2 = area_normalize(Y2);

% Get contour sizes
n1 = length(Y1(:,1));
n2 = length(Y2(:,1));

% Switch contours if n1 > n2
if n1 > n2
    temp_n = n1;
    n1 = n2;
    n2 = temp_n;
    temp_Y = Y1;
    Y1 = Y2;
    Y2 = temp_Y;
end

% Compute cost of dummy matching
dummy_cost = S.max_value + 1;

% Create new dissimilarity matrix, if lengths are different
if n1 ~= n2
    [m, n] = size(S.value);
    newS = zeros(max(m,n), max(m,n)) + dummy_cost;
    newS(1:size(S.value,1), 1:size(S.value,2)) = S.value;
    S.value = newS;
end

% Perform order preserving assignment
% Syntax for old code
%[P, best_cost, iter] = copap(S.value, n1, n2, dummy_cost);
% Syntax for new code
[P, best_cost] = copap(S.value, n1, dummy_cost, 1);

% Create new matching matrix
K = zeros(size(P,2), 2);
% P is switched since rows in S correspond to elements in Y1 and columns
% to elements in Y2
% P(i) is the row assigned to the i-th column
K(:,2) = 1:size(P,2);
K(:,1) = P';
