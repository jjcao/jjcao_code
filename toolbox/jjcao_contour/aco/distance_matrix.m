%
% Compute the pairwise distances between vertices/points of a 2D shape
%
% D = distance_matrix(Y)
%
% Input -
%   - Y: input shape (contour or general point set). Y is a matrix of
%   dimensions <n x d>, where 'n' is the number of vertices/points in
%   the shape and d = 2 is the dimension of the points
%
% Output -
%   - D: structure containing the pairwise distances for the shape. The
%   structure has two fields: 'value' and 'max_value'.  'value' is a
%   matrix of dimensions <n x n>, where value(i, j) is the distance
%   between vertices/points 'i' and 'j' of the shape.  'max_value'
%   stores the maximum distance. If the global option 'contours' is set
%   to 1, pairwise geodesic distances (distances along the contour) are
%   computed. Otherwise, pairwise Euclidean distances are computed,
%   assuming that we have an arbitrary set of 2D points without
%   connectivity. The global variable 'open_contour' specifies whether
%   the contour is open or closed (connected from vertex 1 to vertex
%   'n' or not).
%
function D = distance_matrix(Y)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global contours;
global open_contour;

% Compute geodesic distance along the shape for contours, or Euclidean
% distance for general point sets
if contours
    % Compute geodesic distance
    M = pairwise_geodesic_dist(Y, open_contour);

    % Create distance structure
    D = struct;
    D.max_value = max(max(M));
    D.value = M;
else
    % Compute pairwise Euclidean distances
    M = zeros(size(Y, 1), size(Y, 1));
    for i = 1:size(Y, 1)
        for j = 1:size(Y, 2)
            M(i, j) = sqrt(sum((Y(i,:) - Y(j,:)).^2));
        end
    end

    % Create distance structure
    D = struct;
    D.max_value = max(max(M));
    D.value = M;
end
