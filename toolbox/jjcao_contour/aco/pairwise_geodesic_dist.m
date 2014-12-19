%
% Compute pairwise geodesic distances for an open or closed 2D contour
%
% G = pairwise_geodesic_dist(Y, open_contour)
%
% Input -
%   - Y: input contour. Y is a matrix of dimensions <n x d>, where 'n'
%   is the number of vertices on the contour and d = 2 is the dimension
%   of the points
%   - open_contour: if open_contour = 1, we assume that we have an open
%   contour (vertices 1 and 'n' are disconnected). If open_contour = 0,
%   we assume that we have a closed contour (vertex 1 is connected to
%   vertex 'n')
%
% Output -
%   - G: matrix of dimensions <n x n> containing the pairwise geodesic
%   distances, where G(i, j) is the distance between vertices 'i' and
%   'j' on the contour
%
function G = pairwise_geodesic_dist(Y, open_contour)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get contour length
n = size(Y, 1);

% Init matrix 'G' and temporary distance vector 'D'
G = zeros(n, n);
D = zeros(n, 1);

% Compute geodesic distance from first vertex to all others
s = 0.0;
D(1) = 0.0;
for j = 1:(n-1)
    s = s + norm(Y(j,:) - Y(j+1,:), 2);
    D(j+1) = s;
end
s = s + norm(Y(n,:) - Y(1,:), 2);

% Compute pairwise distances
for i = 1:n
    for j = i:n
        if i ~= j
            dist1 = D(j) - D(i);
            if open_contour
                G(i, j) = dist1;
            else
                dist2 = s - dist1;
                G(i, j) = min(dist1, dist2);
            end
            G(j, i) = G(i, j);
        end
    end
end

% Normalize distances between 0 and 1
G = G/max(max(G));
