% 
% This function computes the cost for a given matching
%
% Input -
%   - G: bipartite graph containing pheromone levels: G(i,j) denotes the
%   level of pheromone on the edge between vertex i on the first contour
%   and vertex j on the second contour
%   - S: descriptor dissimilarity matrix: S(i, j) denotes the distance
%   between the descriptors of vertex i on contour 1 and vertex j on
%   contour 2
%   - Dist1, Dist2: square matrices of geodesic or pairwise distances:
%   Dist1.value(i, j) denotes the distance between the vertices i and j
%   of contour/set 1
%   - matching: the constructed matching: vertex i on the first
%   contour/set matches to vertex matching(i) on the second contour/set
%
% Output -
%   - cost: value representing the cost of the matching
%
function cost = evaluate_matching(G, S, Dist1, Dist2, matching);
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global descriptor_influence;
global descriptor_sim_gaussian;
global distance_sim_gaussian;

% Get contour sizes
n1 = size(G, 1);
n2 = size(G, 2);

% Get similarity gaussians
dist_sigma = distance_sim_gaussian;
%dist_sigma = ((Dist1.max_value+Dist2.max_value)/2) * distance_sim_gaussian;
%((Dist1.max_value+Dist2.max_value)/2) is 1, since the distance are
%normalized
desc_sigma = S.max_value * descriptor_sim_gaussian;

% Compute the gaussian exponent denominators
dist_den = 2*(dist_sigma^2);
desc_den = 2*(desc_sigma^2);

% Compute difference of descriptors (descriptor cost)
indices = sub2ind(size(S.value), 1:size(matching, 1), matching');
desc_cost = sum((exp(-((S.value(indices)).^2./desc_den))))/n1;
desc_cost = desc_cost;

% Compute distance of pairs (distance cost)
dist_cost = 0;
for i = 1:size(matching, 1)
    for j = 1:size(matching, 1)
      % Distances are already normalized between [0..1]
      dist = abs(Dist1.value(i, j) - Dist2.value(matching(i), matching(j)));
      dist_cost = dist_cost + dist*(exp(-((Dist1.value(i, j)^2)/dist_den)));
    end
end
dist_cost = dist_cost / ((n1*(n1-1))/2);
dist_cost = dist_cost;

% Compute total cost
cost = descriptor_influence*(1.0 - desc_cost) + ...
       (1.0 - descriptor_influence)*(dist_cost);
