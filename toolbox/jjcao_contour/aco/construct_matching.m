%
% This function constructs a matching for the ACO algorithm
%
% Input -
%   - G: bipartite graph containing pheromone levels: G(i,j) denotes the
%   level of pheromone on the edge between vertex i on the first
%   contour/set and vertex j on the second contour/set
%   - S: descriptor dissimilarity matrix: S(i, j) denotes the distance
%   between the descriptors of vertex i on contour/set 1 and vertex j on
%   contour/set 2
%   - Dist1, Dist2: square matrices of geodesic or pairwise distances:
%   Dist1.value(i, j) denotes the distance between the vertices i and j
%   of contour/set 1
%
% Output -
%   - matching: the constructed matching: vertex i on the first contour
%   matches to vertex matching(i) on the second contour
%
function matching = construct_matching(G, S, Dist1, Dist2)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global pheromone_influence;
global descriptor_sim_gaussian;
global distance_sim_gaussian;

global order_preserving;
global to_one;

% Get similarity gaussians
dist_sigma = distance_sim_gaussian;
%dist_sigma = ((Dist1.max_value+Dist2.max_value)/2) * distance_sim_gaussian;
%((Dist1.max_value+Dist2.max_value)/2) is 1, since the distance are
%normalized
desc_sigma = S.max_value * descriptor_sim_gaussian;

% Compute the gaussian exponent denominators
dist_den = 2*(dist_sigma^2);
desc_den = 2*(desc_sigma^2);

% Get contour sizes
n1 = size(G, 1);
n2 = size(G, 2);

% Init new matching
matching = zeros(n1, 1);

% Reset number of matched vertices in first contour
count = 0;

% Init probabilities
p = zeros(n2, 1);

% Init matching history
i_2 = 0;
j_2 = 0;
i_3 = 0;
j_3 = 0;

% Generate a permutation of the vertices in the first contour
options = randperm(n1);

% Select first vertex
vertex1 = options(1);

% Init visited information
visited = zeros(n2, 1);

% By default, all vertices are in the valid range
inrange = ones(n2, 1);

% Visit all vertices in the first contour
while count < n1
    % Compute valid vertex range for second contour
    if order_preserving
        if count < 2
            % By default, all vertices are in the valid range
            inrange = ones(n2, 1);
        else
            inrange = valid_range(vertex1, matching, n2);
        end
    end

    % Compute probability of visiting each vertex in the second contour
    if to_one
        condition = (inrange == 1) & (visited == 0);
        if (sum(condition) == 0)
            condition = (inrange == 1);
        end
        len = sum(condition);
    else
        condition = (inrange == 1);
        len = sum(condition);
    end

    % Compute distance similarity
    if count > 0
        dist_sim1 = abs(Dist1.value(vertex1, i_2) - Dist2.value(condition, j_2)) * exp(-(Dist1.value(vertex1, i_2)^2)/dist_den);
    else
        dist_sim1 = zeros(len, 1);
    end
    dist_sim1 = (1.0 - dist_sim1);

    if count > 1
        dist_sim2 = abs(Dist1.value(vertex1, i_3) - Dist2.value(condition, j_3)) * exp(-(Dist1.value(vertex1, i_3)^2)/dist_den);
    else
        dist_sim2 = zeros(len, 1);
    end
    dist_sim2 = (1.0 - dist_sim2);
    
    % Compute descriptor similarity
    desc_sim = exp(-((S.value(vertex1, condition)).^2)./desc_den)';
    
    % Compute final heuristic value
    heuristic = dist_sim1.*dist_sim2.*desc_sim;

    %%%% Compute traversal probability
    p(~condition) = 0;
    p(condition) = (1.0 - pheromone_influence)*heuristic + ...
                   pheromone_influence*(G(vertex1, condition)');
    
    % Probability should never be zero
    %p(p == 0.0) = 0.000001;
    % This does not make sense anymore when order preservation is used
    
    % Normalize probabilities between [0..1]
    total = sum(p);
    p = p / total;
    
    % Select matched vertex according to probabilities
    s = cumsum(p);
    % See where the random number falls in the probability intervals
    r = rand();
    vertex2 = find(s > r, 1);
    % Check whether the selected vertex is in the valid range
    while ~inrange(vertex2)
        r = rand();
        vertex2 = find(s > r, 1);
    end
    
    % Assign new matching
    matching(vertex1) = vertex2;
    visited(vertex2) = 1;

    % Increment matching size
    count = count + 1;
    
    % Update matching history
    i_3 = i_2;
    j_3 = j_2;
    i_2 = vertex1;
    j_2 = vertex2;

    % Select next vertex in first contour
    if count < n1
        vertex1 = options(count+1);
    end
end
