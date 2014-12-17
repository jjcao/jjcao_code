%
% This function updates the pheromone matrix according to a set of
% matchings and their costs
% 
% Input -
%   - G: bipartite graph containing pheromone levels: G(i,j) denotes the
%   level of pheromone on the edge between vertex i on the first
%   contour/set and vertex j on the second contour/set
%   - matchings: a set of matchings: matchings(i, j) represents the
%   matching of the second contour/set of vertex i on the first contour.
%   j is the index of the matching
%   - costs: cost of each matching: cost(1, j) represents the cost of
%   matching j
%
% Output -
%   - Gout: the new pheromone matrix
%
function Gout = update_pheromones(G, matchings, costs)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global pheromone_persistency;
global pheromone_deposit;
global minimum_pheromone;

% Get contour size
n1 = size(G, 1);

% Evaporate pheromone
G = G * pheromone_persistency;

% Deposit new pheromone
for ant = 1:size(costs, 2)
    pheromone_delta = pheromone_deposit/costs(1, ant);
    indices = sub2ind(size(G), 1:size(matchings, 1), matchings(:,ant)');
    G(indices) = G(indices) + pheromone_delta;
end

% Control minimum level of pheromone
G(G < (minimum_pheromone/n1)) = (minimum_pheromone/n1);

% Return modified graph
Gout = G;
