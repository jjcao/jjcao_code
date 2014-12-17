%
% Compute a matching between two shapes using the ACO algorithm
%
% Input -
%   - Y1, Y2: input contours or general sets of points to be matched.
%   Y1 is a matrix of dimensions <n1 x d> and Y2 is a matrix of
%   dimensions <n2 x d>, where 'n1' and 'n2' are the number of
%   vertices/points in each shape and d = 2 is the dimension of the
%   points
%   - Dist1, Dist2: structures containing the pairwise distance matrices
%   computed for each contour or set of points. Please refer to the help
%   of function distance_matrix() for the description of these
%   structures
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
function [K, best_cost] = aco_matching(Y1, Y2, Dist1, Dist2, S)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global num_iterations;
global num_ants;
global initial_pheromone;
global pheromone_influence;
global pheromone_persistency;
global pheromone_deposit;
global minimum_pheromone;
global descriptor_influence;
global distance_sim_gaussian;
global descriptor_sim_gaussian;

global order_preserving;
global contours;
global open_contour;
global to_one;

global verbose;
global show_progress;
global show_results;

% Check if user called set_global first
if isempty(initial_pheromone)
    error('Variable initial_pheromone is empty. Did you run set_global?');
end

% Normalize the two contours with respect to their enclosed areas
%Y1 = area_normalize(Y1);
%Y2 = area_normalize(Y2);
% Not needed, since we already did it in the function shape_matching()

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
    if verbose
        disp('Contours were switched');
    end
end

% Print out ACO and other parameters
if verbose
    disp(['Contour lengths: ' num2str(n1) ' ' num2str(n2)]);
    disp('ACO matching parameters:');
    disp(['Number of iterations: ' num2str(num_iterations)]);
    disp(['Number of ants: ' num2str(num_ants)]);
    disp(['Initial pheromone: ' num2str(initial_pheromone)]);
    disp(['Pheromone influence: ' num2str(pheromone_influence)]);
    disp(['Pheromone persistency: ' num2str(pheromone_persistency)]);
    disp(['Pheromone deposit: ' num2str(pheromone_deposit)]);
    disp(['Minimum pheromone: ' num2str(minimum_pheromone) ' (' num2str(minimum_pheromone/n1) ')']);
    disp(['Descriptor influence: ' num2str(descriptor_influence)]);
    disp(['Descriptor similarity gaussian: ' num2str(descriptor_sim_gaussian) ' (' num2str(S.max_value*descriptor_sim_gaussian) ')']);
    disp(['Distance similarity gaussian: ' num2str(distance_sim_gaussian)]);
    disp(['Order preserving matching: ' num2str(order_preserving)]);
    disp(['Matching contours: ' num2str(contours)]);
    disp(['Open contours: ' num2str(open_contour)]);
    disp(['To-one matching: ' num2str(to_one)]);
end


%%%% ACO algorithm

% Create bipartite graph representing the matching
% Init graph as a bipartite graph "adjacency matrix"
G = zeros(n1, n2);
% Deposit initial pheromone on all edges
G = G + initial_pheromone;

% Init best matching and best matching cost
best_matching = zeros(n1, 1);
best_cost = inf;
% Init cost history
cost_history = zeros(num_iterations*num_ants, 1);
cost_history_count = 1;

% For all iterations
for i = 1:num_iterations
    % Init matchings and costs for current iteration
    matchings = zeros(n1, num_ants);
    costs = zeros(1, num_ants);
    % For all ants
    for j = 1:num_ants
        % Construct a new matching
        matching = construct_matching(G, S, Dist1, Dist2);
        % Evaluate the new matching according to the descriptors
        cost = evaluate_matching(G, S, Dist1, Dist2, matching);
        % Store matching and cost for current ant
        matchings(:, j) = matching;
        costs(1, j) = cost;
        % Store cost history
        cost_history(cost_history_count) = cost;
        cost_history_count = cost_history_count + 1;
        % Compare against best matching so far
        if cost < best_cost
            best_matching = matching;
            best_cost = cost;
        end
    end

    % Update pheromones
    G = update_pheromones(G, matchings, costs);

    % Plot progress, if requested
    if show_progress
        % Plot pheromone matrix
        W = G/max(max(G));
        figure(1);
        imagesc(W, [0 1]);
        colormap gray;
        %title('Pheromones');
        title(['Iteration ' num2str(i) ' Cost of best matching so far ' num2str(best_cost)])
        drawnow;
    end
end

% Translate best matching to output matching format
K = zeros(n1, 2);
K(:, 1) = 1:n1;
K(:, 2) = best_matching;


%%%% Plot results
if show_results
    % Plot pheromone matrix
    if ~show_progress
        W = G/max(max(G));
        figure
        imagesc(W, [0 1]);
        colormap gray;
        title('Pheromones (lighter values indicate more pheromones)');
    end

    % Plot history count
    figure;
    len = num_iterations*num_ants;
    blocks = num_iterations/10;
    for i = 1:blocks;
        data(i) = mean(cost_history(((i-1)*(len/blocks)+1):i*(len/blocks)));
        glob(i) = min(cost_history(((i-1)*(len/blocks)+1):i*(len/blocks)));
    end
    plot(1:blocks, data, '-b');
    hold on;
    plot(1:blocks, glob, '-r');
    title('Matching cost');
    xlabel('Block of iterations');
    ylabel('Cost');
    legend('Average for block of iterations', 'Best for block of iterations');
end
