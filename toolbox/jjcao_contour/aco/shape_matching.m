%
% Compute the matching between two 2D shapes (contours or general sets
% of points)
%
% [K, S, best_cost, D1, D2, Dist1, Dist2] = shape_matching(Y1, Y2, [algorithm [, descriptor [, feature [, distance]]]])
%
% Input:
%   - Y1, Y2: input contours or general sets of points to be matched.
%   Y1 is a matrix of dimensions <n1 x d> and Y2 is a matrix of
%   dimensions <n2 x d>, where 'n1' and 'n2' are the number of
%   vertices/points in each shape and d = 2 is the dimension of the
%   points
%   - algorithm (optional): name of the shape matching method. In this
%   implementation, three algorithms are supported: 'aco', 'bipartite',
%   and 'order_preserving'
%   - descriptor (optional): name of the shape descriptor. In this
%   implementation, only 'shape_context' is supported
%   - feature (optional): name of the shape descriptor feature. In this
%   implementation, no feature is used (so, the value '' should be
%   provided).
%   - distance (optional): metric used to compute the distance between
%   two shape descriptors. Two metrics are available: 'euclidean' or
%   'chisquare'. The default value is 'chisquare'.
%
% Output:
%   - K: a matrix of dimensions <n1 x 2> which contains the computed
%   matching: vertex K(i, 1) on shape 1 is matched to vertex K(i, 2) on
%   shape 2
%   - S: a matrix of dimensions <n1 x n2> representing the similarity
%   matrix for the specified descriptor. S(i, j) is the similarity
%   between the i-th vertex/point of shape 1 and the j-th vertex/point
%   of shape 2
%   - best_cost: cost of computed matching (cost of K according to
%   the specified algorithm)
%   - D1, D2: descriptors extracted for each shape. D1 and D2 are
%   matrices of dimensions <n1 x l> and <n2 x l>, respectively, where
%   'l' is the dimensionality of the descriptor, derived from its
%   parameters
%   - Dist1, Dist2 (optional): structures containing the pairwise
%   distance matrices computed for each contour or set of points. Please
%   refer to the help of function distance_matrix() for the description
%   of these structures
%
% This function computes the matching between two 2D contours or two
% sets of 2D points by using the specified algorithm and descriptor. The
% parameters for this function are set in the script 'set_global.m'
%
function [K, S, best_cost, D1, D2, Dist1, Dist2] = shape_matching(Y1, Y2, varargin)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variables
global contours;
global verbose;
global show_progress;
global show_results;

% Check for additional parameters
algorithm = 'aco';
if nargin > 2
    algorithm = varargin{1};
end

descriptor = 'shape_context';
if nargin > 3
    descriptor = varargin{2};
end

feature = 'curvature';
if nargin > 4
    feature = varargin{3};
end

comparison = 'chisquare';
if nargin > 5
    comparison = varargin{4};
end

% Normalize the two contours with respect to their enclosed areas
Y1 = area_normalize(Y1);
Y2 = area_normalize(Y2);

% Get contour sizes
n1 = length(Y1(:, 1));
n2 = length(Y2(:, 1));

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

% Compute shape descriptors for contours
if verbose
    if strcmp(descriptor, 'shape_context')
        disp(['Extracting descriptor: ' descriptor]);
    else
        disp(['Extracting descriptor: ' descriptor ' (' feature ')']);
    end
end
D1 = extract_descriptor(Y1, descriptor, feature);
D2 = extract_descriptor(Y2, descriptor, feature);

% Compute similarity matrix
if verbose
    disp(['Similarity matrix constructed with distance: ' comparison]);
end
if strcmp(comparison, 'chisquare')
    S = simmat_chisquare(D1, D2);
elseif strcmp(comparison, 'euclidean')
    S = simmat_euclidean(D1, D2);
else
    disp(['Distance measure "' comparison '" is not defined!']);
    return
end

% Compute geodesic distance matrices
if strcmp(algorithm, 'aco')
    if verbose
        if contours
            disp('Computing geodesic distances');
        else
            disp('Computing pairwise distances');
        end
    end
    Dist1 = distance_matrix(Y1);
    Dist2 = distance_matrix(Y2);
end


%%%% Run matching algorithm
if strcmp(algorithm, 'aco')
    if verbose
        disp('ACO matching');
    end
    [K, best_cost] = aco_matching(Y1, Y2, Dist1, Dist2, S);
elseif strcmp(algorithm, 'order_preserving')
     if verbose
        disp('Order preserving matching');
    end
    [K, best_cost] = order_preserving_matching(Y1, Y2, S);
elseif strcmp(algorithm, 'bipartite')
    if verbose
        disp('Bipartite matching');
    end
    [K, best_cost] = bipartite_matching(Y1, Y2, S);
else
    disp(['Algorithm "' algorithm '" is not defined!']);
    return
end


%%%% Plot results
if show_results
    % Plot similarity matrix
    W = S.value/S.max_value;
    figure
    imagesc(W, [0 1]);
    colormap gray;
    title('Similarity matrix for descriptor');

    % Visualize matching in four different ways
    figure;
    viz_matching(Y1, Y2, K, 'lines', 'centered');
    title('Best matching');

    figure;
    viz_matching(Y1, Y2, K, 'lines', 'xdisplaced');
    title('Best matching');
    
    figure;
    viz_matching(Y1, Y2, K, 'lines', 'displaced');
    title('Best matching');

    figure;
    viz_matching(Y1, Y2, K, 'labels', 'xdisplaced');
    title('Best matching');
end
