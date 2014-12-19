%
% This function computes the range of valid assignments for a vertex,
% according to order preservation.
%
% Input -
%   - vertex: index of vertex on first contour that should be matched to
%   vertices on the second contour
%   - matching: matching that has been computed so far: vertices that
%   are still unassigned are matched to 0
%   - n2: size of second contour
%
% Output -
%   - inrange: array of valid ranges: inrange(i) is 1 if vertex i of the
%   second contour is in the valid range of assignments for vertex
%
function inrange = valid_range(vertex, matching, n2)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get first contour size
n1 = size(matching, 1);

% Init valid range
inrange = zeros(n2, 1);

% Find left vertex
left_vertex = find(matching(1:(vertex-1)), 1, 'last');
% If not found, loop around
if isempty(left_vertex)
    left_vertex = find(matching((vertex+1):n1), 1, 'last');
    left_vertex = vertex + left_vertex;
end

% Find right vertex
right_vertex = find(matching((vertex+1):n1), 1, 'first');
% If not found, loop around
if isempty(right_vertex)
    right_vertex = find(matching(1:(vertex-1)), 1, 'first');
else
    % Get actual position on the vector
    right_vertex = vertex + right_vertex;
end

% Find corresponding vertices
left_vertex = matching(left_vertex);
right_vertex = matching(right_vertex);

% Create vector with valid range information
if left_vertex == right_vertex
    inrange(left_vertex) = 1;
elseif left_vertex < right_vertex
    % A standard inverval
    inrange((left_vertex+1):right_vertex) = 1;
else
    % A broken interval (contains the origin of the contour)
    inrange((left_vertex+1):n2) = 1;
    inrange(1:right_vertex) = 1;
end
