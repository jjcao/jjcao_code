%
% Compute a shape descriptor for a given shape
%
% D = extract_descriptor(Y, descriptor, feature);
%
% Input -
%   - Y: input shape (contour or general point set). Y is a matrix of
%   dimensions <n x d>, where 'n' is the number of vertices/points in
%   the shape and d = 2 is the dimension of the points
%   - descriptor: name of shape descriptor to be computed. In this
%   release, only 'shape_context' is supported
%   - feature: name of the shape descriptor feature. In this release, no
%   feature is used (so, the value '' should be provided).
%
% Output -
%   - D: matrix of dimensions <n x l> representing the shape descriptor,
%   where 'l' is the dimensionality of the descriptor, derived from its
%   parameters
%
function D = extract_descriptor(Y, descriptor, feature);
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Get global variable
global show_results;

% Extract descriptor
if strcmp(descriptor, 'shape_context')
    D = extract_shape_context(Y);
else
    disp(['Descriptor "' descriptor '" is not defined!']);
    return;
end
