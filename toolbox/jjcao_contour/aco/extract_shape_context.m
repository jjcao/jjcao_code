%
% Compute shape context descriptor
%
% D = extract_shape_context(Y)
%
% Input:
%   - Y: input shape (contour or general point set). Y is a matrix of
%   dimensions <n x d>, where 'n' is the number of vertices/points in
%   the shape and d = 2 is the dimension of the points
%
% Output -
%   - D: matrix of dimensions <n x l> representing the shape descriptor,
%   where 'l' is the dimensionality of the descriptor, derived from its
%   parameters
%
function D = extract_shape_context(Y)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%
    
% Normalize the contour with respect to its enclosed area
%Y = area_normalize(Y);
% No need to normalize because we already do it in the shape_matching()
% function

% Get contour size
n = length(Y(:,1));

% Set shape context parameters
nbins_theta=12;
nbins_r=5;
r_inner=1/8;
r_outer=2;
out_vec = zeros(1, n);

% Extract shape context
[D, mean_dist] = sc_compute(Y', zeros(1, n), [], ...
                        nbins_theta, nbins_r, r_inner, r_outer, out_vec);

% Show shape contexts (change 0 to 1, for debug)
if 0
    figure
    temp = Y;
    temp(end+1,:) = Y(1,:);
    plot(temp(:,1), temp(:,2), 'k+-');
    axis equal;
    hold on;
    for i = 1:size(Y, 1)
        N=256; t=(0:N)*2*pi/N;
        plot(cos(t)*r_outer + Y(i, 1), sin(t)*r_outer + Y(i, 2), 'b');
        set(text(Y(i, 1), Y(i, 2), num2str(sum(D(i,:)))), 'color', 'r');
    end
end
