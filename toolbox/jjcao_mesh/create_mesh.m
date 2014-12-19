function [vertex, faces] = create_mesh(reso)
%
% We create a simple mesh with fine scale details.
% vertex: 3*n, faces: 3*m
%
%   Copyright (c) 200? Gabriel Peyre

%% We generate point on a square.
if nargin < 1
    p = 150;
else
    p = reso;
end

[Y,X] = meshgrid(linspace(-1,1,p),linspace(-1,1,p));
vertex0 = [X(:)'; Y(:)'; zeros(1,p^2)];
n = size(vertex0,2);

%% We generate a triangulation of a square.
I = reshape(1:p^2,p,p);
a = I(1:p-1,1:p-1); b = I(2:p,1:p-1); c = I(1:p-1,2:p);
d = I(2:p,1:p-1); e = I(2:p,2:p); f = I(1:p-1,2:p);
faces = cat(1, [a(:) b(:) c(:)], [d(:) e(:) f(:)])';

%% Width and height of the bumps.
sigma = .03;
h = .35;
q = 8;

%% Elevate the surface using bumps.
t = linspace(-1,1,q+2); t([1 length(t)]) = [];
vertex = vertex0;
for i=1:q
    for j=1:q
        d = (X(:)'-t(i)).^2 + (Y(:)'-t(j)).^2;
        vertex(3,:) = vertex(3,:) + h * exp( -d/(2*sigma^2)  );
    end
end
