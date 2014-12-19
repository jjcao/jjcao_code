% 
% This function computes the area enclosed by a closed 2D contour.
%
% total_contour_area = contour_area(contour)
%
% ---------------------------------------------------------------
% (C) Richard (Hao) Zhang (2005)
%
function total_contour_area = contour_area(c)

n = length(c);

% get centroid of points
cog = sum(c) / n;

% get total area
total_contour_area = signed_triangle_area(cog, c(n,:), c(1,:));
for i=1:n-1
  total_contour_area = total_contour_area + signed_triangle_area(cog, c(i, :), c(i+1, :));
end
