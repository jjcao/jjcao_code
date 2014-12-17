% 
% This function normalizes the area enclosed by a closed 2D contour
% to unit ( = 1).
%
% new_contour = area_normalize(old_contour)
%
% ---------------------------------------------------------------
% (C) Richard (Hao) Zhang (2006)
%
function new_c = area_normalize(old_c)

n = length(old_c);

% get the centroid of the vertices of the contour
cog = sum(old_c) / n;

% get total area
total_contour_area = contour_area(old_c);

% normalize by dividing the square root of the area
ratio = sqrt(abs(total_contour_area));

for i=1:n
  u = old_c(i, :) - cog;
  new_c(i, :) = cog + u/ratio - cog;
end
