%
% This function plots a contour specified by a list of 2D points
% given in the n by 2 array c.
% 
%    show_contour(c, varargin)
% 
% the clr (contour color) argument, e.g., clr = 'blue', is
% optional. The default is black.
%
% ---------------------------------------------------------------
% (C) Richard (Hao) Zhang (2006)
%
function show_contour(c, varargin)

n = length(c);

% get color
if nargin > 1
    clr = varargin{1};
else
    clr = 'black';
end

% plot contour
plot([c(n,1) c(1,1)], [c(n,2) c(1,2)], 'o-');
hold;

for i=1:n-1
   plot([c(i,1) c(i+1,1)], [c(i,2) c(i+1,2)], 'o-');
end

% set color
set(findobj('Type','line'), 'Color', clr)

% playing around with the axis to give better display
x_min = min(c(:,1));
x_max = max(c(:,1));
y_min = min(c(:,2));
y_max = max(c(:,2));
x_step = (x_max - x_min)/10;
y_step = (y_max - y_min)/10;
axis([x_min-x_step x_max+x_step y_min-y_step y_max+y_step]);
axis equal;

hold off;
