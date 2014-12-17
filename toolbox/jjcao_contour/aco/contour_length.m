% 
% This function computes the length of a contour 
%
%     len = contour_length(X)
%         len = returned contour length
%         X = m x 2 matrix representing the (x,y) coordinates of m
%         vertices ordered counterclockwise, which form the contour
%
function len = contour_length(X)
%
% Copyright (c) 2007 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Compute length of X
len = 0;
for i = 1:(size(X)-1)
    len = len + sqrt(sum((X(i,:)-X(i+1,:)).^2));
end
len = len + sqrt(sum((X(end,:)-X(1,:)).^2));
