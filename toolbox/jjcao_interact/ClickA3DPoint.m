% ClickA3DPoint.m
% by Babak Taati
% http://qlink.queensu.ca/~3bt1/
% Queen's University
% May 4, 2005 (revised Oct 30, 2007)
%
% function    ClickA3DPoint( PointCloud );
%
%   This function shows a 3D point cloud and lets the user click select one 
%   of the points by clicking on it. The selected point will be highlighted 
%   and its index in the point cloud will be printed. 
%
% input:
%   PointCloud:  should be a 3*N matrix representing N 3D points
%                (if it's M*N, M > 3, the higher dimensions are ignored)
%
% output:
%   none
%
% other functions required:
%   callbackClickA3DPoint: the mouse click callback function
%   RowNorm: returns norms of each row of a matrix
%
%       
%   To test this function, run these lines:
%
%                 PointCloud = rand(3,100)*100;
%                 ClickA3DPoint( PointCloud );
%
%   now rotate or move the point cloud and try it again.
%   (on the figure View menu, turn the Camera Toolbar on, ...)
%

function    ClickA3DPoint( PointCloud );

plot3(PointCloud(1,:), PointCloud(2,:), PointCloud(3,:), 'c.'); % visualize the point cloud
hold on; % so we can highlight the clicked points without clearing the point cloud

set(gcf,'WindowButtonDownFcn',{@callbackClickA3DPoint,PointCloud}); % set the callbak
