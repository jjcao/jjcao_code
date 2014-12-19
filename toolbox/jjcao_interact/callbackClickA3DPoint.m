% callbackClickA3DPoint.m
% by Babak Taati
% http://qlink.queensu.ca/~3bt1/
% Queen's University
% May 4, 2005
%
%
% this is the callback function for ClickA3DPoint
%   

function callbackClickA3DPoint(src,eventdata, PointCloud);

Point = get(gca,'CurrentPoint');        % get the mouse click position

CamPos = get(gca,'CameraPosition');     % camera position
CamTgt = get(gca,'CameraTarget');       % where the camera is pointing to

CamDir = CamPos - CamTgt;               % Camera direction

CamUpVect = get(gca,'CameraUpVector');  % camera 'up' vector


% build a orthonormal frame based on the viewing direction and the up vector (the "view frame")
Z_axis = CamDir/norm(CamDir);    
UP_axis = CamUpVect/norm(CamUpVect); 
X_axis = cross(UP_axis,Z_axis);
Y_axis = cross(Z_axis, X_axis);


Rot = [ X_axis ; Y_axis ; Z_axis ];     % view rotation 

RotatedPointCloud = Rot * PointCloud;   % the point cloud represented in the view frame
RotatedPointFront = Rot * Point' ;      % the clicked point represented in the view frame


% --- find the nearest neighbour to the clicked point 
%
% project the 3D point cloud onto a plane (i.e. ignore the z coordinate) 
% and find the point that is the closest to the clicked point in that 2D plane.

% % older version
% diff = [ RotatedPointCloud(1,:) - RotatedPointFront(1)   ;  RotatedPointCloud(2,:) - RotatedPointFront(2) ]; % the difference of all the 2D points from the clicked point (2D)
% AllDistances = RowNorm(diff');                  % the distance of all the 2D points from the clicked point (2D)
% [dist, PointCloudIndex] = min(AllDistances);	% find the point that is closest to the clicked point (2D)

% update (June 2, 2008): this does the same thing
PointCloudIndex = dsearchn( RotatedPointCloud(1:2,:)' , RotatedPointFront(1:2) );

% ---- delete the old nearest neighbour & display the new nearest neighbour

h   = findobj(gca,'Tag','pt'); % try to find the old point
SelectedPoint = PointCloud(:,PointCloudIndex); 

if isempty(h) % Check if it's the first click (i.e. no need to delete the previous point)
    
    h = plot3(SelectedPoint(1,:), SelectedPoint(2,:), SelectedPoint(3,:), 'r.', 'MarkerSize', 20); % highlight the selected point
    set(h,'Tag','pt');   % set its Tag property for later use   

else    % if Not the first click

    delete(h);      % delete the previously selected point
    
    h = plot3(SelectedPoint(1,:), SelectedPoint(2,:), SelectedPoint(3,:), 'r.', 'MarkerSize', 20);  % highlight the newly selected point
    set(h,'Tag','pt');  % set its Tag property for later use

end

fprintf('you clicked on point number %d\n', PointCloudIndex);
