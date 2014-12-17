function B = fitPlane(XYZ)
% FITPLANE - solves coefficients of plane fitted to 3 or more points
%
% Usage:   B = fitPlane(XYZ)
%
% Where:   XYZ - 3xNpts array of xyz coordinates to fit plane to.   
%                If Npts is greater than 3 a least squares solution 
%                is generated.
%
% Returns: B   - 4x1 array of plane coefficients in the form
%                b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0
%                The magnitude of B is 1.
%
% Uses command svd
% See also: RANSACFITPLANE
%
% Copyright (c) 2003-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
  
  [rows,npts] = size(XYZ); 
  
  if rows ~=3
    error('data is not 3D');
  end
  
  if npts < 3
    error('too few points to fit plane');
  end
  
  % Set up constraint equations of the form  AB = 0,
  % where B is a column vector of the plane coefficients
  % in the form   b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.
  
  A = [XYZ' ones(npts,1)]; % Build constraint matrix
  
  if npts == 3             % Pad A with zeros
    A = [A; zeros(1,4)]; 
  end

  [u d v] = svd(A);        % Singular value decomposition.
  B = v(:,4);              % Solution is last column of v.
  
function [] = test_fitPlane()
[x,y]=meshgrid(linspace(0,10,20),linspace(0,10,20));
b=[0,0,1,1];
z=-(b(1)*x+b(2)*y+b(4))/b(3);
x=x(:); y=y(:); z=z(:);
%z=z+(randn(length(z),1));
B= fitPlane([x,y,z]'); 
[X,Y]=meshgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
Z=-(B(1)*X+B(2)*Y+B(4))/B(3);
plot3(x,y,z,'r.'); hold on; grid on;
surf(X,Y,Z,'FaceColor','g'); alpha(0.5);
title(['B=' num2str(B')]);view3d rot;