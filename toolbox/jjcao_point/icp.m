function [R, t, corr, data2] = icp(data1, data2, tol)
% [R, t, corr, data2] = icp(data1, data2, tol)
%
% This is an implementation of the Iterative Closest Point (ICP) algorithm.
% The function takes two data sets and registers data2 with data1. It is
% assumed that data1 and data2 are in approximation registration. The code
% iterates till no more correspondences can be found.
%
% Arguments: data1 - 3 x n matrix of the x, y and z coordinates of data set 1
%            data2 - 3 x m matrix of the x, y and z coordinates of data set 2
%            tol   - the tolerance distance for establishing closest point
%                     correspondences
%
% Returns: R - 3 x 3 accumulative rotation matrix used to register data2
%          t - 3 x 1 accumulative translation vector used to register data2
%          corr - p x 2 matrix of the index no.s of the corresponding points of
%                 data1 and data2
%          data2 - 3 x m matrix of the registered data2 
%
% Copyright : This code is written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%              Computer Science, The University of Western Australia. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion this copyright information.

c1 = 0;
c2 = 1;
R = eye(3);
t = zeros(3,1);
tri = delaunayn(data1');
while c2 > c1    
    c1 = c2;
    [corr, D] = dsearchn(data1', tri, data2');
    corr(:,2) = [1 : length(corr)]';
    ii = find(D > tol);
    corr(ii,:) = [];
    [R1, t1] = reg(data1, data2, corr);
    data2 = R1*data2;
    data2 = [data2(1,:)+t1(1); data2(2,:)+t1(2); data2(3,:)+t1(3)];
    R = R1*R;
    t = R1*t + t1;    
    c2 = length(corr);    
end

%-----------------------------------------------------------------
function [R1, t1] = reg(data1, data2, corr)
n = length(corr); 
M = data1(:,corr(:,1)); 
mm = mean(M,2);
S = data2(:,corr(:,2));
ms = mean(S,2); 
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K = Sshifted*Mshifted';
K = K/n;
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(3);
    B(3,3) = det(V*U');
    R1 = V*B*U';
end
t1 = mm - R1*ms;