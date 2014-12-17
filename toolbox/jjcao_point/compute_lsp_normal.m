function [ normal ] = compute_lsp_normal( pointset )
%FINDMSPN Summary of this function goes here
%   Detailed explanation goes here
% 找输入顶点集合的最小二乘平面的单位法向
% pointset 是一个n-by-3矩阵
x = pointset(:,1);
y = pointset(:,2);
z = pointset(:,3);

[A,B,C]=compute_plane_fit(x,y,z);
normal = [A,B,C]/sqrt(A*A+B*B+C*C);

end

