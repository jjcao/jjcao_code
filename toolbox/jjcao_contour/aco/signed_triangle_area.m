%
% This function returns the *signed* area of a triangle
% ABC in 2D. Crosss product is used.
%
% a = signed_triangle_area(A, B, C)
%
% Assume that A, B, and C are given in counterclockwise order.
%
% ------------------------------------------------------
% (C) Richard (Hao) Zhang (2005)
%
function a = signed_triangle_area(A, B, C)
u = B - A;
v = C - A;

cp = cross([u 0], [v 0]);
a = sign(cp(1,3)) * 0.5 * norm(cp);
