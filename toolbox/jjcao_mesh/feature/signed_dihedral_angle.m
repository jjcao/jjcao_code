function [sda, n1, n2] = signed_dihedral_angle(verts)
% verts: 4 verts arranged anti-clockwise; 1st - 3th verts => triangle 1 anti-clockwise; 
%                                         2nd - 4th => triangle 2
%        v1 is the base point.
% sda: signed dihedral angle

% jjcao @ 2014

v1 = verts(1,:);
v2 = verts(2,:);
v3 = verts(3,:);
v4 = verts(4,:);

n1 = cross(v2 - v1, v3 - v2);
n1 = n1./(dot(n1,n1));
n2 = cross(v3 - v1, v4 - v3);
n2 = n2./(dot(n2,n2));

e31 = v1 - v3;

tmp = dot(cross(n1,n2), e31);
sda = acos(dot(n1,n2)) * tmp/abs(tmp);