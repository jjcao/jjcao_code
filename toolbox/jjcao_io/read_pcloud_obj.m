function [verts, normals] = read_pcloud_obj(filename)

% read_point_obj - load a .obj point cloud file.
%
%   [verts,normals] = read_point_obj(inputmesh);
%
% verts  : node coordinatates
%
%(C) JJCAO
% Update history 2009.10.13

fid=fopen(filename);
tline = fgetl(fid);
if strcmp(tline,'# Studio')
    [verts, normals] = read_studio_obj(fid);
else
    [verts, normals] = read_adrea_obj(fid);
end

fclose(fid);

function [verts, normals] = read_studio_obj(fid)
fseek(fid, -20, 'eof');
nverts = fscanf(fid,'%f');
verts = zeros(nverts,3);
normals = [];

frewind(fid);
tline = fgetl(fid);
tline = fgetl(fid);
for i = 1:nverts
    fseek(fid, 2, 'cof');
    verts(i,:) = fscanf(fid,'%f %f %f');
    fgetl(fid);    
end

function [verts, normals] = read_adrea_obj(fid)
i = 0;
frewind(fid);
while ~feof(fid)
    tmp = fgetl(fid);
    i = i+1;
end
frewind(fid);
if strcmp('vn',tmp(1:2))
    nverts = i/2;
else
    nverts = i;
end

[A,cnt] = fscanf(fid,'%s %f %f %f', 4*nverts);
A = reshape(A, 4, length(A)/4);
A = A(2:end,:);
verts = A';
 
if nverts == i/2
    [A,cnt] = fscanf(fid,'%s %f %f %f', 5*nverts);
    A = reshape(A, 5, length(A)/5);
    A = A(3:end,:);
    normals = A';
else
    normals = [];
end