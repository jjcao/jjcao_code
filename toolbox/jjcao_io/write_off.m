function write_off(filename, vertex, face, normal)

% write_off - write a mesh or a point cloud with or without normals to a OFF file
%
%   write_off(filename, vertex, face, normal);
%
%   vertex must be of size [n,3]
%   face must be of size [p,3] or empty
%   you can also use vertex color as normal
%
%   Copyright (c) 2003 Gabriel Peyr?
%   Changed by (2009) JJCAO

if nargin<4
    normal = [];
end

if size(vertex,2)~=3
    vertex=vertex';
end
if size(vertex,2)~=3
    error('vertex does not have the correct format.');
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% header
if isempty(face)
    fprintf(fid, 'NOFF\n');
    fprintf(fid, '%d 0\n', size(vertex,1));
else
    fprintf(fid, 'OFF\n');
    fprintf(fid, '%d %d 0\n', size(vertex,1), size(face,1));
end

% write the points & faces
if isempty(normal)
    fprintf(fid, '%f %f %f\n', vertex');
else
    vertex = [vertex normal];
    fprintf(fid, '%f %f %f %f %f %f\n', vertex');
end

if isempty(face)
    fclose(fid);
    return;
end

if size(face,2)~=3
    face=face';
end
if size(face,2)~=3
    error('face does not have the correct format.');
end
fprintf(fid, '3 %d %d %d\n', face'-1);

fclose(fid);