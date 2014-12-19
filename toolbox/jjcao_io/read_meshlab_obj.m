function [verts,faces,normals] = read_meshlab_obj(inputmesh)

% read_meshlab_obj - load a .obj point cloud file in meshlab format.
% NB: not each vertex has a normal!
%
%   [verts,faces,normals] = read_meshlab_obj(inputmesh);
%
% faces  : list of triangle elements
% verts  : node coordinatates
%
% Copyright (c) 2009 JJCAO
% Update history 2009.10

fid=fopen(inputmesh);
while true
    str = fgetl(fid);
    tmp = strfind(str, 'Vertices');
    if ~isempty(tmp)
        nverts = sscanf(str,'%*s %*s %d');
        str = fgetl(fid);
        nfaces = sscanf(str,'%*s %*s %d');
        break;
    end
end
verts = zeros(nverts, 3);
faces = [];
normals = zeros(nverts, 3);

while str
    str = fgetl(fid);
    tmp = strfind(str, '#');
    if isempty(tmp)
        break;
    end
end

i = 0;
while str
    i = i+1;
    tmp = strfind(str, '#');
    if ~isempty(tmp)
        break;
    end
    
    tmp = strfind(str, 'vn');
    if isempty(tmp)
        verts(i,:) = sscanf(str, '%*s %f %f %f');
        str = fgetl(fid);
        continue;
    else
        normals(i,:) = sscanf(str, '%*s %f %f %f');
        str = fgetl(fid);     
        verts(i,:) = sscanf(str, '%*s %f %f %f');
        str = fgetl(fid);
        continue;
    end    
end

fclose(fid);