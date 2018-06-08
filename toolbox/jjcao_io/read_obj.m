function [vertices,faces,veterx_normal,vertex_texture,face_texture,face_normal] = read_obj(filename)

% read_obj - load a .obj file.
%
%   [vertices,faces,veterx_normal,vertex_texture,face_texture,face_normal] = read_obj(filename);
%
%   vertices  : coordinates array
%   faces    : index of vertices
%   veterx_normal: normal of vertices
%   vertex_texture: texture of vertices
%   face_texture: index of texture
%   face_normal: index of normal
%
%   refer to: https://en.wikipedia.org/wiki/Wavefront_.obj_file#Face_elements
%   Changed by Junjie Cao @ 2012, and Shuangtao Wang @ 2018
%   Copyright (c) 2008 Gabriel Peyre

fid = fopen(filename);
if fid<0
    error(['Cannot open ' filename '.']);
end

frewind(fid);
a = fscanf(fid,'%c',1);
if strcmp(a, 'P')
    % This is the montreal neurological institute (MNI) specific ASCII facesangular mesh data structure.
    % For FreeSurfer software, a slightly different data input coding is
    % needed. It will be provided upon request.
    fscanf(fid,'%f',5);
    n_points=fscanf(fid,'%i',1);
    vertices=fscanf(fid,'%f',[3,n_points])';
    veterx_normal=fscanf(fid,'%f',[3,n_points])';
    n_faces=fscanf(fid,'%i',1);
    fscanf(fid,'%i',5+n_faces);
    faces=fscanf(fid,'%i',[3,n_faces])+1;
    fclose(fid);
    return;
end

frewind(fid);
vertices = [];
faces = [];
veterx_normal = [];
vertex_texture = [];
face_texture = [];
face_normal = [];
while 1
    s = fgetl(fid);
    if ~ischar(s) 
        break;
    end
    if ~isempty(s) && strcmp(s(1), 'f')% face
        idx = strfind(s, '/');
        if isempty(idx)
            faces(end+1,:) = sscanf(s(3:end), '%d %d %d');
        else
            s(idx) = ' ';
            tmp = sscanf(s(3:end), '%d %d %d %d %d');
            faces(end+1,:) = tmp(1:3:end);
            face_texture(end+1,:) = tmp(2:3:end);
            face_normal(end+1,:) = tmp(3:3:end);
        end
    elseif ~isempty(s) && strcmp(s(1:2), 'v ') % vertices
        vertices(end+1,:) = sscanf(s(3:end), '%f %f %f');
    elseif ~isempty(s) && strcmp(s(1:2), 'vn') % vertices normal
        veterx_normal(end+1,:) = sscanf(s(3:end), '%f %f %f');
    elseif ~isempty(s) && strcmp(s(1:2), 'vt') % texture uv
        vertex_texture(end+1,:) = sscanf(s(3:end), '%f %f');        
    end
end
fclose(fid);

