function [] = write_pwn(filename, verts, normals)

% write_pwn - write a mesh to a pwn file for MPU reconstruction
%
%   write_pwn(filename, vertex, normals);
%
%   vertex must be of size [n,3]
%   normals must be of size [p,3]
%
% % example
% filename = '../data/fertility_4n.off';
% [M.verts, M.faces, M.normals] = read_off(filename);
% write_pwn([filename(1:(end-3)) 'pwn'], M.verts, M.normals);
%
%   Copyright (c) 2009 JJCAO

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

% header
fprintf(fid, '%d\n', size(verts,1));

% write the points & normals
fprintf(fid, '%f %f %f\n', verts');
fprintf(fid, '%f %f %f\n', normals');

fclose(fid);