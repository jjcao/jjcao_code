function write_point_obj(filename, verts)

% write_point_obj - write a point cloud to an OBJ file
%
%   write_point_obj(filename, verts)
%
%   verts must be of size [n,3]
%
%   Copyright (c) 2009 JJCAO
% filename = '../data/1338_Galaad_clean.obj';
% [M.verts,M.faces] = read_mesh(filename);
% write_point_obj('test1.obj',M.verts);

nverts = size(verts, 1);

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

fprintf(fid, '# Studio\n');
fprintf(fid, 'g Point_Model_1\n');

% vertex position
for i = 1:nverts
    fprintf(fid, 'v %f %f %f\np %d\n', verts(i,:), i);
end
fprintf(fid, '# end of file\n');

fclose(fid);