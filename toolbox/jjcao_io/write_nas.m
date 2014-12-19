function write_nas(filename, vertex, face);

%   write_nas - write a mesh to a NAS file
%
%   write_nas(filename, vertex, face);
%
%   vertex must be of size [n,3]
%   face must be of size [p,3]
%   Copyright (c) 2003 Gabriel Peyr?

if size(vertex,2)~=3
    vertex=vertex';
end
if size(vertex,2)~=3
    error('vertex does not have the correct format.');
end

if size(face,2)~=3
    face=face';
end
if size(face,2)~=3
    error('face does not have the correct format.');
end

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

for i=1:size(vertex,1)
    pos = vertex(i,:);
 fprintf(fid, 'GRID   *%16d               0%16f%16f  N%-5d\n*N%-6d%16f\n',i, pos(1), pos(2),i,i, pos(3));
%     fprintf(fid, 'GRID   *%16d%160\t\t%f\t\t%f  N%d    \n*N%d\t\t%f\n',i, pos(1), pos(2),i,i, pos(3));
end

for i=1:size(face,1)
    f = face(i,:);    % index start at ?
    fprintf(fid, 'CTRIA3  %8d       0%8d%8d%8d\n',i, f(1), f(2), f(3));   
%     CTRIA3       230       0      90     111     112
end

fclose(fid);