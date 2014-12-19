function sdf = read_sdf(filename, nface)
% sdf: shape diameter function
% Copyright (c) 2012 Junjie Cao
[pathstr, name, ext] = fileparts(filename);
fid = fopen([pathstr '/' name '.sdf.txt'],'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

[sdf,cnt] = fscanf(fid,'%f', nface);
if cnt~=nface
    warning('Problem in reading sdf values .');
end
fclose(fid);





