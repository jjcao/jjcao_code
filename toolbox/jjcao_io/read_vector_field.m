function vf = read_vector_field(filename)
%format:
%24 3
%0.999987 0.00013041 -0.00361332 0.00013041 0.998699 0.036044 -0.00361332 0.036044 0.00131395
%0.781284 -0.0618677 0.40872 -0.0618677 0.9825 0.115613 0.40872 0.115613 0.236217
%...
%
%   Copyright (c) 2008 JJCAO

%% open file
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

%% read num of vectors in a vector field and num of vector fields
str = fgets(fid);
[a,str] = strtok(str); nvect = str2num(a);%num of vectors in a vector field
[a,str] = strtok(str); nvf = str2num(a);%num of vector fields

%% read vector fields
[vf,cnt] = fscanf(fid,'%f %f %f',3*nvect*nvf);
if cnt~=3*nvect*nvf
    warning('Problem in reading vector field.');
end
vf = reshape(vf, 3*nvf,nvect);
vf=vf';

%% close file and return
fclose(fid);
return;