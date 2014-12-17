function vf = read_constraints(filename)
%format: 
% number
% id x y z
%example:
% 24
% 11 0.999987 0.00013041 -0.00361332
% 3 0.781284 -0.0618677 0.40872
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
[a,str] = strtok(str); num = str2num(a);%num of vectors in a vector field

%% read vector fields
[vf,cnt] = fscanf(fid,'%d %f %f %f',4*num);
if cnt~=4*num
    warning('Problem in reading vector field.');
end
vf = reshape(vf, 4,num);
vf=vf';

%% close file and return
fclose(fid);
return;