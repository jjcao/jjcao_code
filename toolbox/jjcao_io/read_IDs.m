function IDs = read_IDs(filename)

% read_IDs -- read IDs
%
% format of the file
% source mesh file name
% target mesh file name
% number of constraint pairs
% index of constraint pairs, such as 12 43
%
% IDs is a vector:
%
%   Copyright (c) 2010 JJCAO

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);
a=strtok(str); nID = str2num(a);
[IDs,cnt] = fscanf(fid,'%d', nID );
if (cnt~=nID)
    warning('Problem in reading ID. Number of ids is not right, maybe.');
end
fclose(fid);
return;