function [vertex,face] = read_nas(filename)

% read_nas - read data from NAS file.
%
%   [vertex,face] = read_nas(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2010 JJCAO

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can not open the file.');
    return;
end
nv = 0;
nf = 0;
vertex = [];
face = []; 

while ~feof(fid)
    str = fgets(fid);
    if strcmp(str(1), 'G')
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        x = str2num(token);
        [token, str] = strtok(str);
        y = str2num(token);
        str = fgets(fid);
        [token, str] = strtok(str);
        z = str2num(str);
        nv = nv + 1;
        vertex(nv,:) = [x, y, z];
    end
    if strcmp(str(1), 'C')
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        [token, str] = strtok(str);
        fv1 = str2num(token);
        [token, str] = strtok(str);
        fv2 = str2num(token);
        fv3 = str2num(str);
        nf = nf + 1;
        face(nf,:) = [fv1, fv2, fv3]; 
    end
       
end
        
   fclose(fid);         
    
    


        

