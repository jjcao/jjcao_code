function write_mesh(filename, vertex, face, normal)

% write_mesh - read data to OFF, PLY, SMF or WRL filename.
%
%   write_off(filename, vertex, face);
%   write_off([], vertex, face); % open a dialog box
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%   'normal' is a 'nb.vert x 3' array specifying vertex normals.
%
%   Add normal by jjcao, 2009, 2012
%   Copyright (c) 2005 Gabriel Peyr?

if isempty(filename)
    currenttime=clock;
    defaultname=strcat(num2str(currenttime(1)),'_',num2str(currenttime(2)),'_',...
        num2str(currenttime(3)),'(',num2str(currenttime(4)),'_',...
        num2str(currenttime(5)),'_',num2str(currenttime(6)),')','.off');

    filterspec = {'*.off', 'off'; ...
                  '*.obj', 'obj'; ...
                  '*.xyz', 'xyz';...
                  '*.ply', 'ply';...
                  '*.nas', 'nas';};
    [filename,pathname,filterindex] = uiputfile(filterspec,'Save as',defaultname);
    if isequal(filename,0) || isequal(pathname,0)
       disp('User selected Cancel');
       return;
    end
    filename = fullfile(pathname,filename);
end

ext = filename(end-2:end);
ext = lower(ext);
if isempty(ext)
    ext = filterspec{filterindex,2};
end

if strcmp(ext, 'off')
    write_off(filename, vertex, face, normal);
elseif strcmp(ext, 'xyz')    
    write_xyz(filename, vertex, normal);
elseif strcmp(ext, 'ply')
    write_ply(filename, vertex, face);
elseif strcmp(ext, 'smf')
    write_smf(filename, vertex, face);
elseif strcmp(ext, 'wrl')
    write_wrl(filename, vertex, face);
elseif strcmp(ext, 'obj')
    write_point_obj(filename, vertex);
else
    error('Unknown extension.');    
end