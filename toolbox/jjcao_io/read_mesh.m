function [vertex,face,normal, uv, sphparam] = read_mesh(file)

% read_mesh - read data from OFF, PLY, SMF or WRL file.
%
%   [vertex,face] = read_mesh(filename);
%   [vertex,face] = read_mesh;      % open a dialog box
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Supported file extensions are: .off, .ply, .wrl, .obj, ,ted, nas .m, .gim.
%
%   Changed by jjcao, 2011, 2012
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin==0
    [f, pathname] = uigetfile({'*.off;*.obj,*.ply;*.wrl;','*.off,*.ply,*.wrl,*.smf,*.nas Files'},'Pick a file');
    file = [pathname,f];
end

i = strfind(file,'.');
ext = file(i(length(i))+1:end);

switch lower(ext)
    case 'off'
        [vertex,face,normal] = read_off(file); % adapted by jjcao
    case 'ply'
        [vertex,face] = read_ply(file);
    case 'smf'
        [vertex,face] = read_smf(file);
    case 'wrl'
        [vertex,face] = read_wrl(file);
    case 'obj'
        [vertex,face,normal,uv] = read_obj(file);
    case 'tet'
        [vertex,face] = read_tet(file);        
    case 'm'
        if isfield(options, 'type')
            type = options.type;
        else
            type = 'gim';
        end
        [vertex,face,normal, uv, sphparam] = read_mfile(file, type);
    case 'gim'
        sub_sample = 1;
        [M,Normal] = load_gim(name, options);
        [vertex,face] = convert_gim2mesh(M, sub_sample);
        normal = convert_gim2mesh(Normal, sub_sample);
    case 'nas' % added by jjcao
        [vertex,face] = read_nas(file);        
    otherwise
        error('Unknown extension.');
end