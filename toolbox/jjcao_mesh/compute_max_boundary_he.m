function boundary = compute_max_boundary_he(Mhe, options)

%   If options.type = 'max_vnum' then 
%       choose boundary with most vertices as boundary
%   Copyright (c) 2009 JJCAO

if isfield(options, 'type')
    type = options.type;
else
    type = 'max_vnum';
end

if STRCMP(type, 'max_vnum')
    boundary = Mhe.boundary_vertices;
    if iscell(boundary)
        tmp = size(boundary);
        for i = 1:length(boundary)
            tmp(i) = length(boundary{i});
        end
        [m I]  = max(tmp);
        boundary = boundary(I);
    end
    boundary = boundary';
else
    error('Does not work with other type, still!');    
end