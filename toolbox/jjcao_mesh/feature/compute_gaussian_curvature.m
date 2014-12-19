function Cgauss = compute_gaussian_curvature(vertices, faces, options)
% compute_gaussian_curvature - compute Gaussian curvature of a mesh
%
%   Cgauss = compute_gaussian_curvature(vertices, faces, options);
%
%   options.type = 'angle_defect'
%   options.type = 'normal_cycles'
%   options.boundary: id of all boundary vertices
%  
%   Copyright (c) 2009 JJCAO

options.null = 0;
if ~isfield(options, 'type')
    options.type = 'angle_defect';
end

switch lower(options.type)
    case {'normal_cycles','normal cycles', 'normal_cycle', 'normal cycle'}
        [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertices,faces,options);        
    case {'angle_defect', 'angle defect'}
        Cgauss = compute_angle_defect(vertices, faces, options);
    otherwise
        disp('Unknown method.')
end