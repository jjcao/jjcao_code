function fid = compute_mesh_harmonic_field(verts, faces, constraint_id, constraint_value,type, options)
% 
% Copyright (c) 2013 Junjie Cao
options.null = 0;
L = compute_mesh_laplacian(verts,faces,type,options);

options.solver = getoptions(options, 'solver', 1);
options.method = 'hard';% 'hard', 'soft'
b = zeros(size(L,1),1);
fid = compute_least_square_system(L, b, constraint_id, constraint_value,options);