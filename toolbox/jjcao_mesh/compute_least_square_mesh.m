function [newVertices, A] = compute_least_square_mesh(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options)
% compute_least_square_mesh - compute a least square mesh
%
% [newVertices A] =
% compute_least_square_mesh(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options)
%
% options.method can be: 'hard', 'soft', 'bi-harmonic'
% A is the final matrix (Ax=b);
%
% modified by jjcao, 2012
% Copyright (c) 2009 JJCAO

method = getoptions(options, 'method', 'bi-harmonic');
if ~isfield(options, 'rings')
    options.rings = compute_vertex_face_ring(faces);
end

if ~isfield(options, 'normalize')
    options.normalize = 0; 
end
if ~isfield(options, 'symmetrize')
    options.symmetrize = 1; 
end
if ~isfield(options, 'constraint_pos')
    options.constraint_pos = vertices(constraint_id,:);
end

switch lower(method) 
    case {'bi-harmonic'}
        [newVertices A] = compute_lsm_bh(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options);
    case 'hard'
        [newVertices A] = compute_lsm_hard(vertices,faces,constraint_id,delta_coords,options);
    case {'soft'}
        [newVertices A] = compute_lsm_soft(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options);
    otherwise
        error('Unknown method.');
end

function [newVertices, A] = compute_lsm_bh(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options)
% compute_ls_soft - compute a least square mesh according to "Least-squares
% Meshes" by Olga Sorkine and Daniel Cohen-Or, 2004.
%
% A is the final matrix (Ax=b);
%
% Copyright (c) 2009 JJCAO
nverts = length(vertices);
ncvs = length(constraint_id);

% Build matrix ////////////////////////////////////////////////////////////
sprintf('Build Laplacian matrix, soft:')
tic;
type = getoptions(options, 'type', 'dcp');%  'combinatorial', 'distance', spring, 'Laplace-Beltrami', mvc, 
S = sparse(1:nverts,1:nverts,s_weights, nverts, nverts);
L = compute_mesh_laplacian(vertices,faces,type,options);

B = zeros(nverts+ncvs,size(vertices,2));
for coord = 1:size(vertices,2)
    % compute right hand side
    b = S*delta_coords(:,coord);
    b(nverts+1:nverts+ncvs,1)=c_weights.*options.constraint_pos(:,coord);
    B(:,coord) = b;
end
toc;

% solve ///////////////////////////////////////////////////////////////////
C = sparse(1:ncvs, constraint_id', c_weights, ncvs,nverts);
A = S*L;
sprintf('solve Laplacian matrix, bi-harmonic')
tic
newVertices = (A'*A + C'*C)\([A;C]'*B);
A = [A;C];
toc


function [newVertices, A] = compute_lsm_soft(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options)
% compute_ls_soft - compute a least square mesh according to "Dynamic Harmonic Fields for Surface Processing" 2009.
%
% A is the final matrix (Ax=b);
%
% Copyright (c) 2009 JJCAO
nverts = length(vertices);

% Build matrix ////////////////////////////////////////////////////////////
sprintf('Build Laplacian matrix, soft:')
tic;
type = getoptions(options, 'type', 'dcp');%  'combinatorial', 'distance', spring, 'Laplace-Beltrami', mvc, 
S = sparse(1:nverts,1:nverts,s_weights, nverts, nverts);
L = compute_mesh_laplacian(vertices,faces,type,options);

D = zeros(nverts,1);
D(constraint_id) = c_weights;
P = sparse(1:nverts,1:nverts, D, nverts,nverts);

B = zeros(nverts,size(vertices,2));
for coord = 1:size(vertices,2)
    % compute right hand side
    b = delta_coords(:,coord);
    b(constraint_id)= options.constraint_pos(:,coord);  
    B(:,coord) = (S+P)*b;
end
toc;

% solve ///////////////////////////////////////////////////////////////////
sprintf('solve Laplacian matrix, soft:') 
tic;
A=(S*L+P);
newVertices = A\B;
toc;

function [newVertices L] = compute_lsm_hard(vertices,faces,constraint_id,delta_coords,options)
type = getoptions(options, 'type', 'dcp');%  'combinatorial', 'distance', spring, 'Laplace-Beltrami', mvc,

% Build matrix ////////////////////////////////////////////////////////////
sprintf('Build Laplacian matrix, hard:')
tic;
L=compute_mesh_laplacian(vertices,faces,type,options);
B=delta_coords;
n = length(vertices);
B(constraint_id,:) = options.constraint_pos;
L(constraint_id,:) = 0;
L = L + sparse(constraint_id, constraint_id, 1, n, n);
toc;

% solve ///////////////////////////////////////////////////////////////////
sprintf('solve Laplacian matrix, hard:')
tic
newVertices = L\B;
toc

