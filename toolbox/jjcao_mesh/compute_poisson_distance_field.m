function [Poisson_dist A]= compute_poisson_distance_field(vertices,faces,constraint_values,delta_values,options)
%
%  constraint_values: col 1 is id, col 2 is constraint values
%  options.weights
%  Copyright (c) 2008 JJCAO

method = getoptions(options, 'method', 'hard');

if size(vertices,2)>size(vertices,1)
    vertices = vertices';
end
if size(faces,2)>size(faces,1)
    faces = faces';
end

switch lower(method)    
    case 'hard'
        [Poisson_dist A] = compute_poisson_distance_field_hard(vertices,faces,constraint_values,delta_values,options);
    case {'soft','bi-harmonic'}
        [Poisson_dist A] = compute_poisson_distance_field_soft(vertices,faces,constraint_values,delta_values,options);
    otherwise
        error('Unknown method.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pdf L]= compute_poisson_distance_field_hard(vertices, faces, constraint_values, delta_values,options)
%
%
%   Copyright (c) 2008 JJCAO
global logstr;
nverts = size(vertices,1);
ncvs = size(constraint_values,1);
if isfield(options, 'type')
    type = options.type;
else
    type = 'conformal';
end

% setup the system
oldTime=cputime;
L = compute_mesh_laplacian(vertices,faces, type,options);
b=delta_values;
for i = 1:ncvs
    b(constraint_values(i,1),1)=constraint_values(i,2);
    L(constraint_values(i,1),:) = zeros(1,nverts);
    L(constraint_values(i,1),constraint_values(i,1)) = 1;
end
[m n] = size(L);
s=sprintf('compute laplacian matrix: %d*%d in %f seconds', m, n, cputime - oldTime );
logstr = sprintf('%s%s\n',logstr,s);

% solve
oldTime=cputime;
pdf = L\b;
s=sprintf('solve in %f seconds', cputime - oldTime );
logstr = sprintf('%s%s\n',logstr,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pdf L]= compute_poisson_distance_field_soft(vertices, faces, constraint_values, delta_values,options)
%  
%  slow a little than hard
%  Copyright (c) 2009 JJCAO
global logstr;
nverts = size(vertices,1);
ncvs = size(constraint_values,1);
if isfield(options, 'type')
    type = options.type;
else
    type = 'conformal';
end

% setup the system
oldTime=cputime;
L = compute_mesh_laplacian(vertices,faces, type,options);
% set constraints
n = size(vertices,1);
D = zeros(n,1);
D(constraint_values(:,1)) = options.c_weights;
P = sparse(1:n, 1:n, D, n,n);
delta_values(constraint_values(:,1)) = constraint_values(:,2);
b = P*delta_values;

s=sprintf('compute laplacian matrix: %d*%d in %f seconds', n, n, cputime - oldTime );
logstr = sprintf('%s%s\n',logstr,s);

% solve
if strcmpi(options.method,'bi-harmonic')
    oldTime=cputime;
    L = L'*L;
    s=sprintf('multiply in %f seconds', cputime - oldTime );
    logstr = sprintf('%s%s\n',logstr,s);
end

oldTime=cputime;
pdf = (L+P)\b;
s=sprintf('solve in %f seconds', cputime - oldTime );
logstr = sprintf('%s%s\n',logstr,s);

% ////////////////////////////////////////////////////////////////////////
% % solve
% LtL=L'*L;
% pdf = LtL\(L'*b);
% %pdf = cholmod2 (LtL, L'*b); %比 LtL\(L'*b);慢一些
% %pdf = L\b;%10多倍的慢！！！！！！！！！！！！！！！！

