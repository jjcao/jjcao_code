function fid = compute_least_square_system(A, b, constraint_id, constraint_value,options)
% solve Ax=b in least square way, with Dirichlet boundary conditions (constraints)
%options.solver  1: simplest way; 2: decompose matrix LU; 3: decompose matrix naive
% Copyright (c) 2013 Junjie Cao, Zhongguang Yu
options.null = 0;
method = getoptions(options, 'method', 'hard'); % 'hard', 'soft'
solver = getoptions(options, 'solver', 1); % 1: simplest way;
switch lower(method) 
    case 'hard'
        fid = compute_lss_hard(A, b, constraint_id, constraint_value,solver);
    case {'soft'}
        error('to be supported!');
        %fid = compute_lss_soft(vertices,faces,constraint_id,delta_coords,s_weights,c_weights,options);
    otherwise
        error('Unknown method.');
end

function fid = compute_lss_hard(A, b, constraint_id, constraint_value,solver)
% jjcao
switch solver
    case 1 % worked both for b with 1 column or multi-column
        b(constraint_id,:) = constraint_value;
        A(constraint_id,:) = 0;
        n = size(A,1);
        A = A + sparse(constraint_id, constraint_id, 1, n, n);
        
        fid = A\b;
    case 2 % Kim G, Xing E P, Fei-Fei L, et al. Distributed cosegmentation via submodular optimization on anisotropic diffusion (ICCV), 2011
        uind = 1:size(A,1);
        uind(constraint_id) = [];
        laplace_mat = A;
        A = laplace_mat(uind, uind);
        tmp = laplace_mat(uind, constraint_id)*constraint_value';
        b = -tmp+b(uind);
        
        if det(laplace_mat(uind,uind))==0 %|| c < eps
            % incomplete LU decomposition
            % Parameters for 'ilu'
            ilu_setup=[];
            ilu_setup.type = 'ilutp';
            [L U P] = ilu(sparse(laplace_mat(uind,uind)), ilu_setup) ;
            [p, foo] = find(P); p(p) = 1:length(p) ;
            fid(uind) = U\(L\(b(p,:)));
        else
            fid(uind)= inv(A)*b;
        end
        
        fid(constraint_id) = constraint_value;
        fid = fid';
    case 3 % worked just for b with 1 column
        uind = 1:size(A,1);
        uind(constraint_id) = [];
        A_uu = A(uind,uind);
        b_L = -A(uind, constraint_id)*constraint_value'+b(uind);
        
        fid(uind)= A_uu\b_L;
        fid(constraint_id) = constraint_value;
        fid = fid';        
    otherwise
        error('Unknown method.');
end
