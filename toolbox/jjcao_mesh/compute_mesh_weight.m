function W = compute_mesh_weight(vertices,faces, type,options)

% compute_mesh_weight - compute a weight matrix
%
%   W = compute_mesh_weight(vertices,faces,type,options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'spring': W(i,j) = 1/d_ij where d_ij is distance between vertex
%           i and j.
%
%       'conformal' or 'dcp': W(i,j) = (cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j). (do not offer W(i,j) = cot(alpha_ij)+cot(beta_ij) anymore)
%           Refer to Computing discrete minimal surfaces and their conjugates_93, 
%           Lemma 2 of On the convergence of metric and geometric properties of polyhedral surfaces_06 and 
%           Characterizing Shape Using Conformal Factors_08.
%           Refer to Skeleton Extraction by Mesh Extraction_08, and Intrinsic Parameterizations of Surface Meshes_02.
%
%       'Mean_curvature' or 'Laplace-Beltrami': W(i,j) = (1/area_i)*(cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j), area_i is the area of vertex i's Voroni vicinity. 
%           Refer to Discrete Differential-Geometry Operators_for triangulated 2-manifolds_02
%       'Manifold-harmonic': W(i,j) = (1/sqrt(area_i*area_j))*(cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j). 
%           Refer to Spectral Geometry Processing with Manifold Harmonics_08
%
%       'mvc': W(i,j) = [tan(/_kij/2)+tan(/_jil/2)]/d_ij where /_kij and /_jil are angles at i
%
%        'Mean_curvature' or 'Laplace-Beltrami'or 'Manifold-harmonic' can
%        be built on result of 'conformal' or 'dcp'. so we do not present
%        them here.
%   If options.ring is offered, the computation of it can be avoided.
%
%   Add spring and mvc weight (JJCAO, 2009)
%   Add options.ring (JJCAO, 2009)
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
[vertices,faces] = check_face_vertex(vertices,faces);

n = max(max(faces));

if isfield(options, 'verb')
    verb = options.verb;
else
    verb = n>5000;
end

if nargin<3
    type = 'dcp';
end

switch lower(type)
    case 'combinatorial'
        W = triangulation2adjacency(faces);
    case 'distance'
        W = my_euclidean_distance_2(triangulation2adjacency(faces),vertices);
        W(W>0) = 1./W(W>0);
    case 'spring'
        W = sqrt(my_euclidean_distance(triangulation2adjacency(faces),vertices));
        W(W>0) = 1./W(W>0);
    otherwise       
        if isfield(options, 'rings')
            rings = options.rings;
        else
            rings = compute_vertex_face_ring(faces);
        end
        switch lower(type)
        case {'conformal','dcp'} % conformal laplacian              
            W = compute_mesh_weight_dcp(vertices, faces, rings,verb);                 
        case 'mvc'% mvc laplacian
            W = sparse(n,n);
            for i = 1:n
                if verb
                    progressbar(i,n);
                end
                for b = rings{i}
                    % b is a face adjacent to a
                    bf = faces(:,b);
                    % compute complementary vertices
                    if bf(1)==i
                        v = bf(2:3);
                    elseif bf(2)==i
                        v = bf([1 3]);
                    elseif bf(3)==i
                        v = bf(1:2);
                    else
                        error('Problem in face ring.');
                    end
                    j = v(1); k = v(2);
                    vi = vertices(:,i);
                    vj = vertices(:,j);
                    vk = vertices(:,k);
                    % angles
                    alpha = myangle(vi-vk,vi-vj);
                    % add weight
                    W(i,j) = W(i,j) + tan( 0.5*alpha )/sqrt(sum((vi-vj).^2));
                    W(i,k) = W(i,k) + tan( 0.5*alpha )/sqrt(sum((vi-vk).^2));
                end
            end     
        otherwise
            error('Unknown type.')
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = compute_mesh_weight_dcp(vertices, faces, rings, verb)
n = length(vertices);
W = sparse(n,n);

for i = 1:n
    if verb
        progressbar(i,n);
    end
    for b = rings{i}
        % b is a face adjacent to a
        bf = faces(:,b);
        % compute complementary vertices
        if bf(1)==i
            v = bf(2:3);
        elseif bf(2)==i
            v = bf([1 3]);
        elseif bf(3)==i
            v = bf(1:2);
        else
            error('Problem in face ring.');
        end
        j = v(1); k = v(2);
        vi = vertices(:,i);
        vj = vertices(:,j);
        vk = vertices(:,k);
        
        u = vk-vi; v = vk-vj;
        W(i,j) = W(i,j) + dot(u,v)/(eps + norm(cross(u,v)));
        u = vj-vi; v = vj-vk;
        W(i,k) = W(i,k) + dot(u,v)/(eps + norm(cross(u,v)));
%         % old
%         % angles
%         alpha = myangle(vk-vi,vk-vj);
%         beta = myangle(vj-vi,vj-vk);
%         % add weight
%         W(i,j) = W(i,j) + cot( alpha );
%         W(i,k) = W(i,k) + cot( beta );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v)
du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance_2(A,vertex)
% square euclidean distance
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2);
W = sparse(i,j,d);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance(A,vertex)
% euclidean distance
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2).^0.5;
W = sparse(i,j,d);  