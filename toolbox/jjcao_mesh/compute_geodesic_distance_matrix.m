function [landmark_geodesic_matrix,landmark2all_geodesic_matrix] = compute_geodesic_distance_matrix(verts, faces, cids)
% compute geodesic distance matrix of landmark vertices specified by cids, of mesh (verts, faces)
% 
%
% (C) Copyright jjcao, 2012

n=length(cids);
landmark_geodesic_matrix = zeros(n,n);
sprintf('compute_geodesic_distance_matrix:')

if n == length(verts)
    landmark2all_geodesic_matrix = [];
    matlabpool(3);
    tic;
    parfor i = 1:n
        landmark_geodesic_matrix(i,:) = perform_fast_marching_mesh(verts, faces, cids(i));
    end
    toc;
    matlabpool close;
else
    landmark2all_geodesic_matrix = zeros(n,length(verts));
    for i = 1:n
        landmark2all_geodesic_matrix(i,:) = perform_fast_marching_mesh(verts, faces, cids(i));
        landmark_geodesic_matrix(i,:) = landmark2all_geodesic_matrix(i,cids);
    end
    toc;    
end

landmark_geodesic_matrix=(landmark_geodesic_matrix+landmark_geodesic_matrix')*0.5; % it is necessary to get a symmetric matrix.
