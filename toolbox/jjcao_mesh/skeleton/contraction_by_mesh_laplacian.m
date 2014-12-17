function [cverts] = contraction_by_mesh_laplacian(verts, faces)
% 
% just a init test
%
%
% refer to Skeleton Extraction by Mesh Contraction 08
% mesh could be simplified by MeshLab/Filter/Clustering Decimation, if it
% is too larged.
% 
% inputs:
%   
%
% outputs:
%   cpts: contracted vertices
%
%
% by: JJCAO @ 2013
%

% profile on;

%##########################################################################
%% setting
%##########################################################################
% visual debug conditions
SHOW_CONTRACTION_PROGRESS = true;
Laplace_type = 'conformal';%conformal%combinatorial%spring%mvc
options.normalize = 0; options.symmetrize = 1;
% setting
tc = getoptions(options, 'tc', 0.01); % contract Termination Conditions for total area ratio 0.01 
iterate_time = getoptions(options, 'iterate_time', 20); 

%% init iteration
if SHOW_CONTRACTION_PROGRESS
    figure('name','CONTRACTION'); set(gcf,'color','white');hold on;
    h=trisurf(faces, verts(:,1), verts(:,2), verts(:,3), 'FaceColor', 'cyan',  'faceAlpha', 0.99); axis off; axis equal; mouse3d    
else
    h = [];
end
t = 1; 
cverts = verts;
while t<iterate_time
    [L, D, W] = compute_mesh_laplacian(cverts,faces, Laplace_type, options);%L = -compute_mesh_laplacian(cverts, faces,Laplace_type,options);
    A = D-0.99*W;
    b = cverts;
    cverts = A\b; 
    
    t = t+1
    
    if SHOW_CONTRACTION_PROGRESS
        % 显示前后点云     
        if ~isempty(h), delete(h),end;
        h=trisurf(faces, cverts(:,1), cverts(:,2), cverts(:,3), 'FaceColor', 'cyan',  'faceAlpha', 0.99); axis off; axis equal; mouse3d    
        title(['iterate ',num2str(t),' time(s)']); drawnow; 
    end
end

% profile off;
% profile viewer;
