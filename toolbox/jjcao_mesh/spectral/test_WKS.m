%test_WKS
%
%
% Copyright (c) 2013 Junjie Cao

clc;clear all;close all;
addpath(genpath('../../'));
%% Load mesh
M.filename = 'wolf2_1k.off';
[M.verts,M.faces] = read_mesh(M.filename);
M.nverts = size(M.verts,1);

%% Laplacian eigen
nbasis = 100;
withAreaNormalization = true;
adjustL = false;
[M.eigvector,M.eigvalue, A]=compute_Laplace_eigen(M.verts,M.faces,nbasis, withAreaNormalization, adjustL);

%% WKS
tic
M.wks = WKS(M.eigvector, M.eigvalue);
toc

%% display WKS of some vertices
shading_type = 'interp';
for i = 1:ceil(size(M.wks,2)/6):size(M.wks,2)
% for i = 1:5
    figure('Name', sprintf('eigen function: %d', i));
    h=trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
        'FaceVertexCData', M.wks(:,i), 'FaceColor',shading_type); 
    axis off; axis equal; set(h, 'edgecolor', 'none'); mouse3d;
end

% %% cluster in feature space, too slow, use c++!! todo
% npatch = 20;
% A = ones(M.nverts, M.nverts);
% ind = find(A);
% [row,col] = ind2sub(size(A),ind);
% d1 = sum( (M.wks(row,50:100) - M.wks(col,50:100)).^2, 2);
% w = exp(-d1./mean(d1)); 
% A = sparse(row,col,w);
% 
% % Set diagonal elements to zeros
% A(speye(size(A))~=0)=0 ;
% 
% %% NCut
% sprintf('ncutW: ')
% tic;
% NcutDiscrete = ncutW(A,npatch);
% toc;
% 
% sprintf('cluster: ')
% tic;
% vert_patch = zeros(M.nverts,1);
% cnf = zeros(1,0);
% for j=1:npatch
%     id = find(NcutDiscrete(:,j)); 
%     if isempty(id)
%         cnf = [cnf j]
% %         warning(' cluster contains no face!');
%     else
%         vert_patch(id,:) = j-length(cnf);
%     end
% end
% 
% npatch = npatch - length(cnf);
% if ~isempty(cnf)
%     warning('some cluster contains no face!');
% end
% toc;
% 
% %% show result
% figure('Name','Supervertex by NCut'); movegui('southwest'); set(gcf,'color','white');
% options.face_vertex_color = vert_patch;
% h = plot_mesh(M.verts, M.faces, options);view3d rot;
% colormap(jet(npatch)); lighting none;
% %%
% figure('Name','Supervertex by NCut'); movegui('southeast'); set(gcf,'color','white');
% options.face_vertex_color = vert_patch;
% h = plot_mesh(M.verts, M.faces, options);view3d rot;
% colormap(colorcube(npatch)); lighting none;