function [A,patch_centers,patch_verts,patch_faces,verts_between_patch] = compute_face_patch_graph(faces,face_patch,verts,npatch)

% compute_face_patch_graph - compute the face patch graph of a given triangulation
%
%   [A,vertex1] = compute_face_patch_graph(faces,face_patch,verts);
%
%   'A' is the adjacency matrix of the abstract dual graph, A(i,j) is
%   number of shared verts between patch i & j. 
%       note: We've removed all A(i,j)=1.
%       (recall that this graph link together adjacent faces
%           in the triangulation).
%
%   'vertex' is optional, and if given, the position of the vertex
%   of the dual graph (contained in 'vertex1') will 
%   the centroids of the face's vertex positions.
%   
%
%   Copyright (c) 2012 JJCAO

%%
if nargin<3
    verts = [];
end
if nargin<4
    npatch = length(unique(face_patch));
end
if ~isempty(verts)
    patch_centers = zeros(npatch,3);
else
    patch_centers=[];
end

%%
npatch = length(unique(face_patch));
patch_verts = cell(npatch,1);
patch_faces = cell(npatch,1);

sprintf('compute_face_patch_graph:')
tic;
for i=1:npatch    
    patch_faces{i} = find(face_patch==i);
    ifaces = faces(patch_faces{i},:);
    ifaces = [ifaces(:,1);ifaces(:,2);ifaces(:,3)];
    patch_verts{i} = sort(unique(ifaces));
end
if ~isempty(verts)
    for i=1:npatch
        iverts = patch_verts{i};
        patch_centers(i,:) = sum(verts(iverts,:))/length(iverts);
    end
end
toc;

tic;
A = sparse(npatch, npatch);
verts_between_patch = sparse_cell(npatch, npatch);
for i=1:npatch
    iverts = patch_verts{i};   
    for j=(i+1):npatch
        jverts = patch_verts{j}; 
        tmp = intersect(iverts,jverts);
        if length(tmp)<2 % if shared verts less than 2, ingore it!
            continue;
        else
            verts_between_patch{i,j}=tmp;
            verts_between_patch{j,i}=tmp;
            A(i,j) = length(tmp);
%             i
        end
    end
end
toc;

A = max(A, A');
% verts_between_patch=max(verts_between_patch,verts_between_patch');
return

%%%%%%%%%%
figure('Name','Supervertex by NCut'); set(gcf,'color','white');
options.face_vertex_color = face_patch/npatch;
h = plot_mesh(verts, faces, options);colormap jet;hold on;view3d rot;
pverts = patch_verts{i};
h=scatter3(verts(pverts,1),verts(pverts,2),verts(pverts,3),80,'b','filled');