function A =compute_patch_iso_similarity(patch_adjancy,verts, faces, patch_verts,patch_faces,face_patch)
% compute patch similarity by patch's developbility and patch's boundary normal
% using formular (2) of sgp04_Iso-charts
%
% (C) Copyright jjcao, 2012

DEBUG=1;
if DEBUG
    figure('Name','Segment by SMO');movegui('southwest'); set(gcf,'color','white');
end
%% compute formula (4) of sgp04_iso-charts
% matlabpool open 7
tic;
A = sparse( size(patch_adjancy,1), size(patch_adjancy,2) );
for i=1:length(patch_adjancy)
%     i
    %% isomap of patch i
    iverts = patch_verts{i};
    DI = zeros(length(iverts));    
    ifaces = patch_faces{i};
    ifaces = faces(ifaces,:);
    for j=1:length(iverts)
        start_points = iverts(j);        
        tmp = perform_fast_marching_mesh(verts, ifaces, start_points);
%         tmp = perform_fast_marching_mesh(verts, faces(ifaces,:), start_points); 
        DI(j,:) = tmp(iverts);
    end
    DI=(DI+DI')*0.5; % it is necessary to get a symmetric matrix.
    [vec val] = compute_isomap(DI, 2);
    
    %% isomap of patch i & j
    npatchid = find(patch_adjancy(i,:)>0);
    for j=1:length(npatchid)
        ifaces = [ifaces; faces(patch_faces{npatchid(j)},:)];
    end
    for j=1:length(npatchid)
%         j        
        jverts = patch_verts{npatchid(j)};
        DJ = zeros(length(jverts));
        DJI = zeros(length(jverts),length(iverts));
        for k=1:length(jverts)
            start_points = jverts(k);  
            tmp = perform_fast_marching_mesh(verts, ifaces, start_points); 
            DJI(k,:) = tmp(iverts);
            DJ(k,:) = tmp(jverts);
        end        
        nvec = compute_mds_interpolation(DI, vec, val, DJI);    
        D = [[DI DJI'];[DJI, DJ]];
        D=(D+D')*0.5;
        A(i,npatchid(j)) = sum(compute_geodesic_distortion([vec;nvec], D))
%         A(i,npatchid(j)) = sum(compute_geodesic_distortion(nvec, DJ))
%        %%
        if DEBUG    
            options.face_vertex_color = zeros(size(face_patch));
            options.face_vertex_color(face_patch==i)=10;
            options.face_vertex_color(face_patch==npatchid(j))=20;
            h = plot_mesh(verts, faces, options);
            colormap(jet);view3d rot;lighting none;        
            delete(h);        
        end
    end    
end
% matlabpool close
toc;
