function patch_curvature_hist =compute_patch_patch_curvature_hist(verts, faces, fnormal, face_patch_verts,face_patch_faces, nbins)
%% todo using histogram of 
% conformal factor
% shape diameter function (SDF)
% shape contexts (SC): Shape matching and object recognition using shape contexts
% average geodesic distance (AGD): Topology matching for fully automatic similarity estimation of 3d shapes

%%
DEBUG = 0;
npatch = length(face_patch_verts);
%% compute patch normals, normalized
patch_normal = zeros(npatch,3);
for i=1:npatch
    pfaces = face_patch_faces{i};
    patch_normal(i,:) = mean(fnormal(pfaces,:));
end
d = sqrt( sum(patch_normal.^2,2) ); 
d(d<eps)=1;
patch_normal = patch_normal ./ repmat( d, 1,3 );
%% compute Gauss curvature
if DEBUG
    figure('Name','Current patch');set(gcf,'color','white');view3d rot;
    plot_mesh(verts, faces);hold on;
    options.face_vertex_color=[]; options.face_color=[1,0,0];
end

[Umin,Umax,Cmin,Cmax,Cmean,Cgauss] = compute_curvature(verts, faces);
patch_curvature_hist = zeros(npatch, nbins);

%%
for i=1:npatch
    pverts = face_patch_verts{i};    
    
    if DEBUG    
        ifaces = faces(face_patch==i,:);
        pfaces = face_patch_faces{i};
%         h=scatter3(verts(pverts,1),verts(pverts,2),verts(pverts,3),40,'b','filled');        
        h = plot_mesh(verts(pverts,:), faces(pfaces,:), options);colormap jet;
        delete(h);
    end
    patch_curvature_hist(i,:) = hist(Cgauss(pverts),nbins)/length(pverts);
end

