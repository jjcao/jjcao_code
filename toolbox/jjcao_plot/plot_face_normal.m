function h = plot_face_normal(vertices, faces, normals, options)
%
% JJCAO 2009
n = size(faces,1);
subsample_normal = getoptions(options, 'subsample_normal', min(4000/n,1) );
normal_scaling = getoptions(options, 'normal_scaling', 0.8 );
lineSpeci = getoptions(options, 'line_speci', '-r' ); 
sel = randperm(n); sel = sel(1:floor(end*subsample_normal));
bc = (vertices(faces(:,1),:) + vertices(faces(:,2),:) + vertices(faces(:,3),:))/3;%barycenter of each face

normals = normals.*normal_scaling;

hold on;
h = quiver3(bc(sel,1),bc(sel,2),bc(sel,3),normals(sel,1),normals(sel,2),...
    normals(sel,3),0,lineSpeci);
hold off;