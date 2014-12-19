function plot_mesh_scalar(verts, faces, scalars, options)
%
% jjcao, 2012.10
options.null = 0;
figname = getoptions(options,'figname','noname');
position = getoptions(options,'position','northwest');
bsaturate = getoptions(options,'bsaturate',1);
tau = getoptions(options,'tau',1.2);

if bsaturate
    options.face_vertex_color = perform_saturation(scalars,tau);
else
    options.face_vertex_color = scalars;
end

figure('name',figname);movegui(position);set(gcf,'color','white');hold on;
plot_mesh(verts,faces, options);
shading interp; camlight; colorbar;
% colormap hsv(256);
colormap jet(256);
% saveas(gcf, [rep name '-cmax.png'], 'png');