% test_gaussian_smoothing
%
% jjcao @ 2014
%
clc;clear all;close all;
addpath(genpath('../../'));
diagLength = 800;
epsilon = 0.002 * diagLength * 1; 
t = (0.1:20:100)*epsilon^2; 
% t = 55*epsilon^2; 
dist_const = 2.5;
% dist_const = 15;

M.filename = 'wolf2_v1565.off';
[M.verts,M.faces] = read_mesh(M.filename); M.nverts = size(M.verts,1);
[M.verts, diagLength] = normalize_vertex3d(M.verts,diagLength);
nverts = size(M.verts, 1);

%% 
tree = kdtree_build(M.verts);
figure(1);set(gcf,'color','white'); movegui(1, 'northeast');
figure(2);set(gcf,'color','white'); movegui(2, 'southeast');
figure(3);set(gcf,'color','white'); movegui(3, 'east');
trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
     'FaceColor', 'cyan',  'edgecolor','none','faceAlpha', 0.7); axis equal;axis off;mouse3d; hold on;
h = 0;
for i = t
    [verts1, nneigh] = gaussian_smoothing(M.verts, M.verts, dist_const*sqrt(i)*ones(nverts,1), i*ones(nverts,1), tree);
    figure(1);
    trisurf(M.faces,M.verts(:,1),M.verts(:,2),M.verts(:,3), ...
     'FaceColor', 'cyan',  'edgecolor','none','faceAlpha', 0.5);  hold on;
    trisurf(M.faces,verts1(:,1),verts1(:,2),verts1(:,3), ...
     'FaceColor', 'interp',  'edgecolor','r'); 
    axis equal;axis off;mouse3d; hold off; %light('Position',[1 0 0],'Style','infinite');lighting phong;
    
    figure(2);
    plot(nneigh);
 
    figure(3);
    v = M.verts - verts1;
    dist2 = sum(v.^2,2);
    sprintf('mean diff: %f', mean( dist2) ) 
    if i==1 
        if abs(min(dist2)) > 1.0e-10
            warning('when only 1 neighbor is selected for each vertex, mesh should unmoved!')
        else
            sprintf('The code seems correct, since mesh unmoved when only 1 neighbor is selected.')
        end        
    end
    
    if h,delete(h),end    
    sel = find(dist2 > mean(dist2));
    h = quiver3(M.verts(sel,1),M.verts(sel,2),M.verts(sel,3),v(sel,1),v(sel,2),v(sel,3));   
end

kdtree_delete(tree);