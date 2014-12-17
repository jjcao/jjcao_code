% function test_compute_grad_field(fid)
% %% 
% if nargin < 1
%     fid = 1;
% end
clear all;
fid = 3;

test_file = {'isosceles_righttriangle.off','cube.off','hexagon_1.off'};
[vertices,faces]=read_mesh(test_file{fid});
% [Mhe,Mifs] = to_halfedge(vertices, faces,'facenormal', compute_face_normals(vertices,faces));
% assert(Mifs.face_normals_exist,'unit normals of faces do not exist!');

vfs = compute_grad_field(vertices, vertices, faces);
%% display
figure;hold on;
plot_mesh(vertices, faces); 
for i = 1:3
    options.subsample_normal = 1; options.normal_scaling = 1;
    switch i
        case 1
            options.line_speci = '-g';
        case 2
            options.line_speci = '-b';
        case 3
            options.line_speci = '-r';
    end    
    plot_face_normal(vertices, faces, vfs(:,3*(i-1)+1:3*i),options);
end
view3d rot;

%% verify 1
if fid==1
    tmp = [1 0 0 0 1 0 0 0 0];    
    assert(sum(vfs-tmp,2)==0,'error!');
end
%% verify translation-invariant
vertices = vertices+0.1;
vfs1=compute_grad_field(vertices, vertices, faces);
for i=1:size(vfs,1)
    for j=1:size(vfs,2)
        tmp = abs(vfs(i,j)-vfs1(i,j));    
        assert( tmp<1e-7,sprintf('(%u, %u) error: %d', i,j, tmp));
    end
end
%% verify scale-invariant
vertices = (vertices-2)*2;
vfs1=compute_grad_field(vertices, vertices, faces);
for i=1:size(vfs,1)
    for j=1:size(vfs,2)
        tmp = abs(vfs(i,j)-vfs1(i,j));    
        assert( tmp<1e-6,sprintf('(%u, %u) error: %d', i,j, tmp));
    end
end