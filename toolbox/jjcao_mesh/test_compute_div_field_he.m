%	Halfedge
%	draw divergences
%	draw differential coordinates

clear;close all;clc;
path('toolbox',path);
%% read mesh
test_file = {'data/8_isosceles_righttriangle.off','data/hexagon.off','data/catHead_v131.off'};
fid = 3;
[vertices,faces]=read_mesh(test_file{fid});
figure;
plot_mesh(vertices, faces);

%% gradient field
[vnormal, fnormal] = compute_normal(vertices,faces);
vnormal = vnormal';
fnormal =  fnormal';
[Mhe,Mifs] = to_halfedge(vertices, faces,'facenormal', fnormal);
assert(Mifs.face_normals_exist,'unit normals of faces do not exist!');
%vector_fields=read_vector_field('golfball_v14_grad.txt');% n*nof cell
vector_fields=compute_grad_field_he(vertices, Mhe, Mifs);

%% divergences
divergences = compute_div_field_he(vector_fields,Mhe, Mifs);
hold on;
normal_scaling = 0.0;%0: no auto scale
quiver3(vertices(:,1),vertices(:,2),vertices(:,3),divergences(:,1),divergences(:,2),divergences(:,3),normal_scaling,'r','LineWidth',2);
%% differential coordinates
options.fill_option='ring';%'face', 'ring'
options.laplacian = 'conformal'; % the other choice is 'combinatorial', 'distance', 'conformal', or, 'authalic'
delta_coords = compute_diff_coords(vertices,faces,options);
hold on;
normal_scaling = 0.0;%0: no auto scale
quiver3(vertices(:,1),vertices(:,2),vertices(:,3),delta_coords(:,1),delta_coords(:,2),delta_coords(:,3),normal_scaling,'g','LineWidth',2);

%% verify
for i=1:size(vertices,1)
    for coord=1:3
        tmp=abs(delta_coords(i,coord)-divergences(i,coord));
        assert( tmp<1e-6,sprintf('(%u, %u) delta coord != divergence error: %d', i,coord, tmp));
    end
end