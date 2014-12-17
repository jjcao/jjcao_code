%function eg_segmentation_cf(file,type)
% eg_segmentation_cf - example for mesh segmentation using conformal
% factors
%
%   [] = eg_segmentation_cf(file);
%
%   Example 1: eg_segmentation_cf armadillo_v502.off
%   Example 2: eg_segmentation_cf armadillo_v502.off 'normal cycle'
%
%   Copyright (c) 2009 JJCAO
clear options;

if nargin < 1
    file = 'armadillo_v4326.off';
end
file = 'armadillo_v502.off';%camel,pig1,dancer_v4000_f8000.off

[vertices,faces]=read_mesh(file);
Mhe = to_halfedge(vertices, faces);
boundary = Mhe.boundary_vertices;
if iscell(boundary)
    boundary=cell2mat(boundary)';
end
%% ////////////////////////////////////////////////////////////////////////
% compute original Gaussian curvature and target Gaussian curvature
if nargin<2
    type = 'angle_defect';
end
options.type = type;
options.rings = compute_vertex_face_ring(faces);
options.boundary = boundary;
options.show_error = 0;
%phi = compute_conformal_factor(vertices, faces, options);
phi = csvread('armadillo_v502_cf.txt');
%% find local minimal extrema
[extrema, minv, maxv] = compute_local_extrema(phi, Mhe.vertex_neighbors());
ev = phi(extrema);
threshold = min(ev) + (max(ev)-min(ev))*0.8;
extrema = extrema(logical(ev<threshold));
ev = (ev-minv)/(maxv-minv);
sprintf('num of extrema: %d\n',length(extrema))

%% segmentation by geodesic distance
% perform the front propagation, Q contains an approximate segementation
[D,S,Q] = perform_fast_marching_mesh(vertices, faces, extrema);
atlas = cell(size(extrema));
for i = 1:length(extrema)
    atlas{i} = find(Q==extrema(i));
    D( atlas{i} ) = i;    
end
options.face_vertex_color = D;
%h = plot_mesh(vertices, faces, options);
[Qexact,DQ, voronoi_edges] = compute_voronoi_mesh(vertices,faces, extrema, options);
options.voronoi_edges = voronoi_edges;
plot_fast_marching_mesh(vertices,faces, D, [], options);

%% segmentation by cf
j = 17;
n_elements = histc(phi(atlas{j}),x);
n_elements = n_elements/length(atlas{j});
diff = zeros(size(atlas));
x = -99:1:100;
vrings = compute_vertex_ring(faces);
vj = dilate_selected_vertices(atlas{j},vrings);

for i = 1:length(atlas)
    if i == j, continue, end    
    if isempty( intersect(vj,atlas{i})), continue, end
    
    elem = histc(phi([atlas{j};atlas{i}]),x);
    elem = elem/(length(atlas{j})+length(atlas{i}));
    if size(elem) ~= size(n_elements)
        elem = elem';
    end
    diff(i) = sum( abs(elem - n_elements));
end

for i = 1:length(atlas)
    if diff(i) == 0
        continue
    elseif diff(i)<0.2
        atlas{j} = [atlas{j};atlas{i}];
        atlas{i} = [];
    end
end
%% display 
figure('Name','Segmentation CF');%hold on;
D = zeros( length(vertices),1 ); 
for i = 1:length(extrema)
    D( atlas{i} ) = i;    
end
%D( atlas{j} ) = 1;
%D( [atlas{j};atlas{4}] ) = 1;
options.face_vertex_color = D;
h = plot_mesh(vertices, faces, options);
set(h, 'edgecolor', 'none');
colormap jet(256);

%% display minimal extrema
hold on;
vertices = vertices';extrema = extrema';
h = plot3(vertices(1,extrema),vertices(2,extrema), vertices(3,extrema), 'r.');
set(h, 'MarkerSize', 25);
vertices = vertices';extrema = extrema';
hold off;
