% how to use Halfedge, plot_edge, plot_mesh & scatter3
% more info can be found by searcing "Matlab Halfedge" in google

clear;close all;clc;
path('toolbox',path);
%% read mesh
test_file = {'data/8_isosceles_righttriangle.off','data/hexagon.off','data/catHead_v131.off'};
fid = 3;
[verts,faces]=read_mesh(test_file{fid});
plot_mesh(verts, faces); hold on;
set(gcf,'color','white');
axis off; axis equal; set(gcf,'Renderer','OpenGL');
camorbit(0,0,'camera'); axis vis3d; view3d rot;

%% to halfedge, method 1
Mifs = indexedfaceset(verts,faces);
Mhe = halfedge(Mifs);
% boundary edge
bEdges = Mhe.boundary_edges();
boVerts = Mhe.edge_orig(bEdges);
bdVerts = Mhe.edge_dest(bEdges);
edges = [boVerts;bdVerts];
plot_edges(edges, verts, 3);

% boundary vertices
bVerts = Mhe.boundary_vertices(); % the same with % "Mhe.edge_dest(bEdges);"
scatter3(verts(bVerts,1),verts(bVerts,2), verts(bVerts,3),50,'r','filled');

%% to halfedge with face normal, method 2
[normal fnormal] = compute_face_normal(verts, faces);
[Mhe,Mifs] = to_halfedge(verts, faces,'facenormal', fnormal);
assert(Mifs.face_normals_exist,'unit normals of faces do not exist!');





