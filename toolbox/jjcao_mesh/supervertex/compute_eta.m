function eta = compute_eta(verts, faces, E1, E2, ind, DEBUG)
%
% % bad
% concave_weight = 0.1;
% convex_weight = 0.2;
% normal_weight = 1;

%% good
concave_weight = 1;
convex_weight = 0.8;
normal_weight = 0.1;

%%
eta = ones(length(ind), 1)*normal_weight;
[concave_id,convex_id] = compute_concave_convex_vertices(verts,faces);
if DEBUG
    figure('Name','Concave vertices'); set(gcf,'color','white');hold on;
    h = plot_mesh(verts,faces);
    h=scatter3(verts(concave_id,1),verts(concave_id,2),verts(concave_id,3),80,'b','filled');
    h=scatter3(verts(convex_id,1),verts(convex_id,2),verts(convex_id,3),80,'r','filled');
    view3d zoom;
end

edges = [E1(ind),E2(ind)];
cidx = edges_sharing_vertex_property(edges, faces, concave_id, 1);
conidx = edges_sharing_vertex_property(edges, faces, convex_id, 1);

eta(cidx) = concave_weight;
eta(conidx) = convex_weight;