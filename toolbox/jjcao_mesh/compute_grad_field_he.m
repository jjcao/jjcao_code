function grads = compute_grad_field_he(scalar_field, Mhe, Mifs)
% construct a piecewise-constant gradient field defined on faces from a
% scalar field defined on vertices.
%
% Mhe: the halfedge class;
% Mifs: the indexed face sets;
% 
% Copyright (c) 2008 JJCAO
nscalars =  size(scalar_field,2);
grads=zeros(Mifs.nof(),nscalars*3);
for i = 1:Mifs.nof()    
    area = face_area2(i, Mhe, Mifs);
    
    edges = Mhe.face_edges(i);
    normals = Mifs.face_normals(i); 
    %normals = [normals;normals;normals];
    hodge = edge_rot_90(edges, repmat(normals,3,1), Mhe, Mifs);
    
    vids = Mhe.edge_next(edges);
    vids = Mhe.edge_dest(vids);
    fs = scalar_field(vids,:);

    for j = 1:3
        for coord = 1:nscalars        
            grads(i,3*(coord-1)+1: 3*coord) = grads(i,3*(coord-1)+1: 3*coord) + fs(j,coord)*hodge(j,:);
        end        
    end    
    grads(i,:) = grads(i,:)/area;
    %ind = find(abs(grads(i,:))<1e-10);
    %grads(i,ind) = 0;
end

function area=face_area2(face_id, Mhe, Mifs)
% 2 times area of face
edges = Mhe.face_edges(face_id);
vid1 = Mhe.edge_orig(edges(1));
vid2 = Mhe.edge_dest(edges(1));
vid3 = Mhe.edge_dest(edges(2));

e1=Mifs.vertex_coords(vid2)-Mifs.vertex_coords(vid1);
e2=Mifs.vertex_coords(vid3)-Mifs.vertex_coords(vid2);

area = norm(cross(e1,e2,2));