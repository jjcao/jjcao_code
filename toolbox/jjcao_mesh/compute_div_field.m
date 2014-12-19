function divergences = compute_div_field(vector_field, vs, fs, options)
% compute divergences of each vertex from some piecewise-constant vector
% fields defined on faces of a mesh
%
% $$div(f_i)=\sum_{i}f_i.e_i^90$$ or
% $$div(f_i)=\sum_{i}f_i.(cot(j).e_j+cot(k).e_k)$$
%
% vector_field: vector field defined on each face, 3 columns
% vs: vertices
% fs: faces
% options.rings: vertex faces rings
% options.normals: face normals
% options.areas: 1-ring area of vertex
% Copyright (c) 2008 JJCAO

rings = getoptions(options, 'rings', compute_vertex_face_ring(fs));
normals = getoptions(options, 'normals', compute_face_normals(vs, fs));
areas = getoptions(options, 'areas', area_1_ring(vs, fs, rings));
nvs = size(vs,1);
divergences = zeros(nvs, 1);

for i = 1:nvs
    ring = rings{i};
    nr = length(ring);
    edge = zeros(nr, 2); normal = zeros(nr, 3); vf = zeros(nr, 3);
    for j = 1:nr
        face = fs(ring(j),:);
        if face(1)==i
            edge(j,:) = face(2:3);
        elseif face(2)==i
            edge(j,:) = face([3 1]);
        elseif face(3)==i
            edge(j,:) = face(1:2);
        else
            error('Problem in face ring.');
        end
        
        normal(j,:) = normals(ring(j),:);  
        vf(j,:) = vector_field(ring(j),:);
    end
    hodge = cross(normal,vs(edge(:,2),:) - vs(edge(:,1),:),2);% edge_rot_90 
    
    divergences(i,:) =  sum( dot(vf, hodge, 2) )/areas(i);
end