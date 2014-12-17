% #####################################################
function divergences = compute_div_field_he(vector_field,Mhe, Mifs,sel)
if exist('sel','var') ==0
    sel = 1:Mifs.nov;
end

divergences=zeros(size(sel,2), size(vector_field,2)/3);

for i=1:size(sel,2)
    [divergence, hodge] = div(sel(i), vector_field,Mhe, Mifs);    
    divergences(i,:) = divergence;
end

function [divergence, hodge] = div(vertex_index, vector_field, Mhe, Mifs)
%%
% 
% $$div(f_i)=\sum_{i}f_i.e_i^90$$ or
% $$div(f_i)=\sum_{i}f_i.(cot(j).e_j+cot(k).e_k)$$
% 
%%
[borderV, borderE, borderF] = one_ring_border(vertex_index, Mhe, Mifs);
hodge = edge_rot_90(borderE, Mifs.face_normals(borderF), Mhe, Mifs);
%hodge = e_90(vertex_index, borderF, Mifs);
hodge = 0.5 * hodge;
%%
divergence=zeros(1,size(vector_field,2)/3);
for coord = 1:size(divergence,2)
    vf = vector_field(borderF, 3*(coord-1)+1: 3*coord);
    divs = dot(vf, hodge, 2);
    divergence(coord)=sum(divs);
end

function hodge = e_90(vertex_index, borderF, Mifs)
%%
% return the vector of 
% $$(cot(j).e_j+cot(k).e_k)$$
% 
%%
hodge=zeros(size(borderF,2),3);

i=vertex_index;
for b = 1:size(borderF,2)
    bf = Mifs.face_vertices(borderF(b));
    % compute complementary vertices
    if bf(1)==i
        v = bf(2:3);
    elseif bf(2)==i
        v = bf([3 1]);
    elseif bf(3)==i
        v = bf(1:2);
    else
        error('Problem in face ring.');
    end
    j = v(1); k = v(2);
    vi = Mifs.vertex_coords(i);vj = Mifs.vertex_coords(j);vk = Mifs.vertex_coords(k);

    % angles
    alpha = myangle(vk-vi,vk-vj);
    beta = myangle(vj-vi,vj-vk);
    hodge(b,:) = cot(alpha)*(vi-vj)+cot(beta)*(vi-vk);
end

%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );