function grads = compute_grad_field(scalar_field, vs, fs, normals)
% construct a piecewise-constant gradient field defined on faces from a
% scalar field defined on vertices.
%
% vs: vertices
% fs: faces
% normals: face normals
%
% Copyright (c) 2009 JJCAO

nscalars =  size(scalar_field,2);
nfaces = size(fs, 1);
grads=zeros(nfaces,nscalars * size(vs,2));

if nargin <4
    normals = compute_face_normals(vs, fs);
end

a = cross(vs(fs(:,2),:) - vs(fs(:,1),:), vs(fs(:,3),:) - vs(fs(:,1),:),2);
farea2s = sum(a.^2,2).^0.5;

for i = 1:nfaces    
    face = fs(i,:);
    edges = [[face(2:3), face(1)]', [face(3), face(1:2)]']; 
    normal = repmat(normals(i,:),3,1); 
    hodge = cross(normal,vs(edges(:,2),:) - vs(edges(:,1),:),2);% edge_rot_90 
    
    for coord = 1:nscalars        
        sf = diag(scalar_field(face(1:3),coord));
        grads(i,3*(coord-1)+1: 3*coord) = sum( sf*hodge);
    end
  
    grads(i,:) = grads(i,:)/farea2s(i);   
end
grads( abs(grads)<1e-10)=0;
