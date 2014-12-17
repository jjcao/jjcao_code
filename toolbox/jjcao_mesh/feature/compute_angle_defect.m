function Cgauss = compute_angle_defect(vertices, faces, options)
% compute_angle_defect - compute Gaussian curvature of a mesh
%
%   Cgauss = compute_angle_defect(vertices, faces, options);
%
%   options.type = 'angle_defect', Reffer to Formular 9 of Discrete differential geometry
%   operators for triangulated 2-manifolds_02
%
%   changed by jjcao, 2012
%   Copyright (c) 2009 JJCAO

if nargin < 3
    options = [];
end

if isfield(options, 'rings')
    rings = options.rings;
else
    rings = compute_vertex_face_ring(faces);
end

if isfield(options, 'boundary')
   boundary = options.boundary;
else
   boundary = [];
end
       
if isfield(options, 'conformal_factors')
    Cgauss = compute_angle_defect2(vertices, faces, boundary, rings, options.conformal_factors);
else
    Cgauss = compute_angle_defect1(vertices, faces, boundary, rings);
end

function Cgauss = compute_angle_defect1(vertices, faces, boundary, rings)
n = length(vertices);
is_boundary = zeros(n, 1);
Cgauss = is_boundary;
is_boundary(boundary) = 1;

for i = 1:n
    alpha = 0;
    for b = rings{i}
        % b is a face adjacent to a
        bf = faces(b,:);
        % compute complementary vertices
        if bf(1)==i
            v = bf(2:3);
        elseif bf(2)==i
            v = bf([1 3]);
        elseif bf(3)==i
            v = bf(1:2);
        else
            error('Problem in face ring.');
        end
        j = v(1); k = v(2);
        vi = vertices(i,:);
        vj = vertices(j,:);
        vk = vertices(k,:);
        % angles
        alpha = alpha + myangle(vk-vi,vj-vi);
    end
    
    if is_boundary(i)
        Cgauss(i) = pi - alpha;
    else
        Cgauss(i) = 2*pi - alpha;
    end
end

% compute the cgauss after scaling by conform factor
% 
%
% refer to: Conformal Flattening by Curvature Prescription and Metric Scaling_08
%
%
function Cgauss = compute_angle_defect2(vertices, faces, boundary, rings, conformal_factors)
n = length(vertices);
is_boundary = zeros(n, 1);
Cgauss = is_boundary;
is_boundary(boundary) = 1;

for i = 1:n
    alpha = 0;
    for b = rings{i}
        % b is a face adjacent to a
        bf = faces(b,:);
        % compute complementary vertices
        if bf(1)==i
            v = bf(2:3);
        elseif bf(2)==i
            v = bf([1 3]);
        elseif bf(3)==i
            v = bf(1:2);
        else
            error('Problem in face ring.');
        end
        j = v(1); k = v(2);
        % angles
        alpha = alpha + myangle2(vertices(j,:), vertices(i,:), vertices(k,:),conformal_factors([j;i;k]));
    end
    
    if is_boundary(i)
        Cgauss(i) = pi - alpha;
    else
        Cgauss(i) = 2*pi - alpha;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v)

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );

function beta = myangle2(j,i,k, cf)
% Conformal Flattening by Curvature Prescription and Metric Scaling, say
% the scaling is e^{2cf} and e^{cf}
eij = j-i; eik = k-i; ekj = j-k;
if abs(cf(1)-cf(2))<eps
    sij = exp(cf(1));
else
    sij = ( exp(cf(1)) - exp(cf(2)) )/(cf(1)-cf(2));
end
if abs(cf(3)-cf(2))<eps
    sik = exp(cf(3));
else
    sik = ( exp(cf(3)) - exp(cf(2)) )/(cf(3)-cf(2));
end
if abs(cf(3)-cf(1))<eps
    skj = exp(cf(1));
else
    skj = ( exp(cf(1)) - exp(cf(3)) )/(cf(1)-cf(3));
end
a = sqrt( sum(eij.^2) ) * sij;
b = sqrt( sum(eik.^2) ) * sik;
c = sqrt( sum(ekj.^2) ) * skj;
a = max(a,eps); b = max(b,eps);

beta = acos( 0.5*(a^2+b^2-c^2)/(a*b) );