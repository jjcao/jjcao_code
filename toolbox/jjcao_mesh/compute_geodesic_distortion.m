function GDD = compute_geodesic_distortion(embededVerts, GD)
% GD is geodesic distance matrix of vertices before embedding
% GDD(i) is geodesic distance distortion of embededVerts(i,:)
%
%   Copyright (c) 2012 Junjie Cao
N = length(GD);
ED = L2_distance(embededVerts', embededVerts',1);
% i=1;j=5;
% sqrt(sum((embededVerts(i,:)-embededVerts(j,:)).^2))

GDD = sqrt(sum((ED-GD).^2,2)./(N-1));