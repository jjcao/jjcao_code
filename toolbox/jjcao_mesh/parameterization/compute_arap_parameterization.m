function [vt iterations] = compute_arap_parameterization(verts, faces, vt, tolerance)
%   changed from ARAP.m
%
%   Copyright (c) 2012 Junjie Cao
if nargin < 4
    tolerance = 0.001;
end

%%%%%%%%%%%%%%%%%%%% Pre-Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EV=CalEdgeVectors(verts,faces);   %Calculate edge vectors of all the triangles
C=CalCots(verts,faces);   %Calculate cot weights of each angle in each triangle 
L=laplacian(verts,faces,C);  %Compute cot Laplacian of the mesh
L(1,size(L,1)+1)=1;
L(size(L,1)+1,1)=1;
Linv = inv(L);

%%%%%%%%%%%%%%%%%%%% Local/Global Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=-1;
Epre=0;
iterations=0;
while abs(Epre-E)>tolerance     %Iteration stops when the energy converges 
    iterations=iterations+1;    %to the tolerance, which is specified by
    Epre=E;                     %user.
    R=ARAP_Local(vt,faces,EV,C);
    vt=ARAP_Global(EV,Linv,faces,C,R);
    E=CalRigidEnergy(EV,faces,vt,C,R);
%%%%%%%%%% Plot iteration results %%%%%%%%%%%%%%%%%%%%%%%%
    %figure
    %trimesh(faces,vt(:,1),vt(:,2));
    %axis equal;
    %title1=num2str(iterations);
    %title2=[title1,' Iterations'];
    %title(title2);    
end
