% Compute the rigid energy
% formula (1) of sgp08_A Local Global Approach to Mesh Parameterization

function E=CalRigidEnergy(EV,t,U,C,R)

E=0;
for i=1:size(t,1)
    for j=0:2
        j1=j+1;
        j2=mod(j+1,3)+1;
        j3=mod(j+2,3)+1;
        idx2=t(i,j2);
        idx3=t(i,j3);
        Uij=U(idx3,:)'-U(idx2,:)';
        Eij(1,1)=EV(i,j1,1);
        Eij(2,1)=EV(i,j1,2);
        
        rot(1,1)=R(i,1);
        rot(1,2)=R(i,2);
        rot(2,1)=R(i,3);
        rot(2,2)=R(i,4);
        E=E+C(i,j1)*norm(Uij-rot*Eij)^2;
   end 
end