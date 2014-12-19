%Global step of ARAP iteration

function U=ARAP_Global(EV,Linv,t,C,R)

%%%%%%%%% Initialize the right-hand side of Poisson equation %%%%%%%%%%%
bx=zeros(size(Linv,1),1);
by=zeros(size(Linv,1),1);
for i=1:size(t,1)
    for j=0:2
        j1=j+1;
        j2=mod(j+1,3)+1;
        j3=mod(j+2,3)+1;
        idx2=t(i,j2);
        idx3=t(i,j3);
        Eij(1,1)=EV(i,j1,1);
        Eij(2,1)=EV(i,j1,2);
        
        rot(1,1)=R(i,1);
        rot(1,2)=R(i,2);
        rot(2,1)=R(i,3);
        rot(2,2)=R(i,4);
        bx(idx3) = bx(idx3)+rot(1,:)*Eij*C(i,j1);
        bx(idx2) = bx(idx2)-rot(1,:)*Eij*C(i,j1);
        by(idx3) = by(idx3)+rot(2,:)*Eij*C(i,j1);
        by(idx2) = by(idx2)-rot(2,:)*Eij*C(i,j1);
   end 
end

bx(size(Linv,1)) = 0;   %Using Lagrange multiplication method to solve
by(size(Linv,1)) = 0;   %the linear system under hard positional constraints
                       %which satisfies x(1,:)=[0,0]

%%%%%%%%%%%%%%%%%%%%%%% Solving Poisson equation %%%%%%%%%%%%%%%%%%%%
Ux=Linv*bx;
Uy=Linv*by;
U=[Ux(1:(size(Ux,1)-1),:),Uy(1:(size(Uy)-1),:)];