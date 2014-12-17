%Local step of ARAP iteration

function R=ARAP_Local(u,t,EV,C)

R=[];
for i=1:size(t,1)
    idx1=t(i,1);
    idx2=t(i,2);
    idx3=t(i,3);
    u1=u(idx1,:); 
    u2=u(idx2,:); 
    u3=u(idx3,:);
    uu(1,:)=u3-u2;
    uu(2,:)=u1-u3;
    uu(3,:)=u2-u1;
    xx(1,1)=EV(i,1,1);
    xx(1,2)=EV(i,1,2);
    xx(2,1)=EV(i,2,1);
    xx(2,2)=EV(i,2,2);
    xx(3,1)=EV(i,3,1);
    xx(3,2)=EV(i,3,2);
    cc(1,1)=C(i,1);
    cc(2,2)=C(i,2);
    cc(3,3)=C(i,3);
    CovMat=xx'*cc*uu;    %Compute the coviarance matrix
    
    [s,v,d]=svd(CovMat);
    rot=d*s';
    if det(rot)<0                %Eliminate the flip matrix, which is 
        if v(1,1)<v(2,2)         %[a,b;b,-1] instead of [a,b;-b,a].
            s(:,1)=-1.0*s(:,1);
        else
            s(:,2)=-1.0*s(:,2);
        end
        rot=d*s';
    end
    R(i,1)=rot(1,1);
    R(i,2)=rot(1,2);
    R(i,3)=rot(2,1);
    R(i,4)=rot(2,2);
end
             