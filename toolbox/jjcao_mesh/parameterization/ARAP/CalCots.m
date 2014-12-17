%Calculate cot weights of each angle in each triangle of the mesh
%                                      3
%                                      /\
%                                     /  \
%                                    /    \
%                                 c /  t   \ b
%                                  /        \
%                                 /__________\
%                                1      a     2 
%
%

function C=CalCots(x,t)

C=[];
for i=1:size(t,1)
    idx1=t(i,1);
    idx2=t(i,2);
    idx3=t(i,3);
    v1=x(idx1,:);
    v2=x(idx2,:);
    v3=x(idx3,:);
    a=norm(v1-v2);
    b=norm(v2-v3);
    c=norm(v3-v1);
    angle1=acos((a^2+c^2-b^2)/(2*a*c));
    angle2=acos((a^2+b^2-c^2)/(2*a*b));
    angle3=acos((b^2+c^2-a^2)/(2*b*c));
    C(i,1)=cot(angle1);
    C(i,2)=cot(angle2);
    C(i,3)=cot(angle3);
end