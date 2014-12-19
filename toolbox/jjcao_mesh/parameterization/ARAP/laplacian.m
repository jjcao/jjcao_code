%Compute cot Laplacian of the mesh

function L=laplacian(x, t, C)

n=size(x, 1);
L=sparse(n, n);
for i=1:size(t,1)
    for j=0:2
        j1=j+1;
        j2=mod(j+1,3)+1;
        j3=mod(j+2,3)+1;
        idx1=t(i,j1);
        idx2=t(i,j2);
        idx3=t(i,j3);
        L(idx1,idx1)=L(idx1,idx1)+C(i,j2)+C(i,j3);
        L(idx1,idx2)=L(idx1,idx2)-C(i,j3);
        L(idx1,idx3)=L(idx1,idx3)-C(i,j2);
    end
end

   