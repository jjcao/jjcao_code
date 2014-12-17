function spAdjcMat = build_sp_adjacent_matrix(M,N)
% $Description:
%    -compute the adjacent matrix
% $Agruments
% Input;
%    -M: superpixel label matrix
%    -N: superpixel number 
% Output:
%    -spAdjcMat: adjacent matrix
% copyright @ LU13
% guangyuzhong add 12/10/2013

spAdjcMat = zeros(N,N);
[m n] = size(M);

for i = 1:m-1
    for j = 1:n-1
        if(M(i,j)~=M(i,j+1))
            spAdjcMat(M(i,j),M(i,j+1)) = 1;
            spAdjcMat(M(i,j+1),M(i,j)) = 1;
        end;
        if(M(i,j)~=M(i+1,j))
            spAdjcMat(M(i,j),M(i+1,j)) = 1;
            spAdjcMat(M(i+1,j),M(i,j)) = 1;
        end;
        if(M(i,j)~=M(i+1,j+1))
            spAdjcMat(M(i,j),M(i+1,j+1)) = 1;
            spAdjcMat(M(i+1,j+1),M(i,j)) = 1;
        end;
        if(M(i+1,j)~=M(i,j+1))
            spAdjcMat(M(i+1,j),M(i,j+1)) = 1;
            spAdjcMat(M(i,j+1),M(i+1,j)) = 1;
        end;
    end;
end;    

    