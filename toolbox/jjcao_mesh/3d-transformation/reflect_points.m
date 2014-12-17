function vertsNew = reflect_points(verts, center, normal)

n2 = size(verts,1);
normalMat = repmat(normal,n2,1);
centerMat = repmat(center,n2,1);
len = 2*sum(normalMat.*(centerMat-verts),2);
vertsNew = verts + repmat(len,1,3).*normalMat;