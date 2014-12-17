function [awks,landmark]=compute_awks(verts,faces,nbr_landmarks,DEBUG)

%% Laplacian eigen
nbasis = 300;
withAreaNormalization = true;
adjustL = false;
[M.eigvector,M.eigvalue, A]=compute_Laplace_eigen(verts,faces,nbasis, withAreaNormalization, adjustL);

%% WKS
wks = WKS(M.eigvector, M.eigvalue);

%% compute awks(v)=sum£ûd(v,b_i)*area(b_i)£ý
[face_patch, patch_area,landmark]=supervertex_by_farthest_sampling(verts,faces,nbr_landmarks);

k=max(size(landmark,1),size(landmark,2));
D=zeros(size(verts,1),k);
for i=1:k
    D1=abs(wks-repmat(wks(landmark(i),:),size(wks,1),1));
    D2=abs(wks+repmat(wks(landmark(i),:),size(wks,1),1));
    D(:,i)=sum(D1./D2,2);
end

awks=D*patch_area;  %awks(v)=¶ÔiÇóºÍ£ûd(v,b_i)*area(b_i)£ý
