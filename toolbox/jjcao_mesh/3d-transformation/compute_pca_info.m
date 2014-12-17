function [center, coeff, axis] = compute_pca_info(verts)
center = (sum(verts)/length(verts));
coeff = pca(verts);
vc = center';
for i = 1:3
    axis(2*i-1,:) = vc-coeff(:,i)*0.25;
    axis(2*i,:) = vc+coeff(:,i)*0.25;
end

