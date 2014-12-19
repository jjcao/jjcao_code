% function  GaussWeighcurvature= compute_gaussian_weighted_curvature(vertex,Cmean,delta,tree)
% c++ comments
%

%% matlab implementation
% function  GaussWeighcurvature= compute_gaussian_weighted_curvature(vertex,Cmean,delta,tree)
% 
% row = size(vertex,1); 
% GaussWeighcurvature=zeros(row,2);
% 
% qradi = 2*delta;
% qradii = 2*qradi;
% 
% for i=1:row
%     qpoint = vertex(i,:);  
%     [idxs, dists] = kdtree_ball_query( tree, qpoint, qradi);
%     tmp = exp(-dists.^2/(2*delta^2));
%     GaussWeighcurvature(i,1) = sum(Cmean(idxs).*tmp);
%     GaussWeighcurvature(i,1) = GaussWeighcurvature(i,1)/sum(tmp);
% 
%     [idxss, distss] = kdtree_ball_query( tree, qpoint, qradii);
%     tmp = exp(-distss.^2/(2*qradi^2));
%     GaussWeighcurvature(i,2) = sum(Cmean(idxss,:).*tmp);
%     GaussWeighcurvature(i,2) = GaussWeighcurvature(i,2)/sum(tmp);
% end
  
