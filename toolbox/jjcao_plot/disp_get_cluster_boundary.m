function [BdSegImg] = get_cluster_boundary(clusterImg)
%
%   Draw segment boundary image. 
%
% -------
% INPUTS:
%
%   clusterImg      : Cluster image.
%
% --------
% OUTPUTS:
%
%   BdSegImg        : cluster boundary image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%


% get the binary boundary image (1: boundaries)
[cx,cy] = gradient(clusterImg);
bdImg = (abs(cx)+abs(cy))~=0;
bdImg = bwmorph(bdImg,'skel',Inf) ;
clear cx cy

% thicken the line
BdSegImg = bwmorph(bdImg, 'dilate',1);

end

