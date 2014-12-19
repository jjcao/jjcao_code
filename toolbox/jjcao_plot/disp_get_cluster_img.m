function [clusterImg] = get_cluster_img(spImg, clusterID)
%
%   Get Cluster Img
%
% -------
% INPUTS:
%
%   spImg           : Superpixel image.
%   clusterID       : clusters of SPs
%
% --------
% OUTPUTS:
%
%   clusterImg      : Output segment image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

if nargin<3,
    bDraw = 0 ;
end

[nr nc nch] = size(spImg) ; nch = 3 ;
numSP = max(max(spImg)) ;
numCL = max(clusterID) ;

clusterImg = zeros(nr*nc, 1) ;
for k=1:numSP,
    % 0: no segment is assigned.
    if clusterID(k)==0, continue; end
    hub_reg = find(spImg==k) ;
    clusterImg(hub_reg) = clusterID(k);
end
clusterImg = reshape(clusterImg, [nr, nc]) ;

% fill the boundaries (Optional)
pref_order = 1:numCL ;
[clusterImg] = fill_up_boundary(clusterImg, pref_order) ;

end

% fill the boundaries (Optional)
function [clImg] = fill_up_boundary(clImg, pref_order)

numCL = max(max(clImg)) ;

[nr nc] = size(clImg) ;

for i=1:numCL,
    while 1,
        [rz cz] = find(clImg==0) ;
        if isempty(rz), break; end
        rzm = rz - 1 ; rzp = rz + 1 ;
        rzm(rzm==0) = rzp(rzm==0) ;
        rzp(rzp>nr) = rzm(rzp>nr) ;
        czm = cz - 1 ; czp = cz + 1 ;
        czm(czm==0) = czp(czm==0) ;
        czp(czp>nc) = czm(czp>nc) ;

        rzs = [rzm rzp rz rz] ;
        czs = [cz cz czm czp] ;

        sz = sub2ind([nr nc], rzs, czs) ;
        label_cl = clImg(sz) ;
        
        label_final = zeros(size(sz,1),1) ;
        for j=numCL:-1:1, 
            [rf cf] = find(label_cl==pref_order(j)) ;
            rf = unique(rf) ;
            label_final(rf) = pref_order(j) ;
        end
        clImg(sub2ind([nr nc], rz, cz)) = label_final ;
    end
end

end
