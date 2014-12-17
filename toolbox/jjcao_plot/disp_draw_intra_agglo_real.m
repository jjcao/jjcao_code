function [label_img] = disp_draw_intra_agglo_real(Img, spImg, EvalSet, bDraw, svName)
%
%   draw agglomerative clustering results on real image.    
%
% -------
% INPUTS:
%
%   spImg           : Superpixel image.
%   cntSeg          : center of each SP
%   clID            : Cluster output
%   EvalSet         : SPs to be evaluated.
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   label_img       : Output segment image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%


if nargin<4,
    bDraw = 0 ;
end

% boundary color
bdcolor = [0 0 0] ;

% get the binary boundary image (1: boundaries)
[cx,cy] = gradient(spImg);
bdImg = (abs(cx)+abs(cy))~=0;
bdImg = bwmorph(bdImg,'skel',Inf) ;
clear cx cy

% thicken the line
%bdImg = bwmorph(bdImg, 'dilate');

% get colors for evaluation points.
mcolor = jet(256) ;
mcolor = round(255*mcolor(randperm(256),:)) ;

% get original image
label_img = Img ;

[nr nc nch] = size(Img) ;
numSP = max(max(spImg)) ;
% color the evaluation points
label_img = reshape(label_img, [nr*nc, nch]) ;
for i=1:length(EvalSet),
    hub_reg = find(spImg==EvalSet(i)) ;
    for l=1:nch,
        label_img(hub_reg,l) = mcolor(i,l) ;
    end
end

% draw boundary.
bd_pnt = find(bdImg>0) ;
for c=1:nch,
    label_img(bd_pnt,c) = bdcolor(c) ;
end



label_img = reshape(label_img, [nr nc, nch]) ;
% label_img = reshape(label_img, [nr nc]) ;

if bDraw>0,
    figure; imshow(label_img) ; 
end

% save the temperature images. 
if nargin>4, 
    imwrite(label_img, svName) ;
end

