function [ImgBdSeg] = disp_draw_imgs_clust_boundary(Img, BdSegImg, bDraw, svName)
%
%   Draw Img with segment boundary
%
% -------
% INPUTS:
%
%   spImg           : Superpixel image.
%   BdSegImg        : boundary image
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   ImgBdSeg        : Original image + boundary
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

if nargin<3,
    bDraw = 0 ;
end

% boundary color
bdcolor = [255 0 0] ;

[nr nc nch] = size(Img) ; nch = 3 ;

% boundary image
ImgBdSeg = Img ;
ImgBdSeg = reshape(ImgBdSeg, [nr*nc, nch]) ;
bd_pnt = find(BdSegImg>0) ;
for c=1:nch,
    ImgBdSeg(bd_pnt,c) = bdcolor(c) ;
end
ImgBdSeg = reshape(ImgBdSeg, [nr nc nch]) ;

% draw the output
if (bDraw>0), 
    figure; imshow(ImgBdSeg) ;
end

% save the output 
if nargin>3, 
    imwrite(ImgBdSeg, svName) ;
end

end

