function [ImgSegmHSV] = disp_draw_segment_hsv(Img, clusterImg, BdSegImg, ...
                                         bDraw, svName)
%
%   Draw Segment Img
%
% -------
% INPUTS:
%
%   Img             : original image.
%   clusterImg      : superpixel image.
%   BdSegImg        : boundary image
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   clImg          : Output segment image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

if nargin<4,
    bDraw = 0 ;
end

% ratio between original image and cluster colors. 
RATIO_ORIG = 0.5 ;

% boundary color
bdcolor = [1 1 1] ;

% dog; cow; aeroplane; grass; boat; face, book; car;
% road; tree; sky
clr_disp = [128 192 128; 0 0 128; 192 128 0; 192 0 0; 0 128 0; 192 64 0; 132 63 0; 64 2 130; ...
            130 63 129; 128 128 0; 128 128 128;] ;
num_clr_dist = size(clr_disp,1) ;

% color to HSV
clr_disp = rgb2hsv(clr_disp) ;
clr_disp(:,3) = 1 * clr_disp(:,3) ;
clr_disp = hsv2rgb(clr_disp)/255 ;

%{
cl_dim = 2 ;
while cl_dim<numCL,
    cl_dim = 2 * cl_dim ;
end
clr_disp = flipud(jet(cl_dim)) ; 
clr_disp = uint8(round(255*clr_disp)) ; 
%}


% number of cluster
numCL = max(clusterImg(:)) ;
% image dimension
[nr nc] = size(clusterImg) ; nch = 3 ;

% assign color
SegmImg = zeros(nr*nc, nch) ;
for k=1:numCL,
    cl_pnt = find(clusterImg==k) ;
    cl_k = mod(k,num_clr_dist) ;
    if cl_k==0, cl_k = num_clr_dist; end
    for l=1:nch,
        SegmImg(cl_pnt,l) = clr_disp(cl_k,l);
    end
end
SegmImg = reshape(SegmImg,[nr, nc, nch]) ;

% change original image
ImgSegmHSV = rgb2hsv(Img);
ImgSegmHSV(:, :, 2) = 0.5*ImgSegmHSV(:, :, 2);
ImgSegmHSV = hsv2rgb(ImgSegmHSV);

% set color 
ImgSegmHSV = ImgSegmHSV*RATIO_ORIG + (1-RATIO_ORIG)*SegmImg;  

if ~isempty(BdSegImg), 
    ImgSegmHSV = reshape(ImgSegmHSV, [nr*nc, nch]) ;
    bd_pnt = find(BdSegImg>0) ;
    for c=1:nch,
        ImgSegmHSV(bd_pnt,c) = bdcolor(c) ;
    end
    ImgSegmHSV = reshape(ImgSegmHSV, [nr nc nch]) ;
end

ImgSegmHSV = uint8(round(255*ImgSegmHSV)) ;

% draw the output
if (bDraw>0), figure; imshow(ImgSegmHSV) ; end

% save the output 
if nargin>4, 
    imwrite(ImgSegmHSV, svName) ;
end

end
