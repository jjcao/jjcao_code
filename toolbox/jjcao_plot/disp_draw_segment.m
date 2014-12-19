function [ImgSegm] = disp_draw_segment(clusterImg, bDraw, svName)
%
%   Draw Segment Img
%
% -------
% INPUTS:
%
%   clusterImg      : Superpixel image.
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   ImgSegm        : Output segment image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

if nargin<2,
    bDraw = 0 ;
end

% dog; cow; aeroplane; grass; boat; face, book; car;
% road; tree; sky
clr_disp = [128 192 128; 0 0 128; 192 128 0; 192 0 0; 0 128 0; 192 64 0; 132 63 0; 64 2 130; ...
            130 63 129; 128 128 0; 128 128 128;] ;
num_clr_dist = size(clr_disp,1) ;

% number of cluster
numCL = max(clusterImg(:)) ;
% image dimension
[nr nc] = size(clusterImg) ; nch = 3 ;

% assign color
ImgSegm = zeros(nr*nc, nch, 'uint8') ;
for k=1:numCL,
    cl_pnt = find(clusterImg==k) ;
    cl_k = mod(k,num_clr_dist) ;
    if cl_k==0, cl_k = num_clr_dist; end
    for l=1:nch,
        ImgSegm(cl_pnt,l) = clr_disp(cl_k,l);
    end
end

ImgSegm = reshape(ImgSegm,[nr, nc, nch]) ;

% draw the output
if (bDraw>0), figure; imshow(ImgSegm) ; end

% save the output 
if nargin>2, 
    imwrite(ImgSegm, svName) ;
end

end
