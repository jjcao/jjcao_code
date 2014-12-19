function [CutImg] = disp_draw_imgs_clust_cut(Img, clusterImg, drawClID, bDraw, svName)
%
%   Draw selected segments of an image with gray background
%
% -------
% INPUTS:
%
%   spImg           : Superpixel image.
%   clusterImg      : clusters of SPs
%   drawClID        : cluster ID to be drawn. given by cell
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   CutImg          : Output segment image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

% As default, draw all clusters
if isempty(drawClID),
    numCluster = max(clusterImg(:)) ;
    drawClID = 1:numCluster ;
    drawClID = num2cell(drawClID) ;
end

if nargin<4,
    bDraw = 0 ;
end

bg_color = [255 255 255] ;

nImgDraw = length(drawClID) ;

[nr nc nch] = size(clusterImg) ;
Img = imresize(Img, [nr nc]) ;
Img = reshape(Img, [nr*nc 3]) ;

CutImg = cell(1,nImgDraw) ;
for i=1:nImgDraw,
    CutImg{i} = uint8(repmat(bg_color, [nr*nc 1])) ;

    for j=1:length(drawClID{i}),    
            dInd = find(clusterImg==drawClID{i}(j)) ;
        for l=1:3,
            CutImg{i}(dInd,l) = Img(dInd,l) ;
        end
    end
    CutImg{i} = reshape(CutImg{i}, [nr nc 3]) ;
end



% draw the output
if (bDraw>0), 
    figure; 
    nFigC = min(4,nImgDraw) ;
    nFigR = ceil(nImgDraw/nFigC) ;
    for i=1:nImgDraw;
        subplot(nFigR, nFigC, i); imshow(CutImg{i}) ;
    end
end


% save the output 
if nargin>4, 
    for i=1:nImgDraw;
        imwrite(CutImg{i}, strrep(svName,'.png', ['_' num2str(i) '.png'] )) ;
    end
end

end

