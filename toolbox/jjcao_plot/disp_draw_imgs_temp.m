function [imGray] = disp_draw_imgs_temp(flimg, spImg, Tset, bDraw, svName)
%
%	Draw images (whose paths: 'img_path') on a single big canvas ('Canvas').
%   Each of 'img_path' is resized to 'img_size' and placed on a
%   rectangular grid (length(img_path)/MAX_NUM_H X MAX_NUM_H).
%   
%
% -------
% INPUTS:
%
%   flimg           : The list of images.
%   spImg           : The list of superpixel images.
%   Tset            : The list of temperature
%   bDraw           : 1 = draw the output, 0: do not draw
%   svName          : if it is empty, do not save the output.
%                     if it is a filename, save the output.
%
% --------
% OUTPUTS:
%
%   imGray          : Output temperature image.
%
% -------------
%
%   Gunhee Kim (gunhee@cs.cmu.edu)
%

if nargin<4,
    bDraw = 0 ;
end

[n_sp n_hub] = size(Tset) ;

% gray value of each SP
Tgray = zeros(n_sp, n_hub, 'uint8') ;
for i=1:n_hub
    Tset(:,i)  = Tset(:,i) - min(Tset(:,i)) ; 
    Tgray(:,i) = uint8(round(255*Tset(:,i)/max(Tset(:,i))) + 1) ;
end

% load image
colormaps = gray(256) ;
img = imread(flimg) ;
[nr nc] = size(spImg) ;
nch = size(img, 3) ;
img = imresize(img, [nr nc]) ;

% draw the temperature images by subplot
if (bDraw>0) figure; end
imGray = cell(1,n_hub) ;
for i=1:n_hub,
    imGray{i} = zeros(nr*nc, 1, 'uint8') ;
end
for i=1:n_sp,
    reg_ind = find(spImg==i) ;
    for j=1:n_hub,
        imGray{j}(reg_ind,:) = Tgray(i,j) ;
    end
end


% fill the boundaries
for i=1:n_hub,
    imGray{i} = reshape(imGray{i}, [nr nc]) ;
end
[imGray] = fill_up_boundary(imGray) ;

minC = 20; maxC = 250;
slopeC = (maxC-minC)/255 ;
for i=1:n_hub,
    imGray{i} = colormaps( slopeC*round(imGray{i})+minC,:) ;
    imGray{i} = reshape(imGray{i}, [nr, nc, 3]) ;
    imGray{i} = uint8(round(255*imGray{i})) ;
end

if (bDraw>0) 
    for i=1:n_hub,
        subplot(ceil(n_hub/3), 3, i); imshow(imGray{i}) ; 
    end
end

% save the temperature images. 
if nargin>4, 
    for i=1:n_hub,
        imwrite(imGray{i}, [svName(1:end-4) '_' num2str(i) svName(end-3:end)] ) ;
    end
end

end

% fill the boundaries (Optional)
function [imGray2] = fill_up_boundary(imGray)

n_hub = length(imGray) ;

imGray2 = cell(1,n_hub) ;
for i=1:n_hub,
    [nr nc] = size(imGray{i}) ;
    imGray2{i} = reshape(imGray{i}, [nr nc]) ;
    while 1,
        [rz cz] = find(imGray2{i}==0) ;
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
        imGray2{i}(sub2ind([nr nc], rz, cz)) =max(imGray2{i}(sz),[],2) ;
    end
end

end



