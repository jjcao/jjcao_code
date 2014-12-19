function [cntSeg, npixSeg, contourSeg]= compute_sp_info(spImg, nSeg)
%
% center coordinates & number of pixel of each segment
%
% jjcao @ 2012

cntSeg = zeros(nSeg,2);
npixSeg = zeros(nSeg,1);

if nargout > 2 
    imsize = size(spImg);
    contourSeg = cell(nSeg,1);
end



for j=1:nSeg,
    [r c] = find(spImg==j) ;
    npixSeg(j) = length(r);
    cntSeg(j,:) = round([mean(r) mean(c)]);
    
    if nargout > 2 
        ind = sub2ind(imsize, r, c);
        BW = zeros(imsize);
        BW(ind) = 1;
        B = bwboundaries(BW,'noholes');
        contourSeg{j} = B{1};

    %     imshow(label2rgb(spImg, @jet, [.5 .5 .5]));  hold on;
    %     for k = 1:length(B)
    %         boundary = B{k};
    %         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    %     end
    end
end
    
end