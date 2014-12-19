function [superpixels, spAdjcMat,sp_inds, sp_center, sp_npix] = gene_superpixel(imdir, imname, sp_num_max, spdir, m, n)
%
% sp_center: center of each superpixel
% sp_npix: number of pixel of each superpixel
%
% jjcao @ 2013
%

if ~strcmp(imname(end-3:end),'.bmp')
    imname=[imname(1:end-4) '.bmp'];% the slic software support only the '.bmp' image
end

comm=['SLICSuperpixelSegmentation' ' ' [imdir imname] ' ' ...
    int2str(20) ' ' int2str(sp_num_max) ' ' spdir];
system(comm);
spname=[spdir imname(1:end-4)  '.dat'];
superpixels=read_dat([m,n],spname); % superpixel label matrix
sp_num=max(superpixels(:)); % the actual superpixel number

spAdjcMat = build_sp_adjacent_matrix(superpixels,sp_num);

if nargout > 3
    sp_inds=cell(sp_num,1);
    sp_center = zeros(sp_num,2);
    sp_npix  = zeros(sp_num,1);
    for j=1:sp_num,
        sp_inds{j}=find(superpixels==j);
        [r c] = ind2sub([m,n], sp_inds{j});
        sp_npix(j) = length(r);
        sp_center(j,:) = round([mean(r) mean(c)]);
    end
end