function [sp] = img2superpixel(img,spImg)
% img  : ground truth  pixel
% sp : superpixel 0 or 255
spNum = max(max(spImg));
[m,n] = size(img);

[spcnt, spnpix]= compute_sp_info(spImg, spNum);
for i = 1:spNum
        sp(i,1) = img(spcnt(i,1),spcnt(i,2));
end

% [pos] = find(sp~=0);
% sp(pos) = 1;
sp = double(sp);
