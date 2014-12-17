function cueImgCnt = cue_by_img_center(spCnt, spImg)
%
% reweight superpixel saliency according to its position relative to image center
%
% jjcao @ 2013
%
nSp = length(spCnt);
iw = size(spImg, 2)*0.5;
ih = size(spImg, 1)*0.5;
center = [ih,iw];

dist2 = spCnt - repmat(center, nSp,1);
dist2 = sum(dist2.^2,2);
sigma2 = 2*mean(dist2);

cueImgCnt = exp(-dist2/sigma2);

%% % old
% function cueImgCnt = cue_by_img_center(spCnt, spImg)
% %
% % reweight superpixel saliency according to its position relative to image center
% %
% % jjcao @ 2013
% %
% nSp = length(spCnt);
% 
% %% according to image center
% iw = size(spImg, 2)*0.5;
% tmp1 = spCnt(:,2) - repmat(iw, nSp,1);% x
% ww = exp( -tmp1.^2/(iw^2) );
% 
% ih = size(spImg, 1)*0.5;
% tmp2 = spCnt(:,1) - repmat(ih, nSp,1);% x
% wh = exp( -tmp2.^2/(ih^2) );
% 
% cueImgCnt = ww.*wh;    

