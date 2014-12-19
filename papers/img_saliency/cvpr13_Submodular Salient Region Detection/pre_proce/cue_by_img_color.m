function [cue_by_color, smooth_cue_color] = cue_by_img_color(train_data_dir, superpixels,sp_center,sp_inds,rgb_fea,show_para)
% get color cue red & yellow 
%
% changed by jjcao
%
DEBUG=0;

if size(rgb_fea,1)>size(rgb_fea,2)
    rgb_fea = rgb_fea';
end
[height,width] = size(superpixels);
spNum = max(superpixels(:));
segRGB = double(rgb_fea)./max(max(rgb_fea));

%% 计算regionDist，colorDist
locCollection = cell(spNum,1);
for i=1:spNum
    [locCollection{i}(:,1), locCollection{i}(:,2)] = ...
                                        ind2sub([height, width],sp_inds{i});
end

regionDist = zeros(spNum);%R:超像素重心位置之间的距离,位置归一到[0 1]后计算的结果
meanLoc = sp_center/ max(height, width);
for i = 1:spNum
    tmp = repmat(meanLoc(i,:),spNum,1);
    regionDist(i,:) = sqrt(sum( (tmp-meanLoc).^2, 2));
end


%% Color prior from TFP /2.5
%Center prior
centerPriorMap = ones(height,width);
cue_by_color = zeros(spNum,1);

fileID = fopen(train_data_dir,'rb');
T = fread(fileID, 53*53, 'double');
T = reshape(T, 53,53);
ColorPrior = fread(fileID, 20*20, 'double');
ColorPrior = reshape(ColorPrior, 20,20);
fclose(fileID);


for index = 1:spNum
    mx = sp_center(index,1);
    my = sp_center(index,2);
    variance = sum(std(locCollection{index},0,1));
    prior = centerPriorMap(mx,my)*exp(-variance/50);
    
    nR = segRGB(1,index)/(sum(segRGB(:,index))+(1e-6));
    nG = segRGB(2,index)/(sum(segRGB(:,index))+(1e-6));
    x = min(floor(nR/0.05)+1,20);
    y = min(floor(nG/0.05)+1,20);
    prior = prior*(ColorPrior(x,y)+0.5)/1.5;
    cue_by_color(index) = prior;    
end
if DEBUG    
    [cue_by_color_map] = saliency_sp2im(cue_by_color, ...
                                show_para.sp_inds, show_para.sp_num, show_para.m, ...
                                show_para.n, show_para.w);
    figure; subplot(1,2,1); imshow(cue_by_color_map);
    %     plot_prior_saliency(cue_by_color, show_para);
end

if nargout < 2
    return;
end

%% 先验平滑 TFP: Spatial smoothing on the high-level priors
tmp = exp( -regionDist/0.05 );
smooth_cue_color = tmp * cue_by_color;
smooth_cue_color = smooth_cue_color./sum(tmp,2);


if DEBUG        
    [smooth_cue_color_map] = saliency_sp2im(smooth_cue_color, ...
                            show_para.sp_inds, show_para.sp_num, show_para.m, ...
                            show_para.n, show_para.w);
    subplot(1,2,2); imshow(smooth_cue_color_map);
%     plot_prior_saliency(smooth_cue_color, show_para);
end
