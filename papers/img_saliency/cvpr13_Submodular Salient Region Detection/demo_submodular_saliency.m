% Demo of cvpr13_Submodular Salient Region Detection
%
% There are 3 priors in the paper: face, color and center. We only have
% center and color.
%
%
% by Guangyu Zhong & JJCAO @ 2014

clc;clear all;close all;
addpath(genpath('.'));
addpath(genpath('../../../toolbox/jjcao_common'));
addpath(genpath('../../../toolbox/jjcao_img'));
addpath(genpath('../../../toolbox/jjcao_plot'));
addpath(genpath('../../../toolbox/kdtree'));
DEBUG = 1;

train_data_dir = 'TrainedData';
%% parameters
graphOpts.featMode=1; % 3的效果往往不好！！
if graphOpts.featMode == 3
    graphOpts.feat_dist_opt =  'SDMVC';
end
graphOpts.k = 8; %% k nearest neighbor

priorOpts.enlarge_convex_ratio = 0.1;  % 扩大前景凸包用，默认用最显著的50个角点，实际使用的最显著角点数量 = 50 + 50*enlarge_convex_ratio, 可能需要调整
priorOpts.bseed_percentage =1;        % 背景点占背景的个数的百分比
priorOpts.bseed_thresh = 26;     % 背景点阈值
priorOpts.bseed_sort = 1;        % 背景点怎么取.  1: 凸包外面的点的harris值排序取, 0: 随机取

seed_num = 200;    % seed number
sp_num_max = 900;    % inital superpixel number (upper bound number)
theta = 10; % control the edge weight 

%% gene dir and file name
database ='MSRA1000';%Berkeley300, MSRA1000
[imdir, spdir, saldir,gdir] = gene_dir(database);
imnames=dir([imdir '*' 'bmp']);
gnames = dir([gdir '*' 'bmp']);
%%
% 14, 57
for ii=1:length(imnames)
    close all;    clear w;
    %% read image, generate superpixels, compute features
    imname=imnames(ii).name;
    
    sprintf('%d: %s', ii, imname)
    [input_im] = imread([imdir, imname]);  
    [m,n,~] = size(input_im); 
    w = [m,n,1,m,1,n];
    [m,n,k] = size(input_im);
    
    %    
    [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = gene_superpixel(imdir, imname, sp_num_max, spdir, m, n);
%     [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = run_superpixel(input_im, sp_num_max) ;  % turbopixel虽然质量更好，但是含有border，必须自己写代码去掉才能用         
    sp_num = size(spAdjcMat,1);    
    
    % 
    [sp_fea, rgb_fea]  = gene_feature(input_im, superpixels, sp_center, sp_npix, graphOpts);        
%     min_sp_fea = min(sp_fea(:)); min_rgb_fea = min(rgb_fea(:));
%     sp_fea=(sp_fea-min_sp_fea)/(max(sp_fea(:))-min_sp_fea);
%     rgb_fea=(rgb_fea-min_rgb_fea)/(max(rgb_fea(:))-min_rgb_fea);
        
    %% get pos_label_inds
    % build affinity matrix
    [adj_mat,sigma] = build_affinity_matrix_knearest(sp_fea, graphOpts.k);
    aff_mat = gene_weight( adj_mat, sp_fea, sigma);
    NUM_EVAL_PNT = 100;
    [pos_label_inds aggloClID] = gen_eval_list(sp_fea, adj_mat, ...
    superpixels, spAdjcMat, NUM_EVAL_PNT) ;
     disp_draw_intra_agglo_real(input_im, superpixels, pos_label_inds, DEBUG);
    %% compute prior
    cueImgCnt = cue_by_img_center(sp_center, superpixels);%sp_center, sp_npix
%     cueImgCnt = ones(sp_num,1);    % Berkely300,不加中心先验会好一点点
    show_para = [];
    show_para.sp_inds = sp_inds;
    show_para.sp_num = sp_num;
    show_para.m = m;
    show_para.n = n;
    show_para.w = w; 
    [cue_by_color,smooth_cue_color] = cue_by_img_color(train_data_dir, superpixels,sp_center,sp_inds,rgb_fea,show_para); 
    if DEBUG            
        cue_by_color_im = saliency_sp2im(cue_by_color, sp_inds, sp_num, m, n, w);
        figure('name', 'cue_by_color_im'); subplot(2,2,1); imshow(cue_by_color_im);       
        cue_by_color_im = saliency_sp2im(smooth_cue_color, sp_inds, sp_num, m, n, w);
        subplot(2,2,2); imshow(cue_by_color_im);       
        cue_by_color_im = saliency_sp2im(cueImgCnt.*cue_by_color, sp_inds, sp_num, m, n, w);
        subplot(2,2,3); imshow(cue_by_color_im);       
        cue_by_color_im = saliency_sp2im(cueImgCnt.*smooth_cue_color, sp_inds, sp_num, m, n, w);
        subplot(2,2,4); imshow(cue_by_color_im);   
    end
    prior_mode = 'high_prior';
    switch prior_mode
        case 'high_prior'
            %      prior = [cue_by_color;0];  % 已经在函数内部加完中心先验了
            prior = [cueImgCnt.*cue_by_color; 0];
        case 'ground'
            gname=gnames(ii).name;
            ground_map = imread([gdir, gname]);
            prior = img2superpixel(double(ground_map),superpixels);
            prior = [prior./255;0];
    end     
    
   %% select salient seeds and regions by submodular    
    K = 5;    lambda =0; C_theta = 0.05;
%     pos_label_inds = [superpixels(197,157),superpixels(77,142),superpixels(88,126),superpixels(158,115),superpixels(186,103)];
% pos_label_inds = [451,60,540,101];  
    [fa_location, fa_assig,K] = submodular_salient_region_detection(aff_mat,prior,K, lambda, pos_label_inds,sp_num,C_theta);
        
    %
    outname=[saldir imnames(ii).name(1:end-4) '_ssrd.png'];
    disp_draw_intra_agglo_real(input_im, superpixels, fa_location, DEBUG,outname);
    if DEBUG            
        region_im = saliency_sp2im(fa_assig, sp_inds, sp_num, m, n, w);
        figure('name', 'regions'); imshow(region_im);       
    end
    
    %% 3.4. Saliency Map Construction
    [f, sregion] = compute_saliency_map(fa_location, fa_assig, sp_fea, sp_center);

    if DEBUG
        F = zeros(sp_num,1);
        for i=1:length(f)
            F(sregion{i}) = f(i);
        end
        tmp = saliency_sp2im( F, sp_inds, sp_num, m, n, w);
        figure('name', 'region saliency');
        imshow(tmp);
    end
    %
    sp_saliency = compute_S(fa_location, f, sp_fea, 0.2); % sigma*20
    im_saliency = saliency_sp2im( sp_saliency, sp_inds, sp_num, m, n, w);
     
    %% output
    if DEBUG
        figure('name', 'saliency map');
        imshow(im_saliency);
    end
    outname=[saldir imnames(ii).name(1:end-4) '.bmp'];
    imwrite(im_saliency,outname);
    
end 
% test_evaluate