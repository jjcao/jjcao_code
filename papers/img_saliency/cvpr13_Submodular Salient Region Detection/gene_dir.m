function [imdir, spdir, saldir,gdir] = gene_dir(database)
% imdir -- test image path
% spdir -- the superpixel label file path
% saldir  -- the output path of the saliency map
% gdir -- the groundtruth path of image

imdir=['../img/input/' database '/'];

if ~isdir(imdir)
    error('No image direction')
end

spdir='../img/output/superpixels/';
if ~isdir(spdir)
    mkdir(spdir);
end
saldir=['../img/output/saliencymap/' database '/']; 
if ~isdir(saldir)
    mkdir(saldir);
end

gdir=['../img/G/' database '/']; 
if ~isdir(gdir)
    mkdir(gdir);
end