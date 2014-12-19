function  [step,ori_kseedl,ori_kseeda,ori_kseedb,ori_kseedx,ori_kseedy] = origkseed(imgL,imga,imgb,k)
%% 
% original kseeds  color&position
% Guangyu Zhong 2013.4.15
%%
[width, height] = size(imgL);
N = numel(imgL);
step = fix(sqrt(N/k));
nx = fix(width/step);
ny = fix(height/step);
ori_kseedx = zeros(nx,ny);
ori_kseedy = zeros(nx,ny);
ori_kseedl = zeros(nx,ny);
ori_kseeda = zeros(nx,ny);
ori_kseedb = zeros(nx,ny);
for i = 1:nx-1
    for j = 1:ny-1
        ori_kseedx(i,j) = (i-1)*step+fix(step/2);
        ori_kseedy(i,j) = (j-1)*step+fix(step/2);
        ori_kseedl(i,j) = imgL(ori_kseedx(i,j),ori_kseedy(i,j));
        ori_kseeda(i,j) = imga(ori_kseedx(i,j),ori_kseedy(i,j));
        ori_kseedb(i,j) = imgb(ori_kseedx(i,j),ori_kseedy(i,j));
    end
end
for j = 1:ny-1
    ori_kseedx(nx,j) = fix(width-(width-(nx-1)*step)/2);
    ori_kseedy(nx,j) = (j-1)*step+fix(step/2);
    ori_kseedl(nx,j) = imgL(ori_kseedx(nx,j),ori_kseedy(nx,j));
    ori_kseeda(nx,j) = imga(ori_kseedx(nx,j),ori_kseedy(nx,j));
    ori_kseedb(nx,j) = imgb(ori_kseedx(nx,j),ori_kseedy(nx,j));
end
for i = 1:nx-1
    ori_kseedx(i,ny) = (i-1)*step+fix(step/2);
    ori_kseedy(i,ny) = fix(height-(height-(ny-1)*step)/2);
    ori_kseedl(i,ny) = imgL(ori_kseedx(i,ny),ori_kseedy(i,ny));
    ori_kseeda(i,ny) = imga(ori_kseedx(i,ny),ori_kseedy(i,ny));
    ori_kseedb(i,ny) = imgb(ori_kseedx(i,ny),ori_kseedy(i,ny));
end
ori_kseedx(nx,ny) = fix(width-(width-(nx-1)*step)/2);
ori_kseedy(nx,ny) = fix(height-(height-(ny-1)*step)/2);
ori_kseedl(nx,ny) = imgL(ori_kseedx(nx,ny),ori_kseedy(nx,ny));
ori_kseeda(nx,ny) = imga(ori_kseedx(nx,ny),ori_kseedy(nx,ny));
ori_kseedb(nx,ny) = imgb(ori_kseedx(nx,ny),ori_kseedy(nx,ny));

        