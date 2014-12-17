function [im_saliency] = saliency_sp2im(fsal, sp_inds, sp_num, m, n, w)
% assign the saliency value to each pixel     
     tmp_im_saliency=zeros(m,n);
     for i=1:sp_num
        tmp_im_saliency(sp_inds{i})=fsal(i);
     end
     tmp_im_saliency=(tmp_im_saliency-min(tmp_im_saliency(:)))/(max(tmp_im_saliency(:))-min(tmp_im_saliency(:)));
     %tmp_im_saliency = normalize(tmp_im_saliency);
     
     im_saliency=zeros(w(1),w(2));
     im_saliency(w(3):w(4),w(5):w(6))=tmp_im_saliency;
     im_saliency=uint8(im_saliency*255);  
     
