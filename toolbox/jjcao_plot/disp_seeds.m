function h = saliency_sp2im(pos_label_inds_choosen, sp_inds, m, n, w)
 
 tmp_im_saliency=zeros(m,n);
 for i=1:length(pos_label_inds_choosen)
    tmp_im_saliency(sp_inds{pos_label_inds_choosen(i)})=1;
 end

 im_saliency=zeros(w(1),w(2));
 im_saliency(w(3):w(4),w(5):w(6))=tmp_im_saliency;
 im_saliency=uint8(im_saliency*255);  
 h = imshow(im_saliency);
