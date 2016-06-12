function lossy_compression_image=lossy_comp(image_in,e)

%image_in(image_in < abs(e) & image_in > -abs(e))=0;
image_in(image_in < abs(e))=0;

lossy_compression_image=image_in;

