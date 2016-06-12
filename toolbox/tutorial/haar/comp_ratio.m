function compression_ratio=comp_ratio(image_in,e)


sz=size(find(image_in ~= 0 ),1);

compression_ratio=sz/(sz-size(find(image_in < abs(e) & image_in ~=0),1));

