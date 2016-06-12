%haar.m

function image_out=haar(image_in,compression_ratio)

image_haar=haar_comp(image_in);

ratio=ratio_e(image_haar, compression_ratio);

image_haar_lossy=lossy_comp(image_haar,ratio);

image_out=haar_decomp(image_haar_lossy);