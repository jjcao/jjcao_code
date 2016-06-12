im = imread('buttress.jpg');
row = 1; col = 4;
subplot(row,col,1); imshow(im);
imo=haar(im, 10); subplot(row,col,2); imshow(imo);
imo=haar(im, 50); subplot(row,col,3); imshow(imo);
imo=haar(im, 100); subplot(row,col,4); imshow(imo);