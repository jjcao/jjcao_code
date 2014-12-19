function B = ReadDAT(image_size,data_path)

%读图像superpixel的DAT文件,
%DAT格式:每四个字节存储一个整数
%image_size  图像大小矩阵 [m n] 
%data_path   DAT文件 

row = image_size(1);
colomn = image_size(2);
fid = fopen(data_path,'r');
A = fread(fid, row * colomn, 'uint32');
A = A + 1;
B = reshape(A,[colomn, row]);
B = B';
fclose(fid);