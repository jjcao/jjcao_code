function image_out=haar_decomp(image_haar)

h=[1/8  1/8  1/4    0  1/2    0    0    0; 
   1/8  1/8  1/4    0 -1/2    0    0    0;
   1/8  1/8 -1/4    0    0  1/2    0    0; 
   1/8  1/8 -1/4    0    0 -1/2    0    0; 
   1/8 -1/8    0  1/4    0    0  1/2    0; 
   1/8 -1/8    0  1/4    0    0 -1/2    0;
   1/8 -1/8    0 -1/4    0    0    0  1/2; 
   1/8 -1/8    0 -1/4    0    0    0 -1/2];

invhaar = @(x) inv(h)'*x*inv(h);
image_size=size(image_haar);
rows=image_size(1,1);
columns=image_size(1,2);
image_haar_cell=mat2cell( double(image_haar) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell=cellfun( invhaar, image_haar_cell , 'UniformOutput',false) ;
image_out= uint8(cell2mat(image_out_cell)) ;
