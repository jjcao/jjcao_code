function image_haar=haar_comp(image_in)

h=[1/8  1/8  1/4    0  1/2    0    0    0; 
   1/8  1/8  1/4    0 -1/2    0    0    0;
   1/8  1/8 -1/4    0    0  1/2    0    0; 
   1/8  1/8 -1/4    0    0 -1/2    0    0; 
   1/8 -1/8    0  1/4    0    0  1/2    0; 
   1/8 -1/8    0  1/4    0    0 -1/2    0;
   1/8 -1/8    0 -1/4    0    0    0  1/2; 
   1/8 -1/8    0 -1/4    0    0    0 -1/2];

haar = @(x) h'*x*h;
image_size=size(image_in);
rows=image_size(1,1);
columns=image_size(1,2);
image_in_cell=mat2cell( double(image_in) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell=cellfun( haar, image_in_cell , 'UniformOutput',false) ;
image_haar= cell2mat(image_haar_cell) ;
