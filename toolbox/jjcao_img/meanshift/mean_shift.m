function output=mean_shift(input_image,Hs,Hr,Th)
%Author  - Suraj Vantigodi, Video Analytics Lab, IISc Bangalore
%inputs- 
%Input_image- an RGB image
%Hs         - spatial range to consider while computing mode
%Hr         - RGB range 
%Th         - Threshold for convergence
%Output-
%meanshift segmented or clustered image
%
% demo 
% x = imread('test.png');output=mean_shift(x,40,3,3);

input1=input_image;
if (size(input1, 3) ~= 3)
    error('rgbhist:numberOfSamples', 'Input image must be RGB.')
end

input1=imresize(input1,[256,256]);
input=input1;


spatial_bandwidth=Hs;
color_bandwidth=Hr;
convergence_fact=Th;
%%%%%%%%%%%%%%%%% Color Histogram %%%%%%%%%%%%%%
fprintf('\n Computing Color histogram');
I=input;
nbins=10;

H1=zeros([nbins nbins nbins]);
ct=0;
for i=1:size(I,1)
    for j=1:size(I,2)
        
        p=double(reshape(I(i,j,:),[1 3]));
        p=floor(p/(256/nbins))+1;
        ct=ct+1;
        arr1(ct,:)=[p(1),p(2),p(3)];
        H1(p(1),p(2),p(3))=H1(p(1),p(2),p(3))+1;            
        pixels1(ct,:)=[I(i,j,:)];    
    end
end

H1=H1(:);
H1=H1./sum(H1);
H_all=reshape(H1,[nbins,nbins,nbins]);

%%%%%%%%%%%%%%%%%%% Histogram done %%%%%%%%%%%%%%%%%%%%

ht=size(input,1);
wd=size(input,2);
output=input;
fprintf('\n Starting meanshift');
tic
for i=1:size(input,1)
    for j=1:size(input,2)
    
    p_y=i;
    p_x=j;
       
    r1 = p_x - spatial_bandwidth; r2 = p_x + spatial_bandwidth;
    c1 = p_y - spatial_bandwidth; c2 = p_y + spatial_bandwidth;
    
    % check boundaries of the region
    if (r1<1) r1 = 1 ; end
    if (c1<1) c1 = 1 ; end
    if (r2>ht) r2 = ht ; end
    if (c2>wd) c2 = wd ; end 
    
    
    patch=input(r1:r2,c1:c2,:); 
    
    R=input(p_y,p_x,1); G=input(p_y,p_x,2); B=input(p_y,p_x,3);
    factor=256/nbins;
    bin_r=ceil(double(R)/factor); bin_g=ceil(double(G)/factor); bin_b=ceil(double(B)/factor);
    old_bins=[bin_r bin_g bin_b];    
           
  %%%%%%%%% Meanshift part %%%%%%%%%%%
    
    dist=convergence_fact+1;
    
    while(dist>convergence_fact)
        
        hr=min(nbins,(bin_r+color_bandwidth)); lr=max(1,(bin_r-color_bandwidth));
        hg=min(nbins,(bin_g+color_bandwidth)); lg=max(1,(bin_g-color_bandwidth));
        hb=min(nbins,(bin_b+color_bandwidth)); lb=max(1,(bin_b-color_bandwidth));
        s_r=0; s_b=0; s_g=0;
        weight=0;
    
    for k=lr:hr
        for l=lg:hg
            for m=lb:hb
                
                s_r=s_r+k*H_all(k,l,m);
                s_g=s_g+l*H_all(k,l,m);
                s_b=s_b+m*H_all(k,l,m);
                weight=weight+H_all(k,l,m);
              
            end
        end
    end
        s_r=s_r/weight; s_g=s_g/weight; s_b=s_b/weight;
        rd=(s_r-bin_r); gd=(s_g-bin_g); bd=(s_b-bin_b);
        dist=sqrt(rd^2+gd^2+bd^2);
        bin_r=round(s_r); bin_g=round(s_g); bin_b=round(s_b);
%       dist
    end
   %%%%%%%%%%%%%%%%%%% meanshift done %%%%%%%%%%%%%%%
   
    color_r=bin_r*(256/nbins); %% computing color to be assigned
    color_g=bin_g*(256/nbins);
    color_b=bin_b*(256/nbins);
   
    output(i,j,1)=color_r;
    output(i,j,2)=color_g;
    output(i,j,3)=color_b;    

    end    
end
fprintf('\n time taken for meanshift=%f',toc);

%%%%%%%%%%%%%%% computing joint histogram of output %%%%%%%%
fprintf('\n Computing Color histogram and plotting\n');
I=output;
ct=0;
for i=1:size(I,1)
    for j=1:size(I,2)
        
        p=double(reshape(I(i,j,:),[1 3]));
        p=floor(p/(256/nbins))+1;
        ct=ct+1;
        arr2(ct,:)=[p(1),p(2),p(3)];              
        pixels2(ct,:)=I(i,j,:);
    end
end

figure(1),subplot(2,2,1),imshow(input); title('input image');
subplot(2,2,2),imshow(output); title('meanshift segmented image')
subplot(2,2,3),plot3(pixels1(:,1),pixels1(:,2),pixels1(:,3),'o');title('color distribution of input image');
subplot(2,2,4),plot3(pixels2(:,1),pixels2(:,2),pixels2(:,3),'o');title('color distribution of output');



