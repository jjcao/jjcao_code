function ratio_e=ratio_e(image_in,r)

index=1;
while(comp_ratio(image_in,index)<r)
index=index+1;
end

index=index-1;

while(comp_ratio(image_in,index)<r)
index=index+0.1;
end

index=index-0.1;

while(comp_ratio(image_in,index)<r)
index=index+0.01;
end

ratio_e=index;
