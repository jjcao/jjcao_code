function  I = saliency(pzw,img,row,col,R,G,B,topic)

I = zeros(row,col);
for i = 1:numel(R)
    for j = 1:row*col
        if(img(j)==R(i)&&img(j+row*col)==G(i)&&img(j+row*col*2)==B(i));
            I(j) = pzw(topic,i);
        end
    end
end