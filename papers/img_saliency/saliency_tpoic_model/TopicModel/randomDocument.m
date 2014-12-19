function X = randomDocument(documentNum,sp, R, G, B,row, col,depth)

img = reshape(sp,row,col,depth);
X = zeros(numel(R),documentNum);
size = row*col;
%% each color num in the whole pic
All = zeros(1,numel(R));
for i = 1:numel(R)
    [rnum] = find(sp(1:size)==R(i));
    [gnum] = find(sp((size+1):size*2)==G(i));
    [bnum] = find(sp((size*2+1):size*3)==B(i));
    [alla] = intersect(rnum,gnum);
    [all] = intersect(alla,bnum);
    All(i) = numel(all);
end

for j = 1:documentNum

%% find seed patch
seedX = round(rand(1)*(row-1))+1;
seedY = round(rand(1)*(col-1))+1;;
seedrow = round(rand(1)*(min(seedX,row-seedX)-1))+1;
seedcol =  round(rand(1)*(min(seedY,col-seedY)-1))+1;
seedlabels = zeros(row,col);
seedlabels( (seedX-fix(seedrow/2)): (seedX+fix(seedrow/2)), (seedY-fix(seedcol/2)):(seedY+fix(seedcol/2))) = 1;
seedlabels = reshape(seedlabels,1,row*col);
[pos] = find(seedlabels==1);
      for i = 1:numel(R)
            num = 0;
            for k = 1:numel(pos)
                if( sp(pos(k) )==R(i)&&sp(pos(k)+size)==G(i)&&sp(pos(k)+size*2)==B(i))
                    num = num+1;
                end
            end
            X(i,j) = All(i)-num;
      end
end
         



        