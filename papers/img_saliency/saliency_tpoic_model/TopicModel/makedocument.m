function X = makedocument(labels, sp, R, G, B,count)
size = numel(labels);
labels = reshape(labels,1,numel(labels));
sp = reshape(sp,1,numel(sp));
X = zeros(numel(R),max(max(labels)));   
All = zeros(1,numel(R));
for i = 1:numel(R)
    [rnum] = find(sp(1:size)==R(i));
    [gnum] = find(sp((size+1):size*2)==G(i));
    [bnum] = find(sp((size*2+1):size*3)==B(i));
    [alla] = intersect(rnum,gnum);
    [all] = intersect(alla,bnum);
    All(i) = numel(all);
end



for j = 1:max(max(labels))
        if(count(j)==0)
            break;
        end
        [pos] = find(labels==j) ;
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

  %             [rnum] = find(sp(pos)==R(i));
%             [gnum] = find(sp(pos+numel(labels))==G(i));
%             [bnum] = find(sp(pos+numel(labels)*2)==B(i));
%              tmp = isempty(rnum)+isempty(gnum)+isempty(bnum);
%              if(tmp>0)
%                  X(i,j) = 0;
%             else        
