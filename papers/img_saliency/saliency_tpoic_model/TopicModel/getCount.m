function count = getCount(labels)
[m,n] = size(labels);
num = max(max(labels));
count = zeros(1,num);
for i = 1:m
    for j = 1:n
       count(labels(i,j)) = count(labels(i,j))+1;
    end
end
end