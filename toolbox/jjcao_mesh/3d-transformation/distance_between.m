function dist = distance_between(v1, v2)
dist = 0;
for i = 1:size(v1,1)
    tmp = repmat(v1(i,:), size(v1,1),1);
    dist = dist + min( sqrt( sum((tmp-v2).^2,2) ));
end